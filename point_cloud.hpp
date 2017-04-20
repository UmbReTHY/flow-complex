#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <cassert>

#include <iterator>
#include <memory>
#include <tuple>
#include <vector>

#include <gflags/gflags.h>
#include <Eigen/Core>
#ifdef _MSC_VER
  #pragma warning(push)
  #pragma warning( disable : 4267)
    #include <nanoflann.hpp>
  #pragma warning(pop)
#else
  #include <nanoflann.hpp>
#endif

#include "utility.hpp"

DEFINE_int32(kdtree_leaf_size, 16, "maximal number of indices in kd-tree "
                                   "leaf nodes");

namespace FC {

template <class PointCloud>
struct NanoflannDataAdaptor {
using map_type = typename PointCloud::eigen_map;
public:
  typedef typename PointCloud::number_type Float;

  NanoflannDataAdaptor(PointCloud const& pc) : _pc(pc) {
  }

  // Must return the number of data points
  std::size_t kdtree_get_point_count() const {
    return _pc.size();
  }
  
  // Must return the Euclidean (L2) distance between the vector "p1[0:size-1]"
  // and the data point with index "idx_p2" stored in the class:
  Float kdtree_distance(Float const* p1, std::size_t const idx_p2,
                        size_t size) const {
    return (map_type(p1, size) - _pc[idx_p2]).squaredNorm();
  }

  // Must return the dim'th component of the idx'th point in the class:
  Float kdtree_get_pt(const size_t idx, int dim) const {
    return _pc[idx][dim];
  }

  // Optional bounding-box computation: return false to default to a
  // standard bbox computation loop.
  template <class BBOX>
  bool kdtree_get_bbox(BBOX &) const {
    return false;
  }
private:
  PointCloud const& _pc;
};

/**
  @brief a wrapper to a vector with additional information about the dimension
         of the point cloud
*/
template <typename _number_type, typename _size_type, bool Aligned>
class point_cloud {
using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
// definitions for nanoflann
using self_t = point_cloud<_number_type, _size_type, Aligned>;
using DataAdaptor = NanoflannDataAdaptor<self_t>;
using DistFunc = nanoflann::L2_Simple_Adaptor<_number_type, DataAdaptor>;
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<DistFunc, DataAdaptor,
                          -1,  // dimension of the dataset: -1 == known at rt
                          _size_type>;
// forward declaration
template <class Iterator> struct RadiusSearchResultSet;
public:
  typedef _number_type                            number_type;
  typedef _size_type                              size_type;
  typedef Eigen::Map<eigen_vector const,
                     eigen_align<Aligned>::value> eigen_map;
private:
  using pt_cont = std::vector<eigen_map>;
public:
  typedef typename pt_cont::const_iterator        iterator;
  /**
    @tparam Iterator when dereferenced, returns a pointer to number_type
  */
  template <typename Iterator, typename dim_type>
  point_cloud(Iterator begin, Iterator end, dim_type dim)
  : _points(), _data_adaptor(*this), _kd_tree(nullptr), diameter_(0.0) {
    auto const size = std::distance(begin, end);
    _points.reserve(size);
    for (auto it = begin; it != end; ++it)
      _points.emplace_back(&((*it)[0]), dim);
    using Params = nanoflann::KDTreeSingleIndexAdaptorParams;
    CHECK(FLAGS_kdtree_leaf_size > 0);
    _kd_tree.reset(new KDTree(dim, _data_adaptor,
                              Params(FLAGS_kdtree_leaf_size)));
    _kd_tree->buildIndex();
    // compute diameter of dataset
    number_type sq_diameter = 0.0;
    #pragma omp parallel
    {
      number_type my_sq_diameter = 0.0;
      #pragma omp for schedule(guided, 100)
      for (int i = 0; i < this->size(); ++i) {
        for (int j = i + 1; j < this->size(); ++j) {
          const auto& u = _points[i];
          const auto& v = _points[j];
          const number_type sq_distance = (u-v).squaredNorm();
          if (sq_distance > my_sq_diameter) my_sq_diameter = sq_distance;
        }
      }
      // workaround: MSVC does not support reduction(max:sq_diameter)
      #pragma omp critical (max_reduction)
      if (my_sq_diameter > sq_diameter) sq_diameter = my_sq_diameter;
    }
    diameter_ = sqrt(sq_diameter);
  }
  
  iterator begin() const noexcept {
    return _points.cbegin();
  }
  
  iterator end() const noexcept {
    return _points.cend();
  }

  template <typename Index>
  eigen_map const& operator[](Index idx) const {
    CHECK(typename pt_cont::size_type(idx) < _points.size());
    return _points[idx];
  }
  
  size_type dim() const noexcept {
    CHECK(!_points.empty() && "YOU CREATED AN EMPTY POINT CLOUD.");
    return convertSafelyTo<size_type>(_points[0].size());
  }
  
  size_type size() const noexcept {
    return convertSafelyTo<size_type>(_points.size());
  }
  
  number_type diameter() const {return diameter_;}
  
  /**
    @brief Finds the nearest neighbor to q. In case of more than one nearest
           neighbor, returns the lowest index.
    @param q query point. pointer to at least dim() elements.
    @return Let r be the returned tuple.
            std::get<0>(r) contains the index of the nearest point to q.
            std::get<1>(r) contains the squared distance of this point to q.
            std::get<2>(r) is set to true if there is more than one nearest
                           neighbor, and false otherwise.
  */
  template <typename Derived>
  std::tuple<size_type, number_type, bool>
  nearest_neighbor(Eigen::MatrixBase<Derived> const& q) const {
    auto r = std::make_tuple(size_type(0),
                             number_type(0),
                             false);
    bool nn_found = false;
    size_type pos = 0;
    for (const auto& p : _points) {
      const number_type tmp = (q - p).squaredNorm();
      if (!nn_found || tmp < std::get<1>(r)) {
        std::get<0>(r) = pos;
        std::get<1>(r) = tmp;
        std::get<2>(r) = false;
        nn_found = true;
      } else if (tmp == std::get<1>(r)) {
        std::get<2>(r) = true;
      }
      ++pos;
    }
    return r;
  }
  
  template <class Derived>
  void k_nearest_neighors(Eigen::MatrixBase<Derived> const& q,
                          size_type k,
                          size_type * indices,
                          number_type * sq_distances) {
    _kd_tree->knnSearch(&q[0], k, indices, sq_distances);
  }
  
  template <class Iterator, class Derived>
  Iterator radius_search(Eigen::MatrixBase<Derived> const& q,
                         number_type squared_radius,
                         Iterator indices_begin) const {
    RadiusSearchResultSet<Iterator> result(indices_begin, squared_radius);
    using SearchParams = nanoflann::SearchParams;
    _kd_tree->findNeighbors(result, &q[0], SearchParams(32, 0, false));
        
    return result.end();
  }

private:
  template <class Iterator>
  struct RadiusSearchResultSet {
  public:
    RadiusSearchResultSet(Iterator begin, number_type sq_radius)
    : _begin(begin), _sq_radius(sq_radius) {
    }
  
    void addPoint(number_type dist, size_type index) {
      assert(dist <= worstDist());
      *(_begin++) = index;
    }
    
    number_type worstDist () const {
      return _sq_radius;
    }
    
    Iterator end() const {
      return _begin;
    }
    
    bool full() const {
			return false;
		}
    
  private:
    Iterator        _begin;
    number_type _sq_radius;
  };

  pt_cont                       _points;
  DataAdaptor             _data_adaptor;
  std::unique_ptr<KDTree>      _kd_tree;
  _number_type                diameter_;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

