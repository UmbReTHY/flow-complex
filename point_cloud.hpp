#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <cassert>

#include <iterator>
#include <memory>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <nanoflann.hpp>

#include "utility.hpp"

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
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<DistFunc, DataAdaptor>;
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
  : _points(), _data_adaptor(*this), _kd_tree(nullptr) {
    auto const size = std::distance(begin, end);
    _points.reserve(size);
    for (auto it = begin; it != end; ++it)
      _points.emplace_back(*it, dim);
    using Params = nanoflann::KDTreeSingleIndexAdaptorParams;
    _kd_tree.reset(new KDTree(dim, _data_adaptor, Params(15)));
    _kd_tree->buildIndex();
  }
  
  iterator begin() const noexcept {
    return _points.cbegin();
  }
  
  iterator end() const noexcept {
    return _points.cend();
  }

  template <typename Index>
  eigen_map const& operator[](Index idx) const {
    assert(idx < _points.size());
    return _points[idx];
  }
  
  size_type dim() const noexcept {
    assert(not _points.empty() && "YOU CREATED AN EMPTY POINT CLOUD.");
    return _points[0].size();
  }
  
  size_type size() const noexcept {
    return _points.size();
  }
  
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
    for (size_type i = 0; i < _points.size(); ++i) {
      number_type tmp = (q - _points[i]).squaredNorm();
      if ((not nn_found) or (tmp < std::get<1>(r))) {
        std::get<1>(r) = tmp;
        std::get<0>(r) = i;
        nn_found = true;
      } else if (tmp == std::get<1>(r)) {
        std::get<2>(r) = true;
      }
    }
    return r;
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
    
  private:
    Iterator        _begin;
    number_type _sq_radius;
  };

  pt_cont                       _points;
  DataAdaptor             _data_adaptor;
  std::unique_ptr<KDTree>      _kd_tree;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

