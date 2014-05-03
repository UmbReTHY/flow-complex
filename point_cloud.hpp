#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <cassert>

#include <iterator>
#include <tuple>
#include <vector>

#include <Eigen/Core>

#include "utility.hpp"

namespace FC {

/**
  @brief a wrapper to a vector with additional information about the dimension
         of the point cloud
*/
template <typename _number_type, typename _size_type, bool Aligned>
class point_cloud {
using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;

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
    point_cloud(Iterator begin, Iterator end, dim_type dim) {
      auto const size = std::distance(begin, end);
      _points.reserve(size);
      for (auto it = begin; it != end; ++it)
        _points.emplace_back(*it, dim);
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
    
    typename eigen_map::Index dim() const noexcept {
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
  
  private:
    pt_cont _points;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

