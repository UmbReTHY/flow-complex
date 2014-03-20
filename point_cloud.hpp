#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <cassert>

#include <limits>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

namespace FC {

/**
  @brief a wrapper to a vector with additional information about the dimension
         of the point cloud
*/
template <typename _number_type, typename _size_type>
class point_cloud {
  using pt_cont = std::vector<_number_type const*>;

  public:
    typedef _number_type number_type;
    typedef _size_type     size_type;
    typedef typename pt_cont::const_iterator const_iterator;
  
    /**
      @tparam Iterator when dereferenced, returns a (reference of a) pointer
                       to number_type
    */
    template <typename Iterator>
    point_cloud(Iterator begin, Iterator end, size_type dim)
      : _points(begin, end), _dim(dim) {
    }
    
    const_iterator cbegin() const noexcept {
      return _points.cbegin();
    }
    
    const_iterator cend() const noexcept {
      return _points.cend();
    }
  
    number_type const* operator[](size_type idx) const {
      assert(idx < _points.size());
      return _points[idx];
    }
    
    size_type dim() const noexcept {
      return _dim;
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
    // TODO check exception safety
    std::tuple<size_type, number_type, bool>
    nearest_neighbor(number_type const* q) const {
      auto r = std::make_tuple(size_type(0),
                               std::numeric_limits<number_type>::infinity(),
                               false);
      using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
      using cmap = Eigen::Map<eigen_vector const>;
      auto q_map = cmap(q, _dim);
      for (size_type i = 0; i < _points.size(); ++i) {
        auto tmp = (q_map - cmap(_points[i], _dim)).squaredNorm();
        if (tmp < std::get<1>(r)) {
          std::get<1>(r) = tmp;
          std::get<0>(r) = i;
        } else if (tmp == std::get<1>(r)) {
          std::get<2>(r) = true;
        }
      }
      
      return r;
    }
  
  private:
    pt_cont _points;
    size_type _dim;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

