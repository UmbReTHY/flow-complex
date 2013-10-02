#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <cassert>

#include <iterator>
#include <limits>
#include <vector>

#include <Eigen/Dense>

namespace fc {

template <typename _number_type, typename _dim_type, typename _size_type>
class point_cloud {
  public:
    typedef _number_type number_type;
    typedef _dim_type dim_type;
    typedef _size_type size_type;

  private:
    using point_container = std::vector<_number_type const*>;

  public:
    typedef typename point_container::const_iterator const_iterator;

    /**
      @tparam point_iterator the expression (*<point_iterator>)[i] must be valid
                             and the returned values are expected to be located
                             in contiguous memory
      @params begin,end iterators defining the range of points that constitute
                        the point cloud.
    */
    template <typename point_iterator>
    point_cloud(dim_type const dim, point_iterator begin, point_iterator end)
      : _dim(dim) {
      _points.reserve(std::distance(begin, end));
      for (;begin != end; ++begin)
        _points.push_back(&(*begin)[0]);
    }

    point_cloud() = delete;
    point_cloud(point_cloud const&) = delete;
    point_cloud & operator=(point_cloud const&) = delete;
    point_cloud(point_cloud const&&) = delete;
    point_cloud & operator=(point_cloud const&&) = delete;

    const_iterator begin() const {
      return _points.cbegin();
    }

    const_iterator end() const {
      return _points.cend();
    }

    /**
      @param idx idx of the point within the point cloud
    */
    number_type const* operator[](size_type const idx) const noexcept {
      return _points[idx];
    }

    dim_type dim() const noexcept {
      return _dim;
    }

    size_type size() const {
      return _points.size();
    }
    
    /**
      @return the indices of the nearest neighbors to q
    */
    std::vector<size_type> nearest_neighbors(number_type const* q) {
      std::vector<size_type> nn;

      if (_points.size() > 0) {
        number_type min_dist = squared_distance(q, _points[0], _dim);
        nn.push_back(0);
        for (size_type i = 1; i < _points.size(); ++i) {
          number_type const dist = squared_distance(q, _points[i], _dim);
          if (dist < min_dist) {
            nn.resize(1);
            nn.front() = i;
            min_dist = dist;
          } else if (dist == min_dist) {
            nn.push_back(i);
          }
        }
      }
      
      return std::move(nn);
    }

  private:
    number_type squared_distance(number_type const* a, number_type const* b,
                                 dim_type dim) {
      using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
      return (Eigen::Map<const eigen_vector>(a, dim) -
              Eigen::Map<const eigen_vector>(b, dim)).squaredNorm();
    }

    dim_type const _dim;
    point_container _points;
};

}  // namespace fc

#endif  // POINT_CLOUD_HPP_

