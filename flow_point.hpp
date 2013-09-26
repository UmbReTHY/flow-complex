#ifndef FLOW_POINT_HPP_
#define FLOW_POINT_HPP_

#include <cassert>

#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "point_cloud.hpp"
#include "affine_hull.hpp"

namespace fc {

template <typename _number_type, typename _size_type, typename _dim_type>
class flow_point {
  using affine_hull_t = affine_hull<_number_type, _size_type, _dim_type>;

  public:
    typedef _number_type   number_type;
    typedef _size_type       size_type;
    typedef point_cloud<number_type, _dim_type, size_type> point_cloud_t;
  
    flow_point(point_cloud_t const& pc, number_type const* loc_ptr)
      : _point_cloud(pc), _location(loc_ptr, loc_ptr + pc.dim()),
        _nn_aff_hull(pc), _is_proxy_at_inf(false), _is_finite_max(false) {
      auto const nns = nearest_neighbors(loc_ptr),
      for (auto const& nn : nns)
        _nn_aff_hull.add_point(nn);
      // degenerate input
      if (_nn_aff_hull.size() > (pc.dim() + 1))
        throw std::invalid_argument("point cloud not in general position");
      // TODO: check for finite max if size = d + 1
    }
    
    bool is_proxy_at_inf() const {
      return _is_proxy_at_inf;
    }
    
    bool is_finite_max() const {
      return _is_finite_max;
    }
    
    // TODO: continue here
    void ascend() {
      assert(!is_proxy_at_inf());
      assert(!is_finite_max());
      
      using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
      auto x = Eigen::Map<eigen_vector>(_location.data());
    }

  private:
    point_cloud_t const& _point_cloud;
    std::vector<_number_type> _location;
    affine_hull_t _nn_aff_hull;
    bool _is_proxy_at_inf;
    bool _is_finite_max;
};

}  // namespace fc

#endif  // FLOW_POINT_HPP_

