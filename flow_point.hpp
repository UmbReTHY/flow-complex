#ifndef FLOW_POINT_HPP_
#define FLOW_POINT_HPP_

#include <cassert>
#include <cstdint>

#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "affine_hull.hpp"
#include "point_cloud.hpp"
#include "utility.hpp"

namespace FC {

template <typename _number_type, typename _size_type>
class flow_point {
  using affine_hull_t = affine_hull<_number_type, _size_type>;

  public:
    typedef _number_type                          number_type;
    typedef _size_type                              size_type;
    typedef point_cloud<number_type, size_type> point_cloud_t;
  
    flow_point(point_cloud_t const& pc)
      : _pc(pc), _location(pc.dim()), _nn_aff_hull(pc),
        _proxy_at_inf_flag(0), _finite_max_flag(0) {
      // 1) generate seed
      std::tuple<_size_type, _number_type, bool> nn;
      // reseed as long as the nearest neighbor is not unique
      do {
        gen_convex_comb(pc.cbegin(), pc.cend(), pc.dim(), _location.data());
        nn = pc.nearest_neighbor(_location.data());
      } while (std::get<2>(nn));
      // add the nearest neighbor
      _nn_aff_hull.add_point(std::get<0>(nn));
    }
    
    flow_point(flow_point const&) = delete;
    flow_point & operator=(flow_point const&) = delete;
    flow_point(flow_point &&) = delete;
    flow_point & operator=(flow_point &&) = delete;
    
    bool is_proxy_at_inf() const {
      return (_proxy_at_inf_flag == 1);
    }
    
    bool is_finite_max() const {
      return (_finite_max_flag == 1);
    }
    
    // TODO: continue here
    void ascend() {
      assert(!is_proxy_at_inf());
      assert(!is_finite_max());
      
      using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
      auto x = Eigen::Map<eigen_vector>(_location.data(), _pc.dim());
      
      // TODO remove
      _proxy_at_inf_flag = 1;
    }

  private:
    point_cloud_t const& _pc;
    std::vector<_number_type> _location;
    affine_hull_t _nn_aff_hull;
    std::uint8_t _proxy_at_inf_flag : 1;
    std::uint8_t   _finite_max_flag : 1;
    std::uint8_t                    : 6;  // unused padding
    // TODO store direct predecessor
};

}  // namespace FC

#endif  // FLOW_POINT_HPP_

