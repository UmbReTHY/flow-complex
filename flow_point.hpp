#ifndef FLOW_POINT_HPP_
#define FLOW_POINT_HPP_

#include <cassert>
#include <cstdint>
#include <cstdio>

#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "affine_hull.hpp"
#include "critical_point.hpp"
#include "point_cloud.hpp"
#include "utility.hpp"

namespace FC {

template <typename _number_type, typename _size_type>
class flow_point {
  using affine_hull_t = affine_hull<_number_type, _size_type>;
  using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;

  public:
    typedef _number_type                          number_type;
    typedef _size_type                              size_type;
    typedef point_cloud<number_type, size_type> point_cloud_t;
    typedef critical_point<number_type, size_type> cp_t;
  
    flow_point(point_cloud_t const& pc)
      : _nn_aff_hull(pc), _location(pc.dim()), _pc(pc), _succ(nullptr),
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
    
    void ascend() {
      assert(!is_proxy_at_inf());
      assert(!is_finite_max());

      // compute the driver
      eigen_vector d(_pc.dim());
      compute_driver(&d);

      eigen_vector ray = _location - d;
      auto nn = nearest_neighbor_along_ray(_location.data(),
                                           ray.data(),
                                           _pc[*_nn_aff_hull.begin()],
                                           _pc.dim(),
                                           // TODO respect dropped indices and members
                                           [](){});
      if (std::get<2>(nn) != NN_ERROR::OK) {
        switch (std::get<2>(nn)) {
          case NN_ERROR::DIV_BY_ZERO : {
            std::printf("division by zero in 'nearest_neighbor_along_ray'. "
                        "Perturb data.\n");
            break;
          }
          case NN_ERROR::MORE_THAN_ONE_NN : {
            std::printf("more than on nearest neighbor found in "
                        "'nearest_neighbor_along_ray'. Perturb data.\n");
            break;
          }
          default: assert(false && "unknown error in "
                                   "nearest_neighbor_along_ray");
        }
        std::exit(EXIT_FAILURE);
      } else {
        if (std::get<1>(nn) == std::numeric_limits<number_type>::infinity()) {
          _proxy_at_inf_flag = 1;
        } else {
          _location += std::get<1>(nn) * ray;
          _nn_aff_hull.add_point(std::get<0>(nn));
          // check for finite max
          if (_nn_aff_hull.size() == _pc.dim() + 1) {
            std::vector<number_type> lambda(_nn_aff_hull.size());
            _nn_aff_hull.project(_location.data(), lambda.data());
            auto m_it = _nn_aff_hull.begin();
            for (auto coeff : lambda) {
              // TODO remember dropped indices
              if (coeff < 0)
                _nn_aff_hull.drop_point(*m_it);
              m_it++;
            }
            if (_nn_aff_hull.size() == _pc.dim() + 1)
              _finite_max_flag = 1;
          }
        }
      }
    }
    
    /**
      @return Can be nullptr. This indicates that this flow point is not a
              descendent of a maximum.
    */
    cp_t const* get_dir_succ() const noexcept {
      return _succ;
    }
    
    void set_dir_succ(cp_t const* succ) noexcept {
      assert(succ != nullptr);
      _succ = succ;
    }

  private:
    void compute_driver(eigen_vector * d) const {
      std::vector<number_type> lambda(_nn_aff_hull.size());
      _nn_aff_hull.project(_location.data(), lambda.data());
      d->setZero();
      auto m_it = _nn_aff_hull.begin();
      for (auto coeff : lambda)
        d += coeff * Eigen::Map<eigen_vector const>(_pc[*m_it++], _pc.dim());
    }
  
    affine_hull_t _nn_aff_hull;
    eigen_vector _location;
    point_cloud_t const& _pc;
    cp_t const* _succ;
    std::uint8_t _proxy_at_inf_flag : 1;
    std::uint8_t   _finite_max_flag : 1;
    std::uint8_t                    : 6;  // unused padding
    std::uint8_t _unused_bytes[7];        // would be added anyway
};

}  // namespace FC

#endif  // FLOW_POINT_HPP_

