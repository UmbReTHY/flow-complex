#ifndef ASCEND_TASK_HPP_
#define ASCEND_TASK_HPP_

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <utility>
#include <random>
#include <tuple>
#include <limits>

#include <Eigen/Core>

#include "descend_task.hpp"
#include "point_cloud.hpp"
#include "affine_hull.hpp"
#include "nn_along_ray.hpp"

namespace FC {

template <typename _point_cloud_type>
class ascend_task {
  using number_type = typename _point_cloud_type::number_type;
  using size_type = typename _point_cloud_type::size_type;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;

  public:
    typedef _point_cloud_type              point_cloud_type;
    typedef descend_task<point_cloud_type> descend_task_type;
    
    ascend_task(point_cloud_type const& pc)
      : _nn_aff_hull(pc), _location(pc.dim()), _ray(pc.dim()) {
      // generate seed
      std::tuple<size_type, number_type, bool> nn;
      // reseed as long as the nearest neighbor is not unique
      do {
        gen_convex_comb(pc, _location);
        nn = pc.nearest_neighbor(_location);
      } while (std::get<2>(nn));
      // add the nearest neighbor
      _nn_aff_hull.add_point(std::get<0>(nn));
      // set the ray
      _ray = pc[*_nn_aff_hull.begin()] - _location;
    }
    
    ascend_task(descend_task_type && dt)
    : _nn_aff_hull(std::move(dt._nn_aff_hull)) {
      _location.swap(dt._location);
      _ray.swap(dt._ray);
    }
    
    /**
      @return true, if finite maximum, false if maximum at infinity
    */
    bool execute() {
      // TODO a DEBUG-only member, to track double calls
      auto const& pc = _nn_aff_hull.pc();
      // TODO dropped point storage does not have to be a member!
      do {
        // TODO respect dropped indices
        size_type curr_idx = 0;  // for lambda below
        auto get_next = [&](std::size_t * idx_ptr) {
          while (_nn_aff_hull.end() != std::find(_nn_aff_hull.begin(),
                                                 _nn_aff_hull.end(), curr_idx))
            ++curr_idx;  // skip affine hull members
          if (curr_idx < pc.size()) { // valid index
            assert(idx_ptr);
            *idx_ptr = curr_idx;
            return &pc[curr_idx++];  // this ++ is important!!!
          } else {
            return static_cast<decltype(&pc[curr_idx])>(nullptr);
          }
        };
        auto nn = nearest_neighbor_along_ray(_location, _ray,
                                             pc[*_nn_aff_hull.begin()],
                                             get_next);
        if (std::get<2>(nn) != NN_ERROR::OK) {
          switch (std::get<2>(nn)) {
            case NN_ERROR::DIV_BY_ZERO : {
              std::printf("division by zero in 'nearest_neighbor_along_ray'. "
                          "Perturb data.\n");
              break;
            }
            case NN_ERROR::MORE_THAN_ONE_NN : {
              std::printf("more than one nearest neighbor found in "
                          "'nearest_neighbor_along_ray'. Perturb data.\n");
              break;
            }
            default: assert(false && "unknown error in "
                                     "nearest_neighbor_along_ray");
          }
          std::exit(EXIT_FAILURE);
        } else {
          if (std::get<1>(nn) == std::numeric_limits<number_type>::infinity()) {
            break;  // EXIT 1
          } else {
            _location += std::get<1>(nn) * _ray;
            _nn_aff_hull.add_point(std::get<0>(nn));
            // check for finite max
            if (_nn_aff_hull.size() == pc.dim() + 1) {
              eigen_vector lambda(_nn_aff_hull.size());
              _nn_aff_hull.project(_location, lambda);
              // drop negative indices
              auto m_it = _nn_aff_hull.begin();
              for (size_type i = 0; i < lambda.size(); ++i) {
                // TODO remember dropped indices
                if (lambda[i] < 0)
                  _nn_aff_hull.drop_point(*m_it);
                m_it++;
              }
              if (_nn_aff_hull.size() == pc.dim() + 1)
                break;    // EXIT 2 - none have been dropped -> finite max
              else
                update_ray();
            } else {
              update_ray();
            }
          }
        }
      } while(true);
      assert(_nn_aff_hull.size() == pc.dim() or
             _nn_aff_hull.size() == pc.dim() + 1);
      return _nn_aff_hull.size() == pc.dim();
    }
    
  private:
    void update_ray() {
      // TODO measure performance for lambda/driver being members: no reallocation
      eigen_vector lambda(_nn_aff_hull.size());
      _nn_aff_hull.project(_location, lambda);
      eigen_vector driver = eigen_vector::Zero(_location.size());
      auto m_it = _nn_aff_hull.begin();
      for (size_type i = 0; i < _nn_aff_hull.size(); ++i)
        driver += lambda[i] * _nn_aff_hull.pc()[*m_it++];
      _ray = _location - driver;
    }
  
    /**
      @brief generates a convex combination of a set of points, where the
             coefficients of the convex combination are strictly greater than 0
      @param begin, end iterators to pointers that hold the coeeficients of the
             points
      @param dim the the size of the point arrays
      @param target a pointer to at least dim elements where the result is placed
             
    */
    void gen_convex_comb(point_cloud_type const& pc, eigen_vector & target) {
      std::random_device rd;
      std::mt19937 gen(rd());
      // it's important for the lower bound not to be 0, otherwise some points
      // might not take part in the convex combination and the seed point would be
      // located on the convex hull, thus getting stuck during subsequent ascends
      // and not reaching a maximum
      std::uniform_real_distribution<number_type> dis(0.5, 1.0);
      number_type sum = 0.0;
      target.setZero();
      number_type tmp;
      for (auto const& p : pc) {
        tmp = dis(gen);
        sum += tmp;
        target = tmp * p;
      }
      target /= sum;
    }
  
    affine_hull<_point_cloud_type> _nn_aff_hull;
    eigen_vector                   _location;
    eigen_vector                   _ray;
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

