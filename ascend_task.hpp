#ifndef ASCEND_TASK_HPP_
#define ASCEND_TASK_HPP_

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <utility>
#include <random>
#include <limits>
#include <exception>

#include <Eigen/Core>

#include "critical_point.hpp"
#include "descend_task.hpp"
#include "point_cloud.hpp"
#include "affine_hull.hpp"
#include "nn_along_ray.hpp"
#include "update_ray.hpp"

namespace FC {

template <typename point_cloud_t>
class ascend_task {
public:
  typedef point_cloud_t                          point_cloud_type;
  typedef typename point_cloud_type::number_type number_type;
  typedef typename point_cloud_type::size_type   size_type;
private:
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using cp_type = critical_point<number_type, size_type>;
public:
  ascend_task(point_cloud_type const& pc)
    : _ah(pc), _location(pc.dim()), _ray(pc.dim()) {
    // generate seed
    std::tuple<size_type, number_type, bool> nn;
    // reseed as long as the nearest neighbor is not unique
    do {
      gen_convex_comb(pc, _location);
      nn = pc.nearest_neighbor(_location);
    } while (std::get<2>(nn));
    // add the nearest neighbor
    _ah.add_point(std::get<0>(nn));
    // set the ray
    _ray = pc[*_ah.begin()] - _location;
  }
  
  ascend_task(affine_hull<point_cloud_type> && ah,
              Eigen::Matrix<number_type, Eigen::Dynamic, 1> && location,
              Eigen::Matrix<number_type, Eigen::Dynamic, 1> && ray)
  : _ah(std::move(ah)) {
    assert(_ah.size() == _ah.pc().dim());  // should be at a (d-1) cp
    _location.swap(location);
    _ray.swap(ray);
  }
  
  /**
    @return true, if finite maximum, false if maximum at infinity
  */
  template <typename ATHandler, typename DTHandler, typename CPHandler>
  void execute(ATHandler & ath, DTHandler & dth, CPHandler & cph) {
    // TODO a DEBUG-only member, to track double calls
    auto const& pc = _ah.pc();
    std::vector<size_type> nnvec(pc.dim());  // for at most d additional nn
    // TODO dropped point storage does not have to be a member!
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    do {
      // TODO respect dropped indices
      size_type curr_idx = 0;  // for lambda below
      auto get_next = [&](std::size_t * idx_ptr) {
        while (_ah.end() != std::find(_ah.begin(),
                                               _ah.end(), curr_idx))
          ++curr_idx;  // skip affine hull members
        if (curr_idx < pc.size()) { // valid index
          assert(idx_ptr);
          *idx_ptr = curr_idx;
          return &pc[curr_idx++];  // this ++ is important!!!
        } else {
          return static_cast<decltype(&pc[curr_idx])>(nullptr);
        }
      };
      auto nn = std::make_pair(nnvec.begin(), number_type(0));
      try {
        nn = nearest_neighbor_along_ray(_location, _ray, pc[*_ah.begin()],
                                        get_next, nnvec.begin(), nnvec.begin() +
                                        (pc.dim() + 1 - _ah.size()));
      } catch(std::exception & e) {
        std::printf("error: %s\n", e.what());
        std::exit(EXIT_FAILURE);
      }
      using dt = descend_task<point_cloud_type>;
      if (nn.second == std::numeric_limits<number_type>::infinity()) {
        assert(_ah.size() == pc.dim());
        dth(dt(std::move(_ah), std::move(_location),
               cph(cp_type(pc.dim())).second));
        break;  // EXIT 1
      } else {
        _location += nn.second * _ray;
        for (auto it = nnvec.begin(); it != nn.first; ++it)
          _ah.add_point(*it);
        // check for finite max
        if (_ah.size() == pc.dim() + 1) {
          eigen_vector lambda(_ah.size());
          _ah.project(_location, lambda);
          // drop negative indices
          auto m_it = _ah.begin();
          for (size_type i = 0; i < lambda.size(); ++i) {
            // TODO remember dropped indices
            if (lambda[i] < 0)
              _ah.drop_point(*m_it);
            m_it++;
          }
          if (_ah.size() == pc.dim() + 1) {
            assert(_ah.size() == pc.dim() + 1);
            number_type dist = (_location - pc[*_ah.begin()]).norm();
            auto r_pair = cph(cp_type(_ah.begin(), _ah.end(), dist));
            if (r_pair.first) {
              auto max_ptr = r_pair.second;
              if (max_ptr->index() > 1) { // don't descend to idx-0 cp-s
                // TODO maybe pass the recently dropped indices
                for (auto it = ++_ah.begin(); it != _ah.end(); ++it) {
                  auto new_ah = _ah;
                  new_ah.drop_point(*it);
                  dth(dt(std::move(new_ah), eigen_vector(_location), max_ptr));
                }
                // reuse this task's affine hull for one of the descent tasks
                _ah.drop_point(*_ah.begin());
                dth(dt(std::move(_ah), std::move(_location), max_ptr));
              }
            }
            break;    // EXIT 2 - none have been dropped -> finite max
          } else {
            update_ray(_ah, _location, lambda, driver, _ray);
          }
        } else {
          update_ray(_ah, _location, lambda, driver, _ray);
        }
      }
    } while(true);
  }
  
private:
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

  affine_hull<point_cloud_type> _ah;
  eigen_vector                  _location;
  eigen_vector                  _ray;
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

