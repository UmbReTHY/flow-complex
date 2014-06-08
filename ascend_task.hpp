#ifndef ASCEND_TASK_HPP_
#define ASCEND_TASK_HPP_

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <utility>
#include <random>
#include <exception>
#include <iostream>

#include <Eigen/Core>

#include "critical_point.hpp"
#include "descend_task.hpp"
#include "affine_hull.hpp"
#include "nn_along_ray.hpp"
#include "update_ray.hpp"
#include "vertex_filter.hpp"
#include "utility.hpp"

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
  
  // TODO reconsider these r-value references with respect to this constructor's usage in dt
  ascend_task(affine_hull<point_cloud_type> ah,  // TODO remove copies here
              Eigen::Matrix<number_type, Eigen::Dynamic, 1> location,
              Eigen::Matrix<number_type, Eigen::Dynamic, 1> ray)
  : _ah(std::move(ah)) {
    assert(_ah.size() == _ah.pc().dim());  // should be at a (d-1) cp
    _location.swap(location);
    _ray.swap(ray);
  }
  
  template <typename DTHandler, typename CPHandler>
  void execute(DTHandler & dth, CPHandler & cph) {
    auto const& pc = _ah.pc();
    std::vector<size_type> nnvec(pc.dim());  // for at most d additional nn
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    auto vf = make_vertex_filter(_ah);
    auto nn = std::make_pair(nnvec.begin(), number_type(0));
    do {
      try {
        nn = nearest_neighbor_along_ray(_location, _ray, pc[*_ah.begin()],
                                        vf, nnvec.begin(), nnvec.begin() +
                                        (pc.dim() + 1 - _ah.size()));
      } catch(std::exception & e) {
        std::cerr << "error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      using dt = descend_task<point_cloud_type>;
      if (nn.first == nnvec.begin()) {  // no nn found -> proxy at inf
        assert(_ah.size() == pc.dim());
        dth(dt(std::move(_ah), std::move(_location),
               cph(cp_type(pc.dim())  // passes a cp at inf to cph
                                      // no insert will happen since there's
                                      // only 1 cp at inf which is already
                                      // part of the fc
                  ).second)  // and returns its address to be used as the
                             // successor of possibly found cps by this dt
           );
        break;  // EXIT 1
      } else {
        _location += nn.second * _ray;
        for (auto it = nnvec.begin(); it != nn.first; ++it) {
          _ah.add_point(*it);
        }
        // check for finite max
        if (_ah.size() == pc.dim() + 1) {
          if (not drop_neg_coeffs(_location, lambda, _ah)) {  // no points dropped
            number_type sq_dist((_location - pc[*_ah.begin()]).squaredNorm());
            auto r_pair = cph(cp_type(_ah.begin(), _ah.end(),
                                      std::move(sq_dist)));
            if (r_pair.first) {  // only spawn descends for new maxima
              auto max_ptr = r_pair.second;
              spawn_sub_descends(dth, _location, _ah.size(), max_ptr, _ah);
              dth(dt(std::move(_ah), std::move(_location), max_ptr));
            }
            break;    // EXIT 2 - none have been dropped -> finite max
          } else {
            update_ray<RAY_DIR::FROM_DRIVER>(_ah, _location, lambda,
                                             driver, _ray);
          }
        } else {
          update_ray<RAY_DIR::FROM_DRIVER>(_ah, _location, lambda,
                                           driver, _ray);
        }
      }
      vf.reset();  // make the vertex filter consider all points
                   // of the point cloud again on subsequent calls
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
    using Float = float;
    std::random_device rd;
    std::mt19937 gen(rd());
    // generates numbers in [0, pc.size() - 1]
    std::uniform_int_distribution<size_type> rand_idx(0, pc.size() - 1);
    std::vector<size_type> sampled_indices;
    sampled_indices.reserve(pc.dim() + 1);
    auto already_sampled = [&sampled_indices] (size_type idx) {
      return sampled_indices.end() != std::find(sampled_indices.begin(),
                                                sampled_indices.end(), idx);
    };
    target.setZero();
    Float const fraction = 1.0 / (pc.dim() + 1);
    for (size_type i = 0; i < pc.dim() + size_type(1); ++i) {
      size_type idx;
      do {
        idx = rand_idx(gen);
        assert(0 <= idx);
        assert(idx < pc.size());
      } while(already_sampled(idx));
      sampled_indices.push_back(idx);
      target += fraction * pc[idx];
    }
  }

  affine_hull<point_cloud_type> _ah;
  eigen_vector                  _location;
  eigen_vector                  _ray;
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

