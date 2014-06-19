#ifndef ASCEND_TASK_HPP_
#define ASCEND_TASK_HPP_

#include <cassert>
#include <cstdio>
#include <cstdint>

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
    std::cout << "***AT-CTOR: " << this << std::endl;
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
    _ray = _location - pc[*_ah.begin()];
  }
  
  // TODO reconsider these r-value references with respect to this constructor's usage in dt
  ascend_task(affine_hull<point_cloud_type> ah,  // TODO remove copies here
              Eigen::Matrix<number_type, Eigen::Dynamic, 1> location,
              Eigen::Matrix<number_type, Eigen::Dynamic, 1> ray)
  : _ah(std::move(ah)), _location(), _ray() {
    std::cout << "***AT-CTOR: " << this << std::endl;
    assert(_ah.size() == _ah.pc().dim());  // should be at a (d-1) cp
    _location.swap(location);
    _ray.swap(ray);
  }
  
  ascend_task (ascend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _ray() {
    std::cout << "***AT-MOVE-CTOR: " << this << std::endl;
    _location.swap(tmp._location);
    _ray.swap(tmp._ray);
  }
  
  ~ascend_task() {
    std::cout << "***AT-DESTRUCT: " << this << std::endl;
  }
  
  ascend_task (ascend_task const&) = delete;
  ascend_task & operator=(ascend_task const&) = delete;
  ascend_task & operator=(ascend_task &&) = delete;
  
  template <typename DTHandler, typename CPHandler>
  void execute(DTHandler & dth, CPHandler & cph) {
    auto const& pc = _ah.pc();
    std::vector<size_type> nnvec(pc.dim());  // for at most d additional nn
    std::vector<size_type> dropvec(pc.dim());  // stores recently dropped pts
    auto dropped_end = dropvec.begin();
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    auto vf = make_vertex_filter(_ah, dropvec.begin(), dropped_end);
    auto nn = std::make_pair(nnvec.begin(), number_type(0));
    static int id = 0;
    std::cout << "ASCEND-TASK-ID: " << ++id << std::endl;
    std::uint32_t hop_count = 0;
    do {
      ++hop_count;
      try {
        std::cout << "HULL MEMBERS: ";
        for (auto it = _ah.begin(); it != _ah.end(); ++it)
          std::cout << *it << ", ";
        std::cout << std::endl;
        std::cout << "LOCATION = " << _location.transpose() << std::endl;
        std::cout << "DRIVER = " << driver.transpose() << std::endl;
        std::cout << "RAY = " << _ray.transpose() << std::endl;
        nn = nearest_neighbor_along_ray(_location, _ray, pc[*_ah.begin()],
                                        vf, nnvec.begin(), nnvec.begin() +
                                        (pc.dim() + 1 - _ah.size()));
      } catch(std::exception & e) {
        std::cerr << "error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (nn.first == nnvec.begin()) {  // no nn found -> proxy at inf
        std::cout << "NO STOPPER FOUND\n";
//        assert(_ah.size() == pc.dim());  // TODO this might not be justified
        if (hop_count > 1) {  // no dt back to idx-(d-1) cp
          std::cout << "SPAWNED DESCEND BACK\n";
          size_type size_before = _ah.size();
          auto * inf_ptr = cph(cp_type(pc.dim())).second;  // addr of cp at inf
          for (auto it = dropvec.begin(); it != dropped_end; ++it)
            _ah.add_point(*it);
          spawn_sub_descends(dth, size_before, std::move(_location),
                             std::move(_ah), inf_ptr);
        }
        break;  // EXIT 1
      } else {
        std::cout << "STOPPER FOUND\n";
        _location += nn.second * _ray;
        for (auto it = nnvec.begin(); it != nn.first; ++it)
          _ah.add_point(*it);
        // check for finite max
        if (_ah.size() == pc.dim() + 1) {
          std::cout << "FINITE MAX SUSPECT\n";
          dropped_end = drop_neg_coeffs(_location, lambda, _ah,
                                        dropvec.begin());
          if (dropvec.begin() == dropped_end) {  // no points dropped
            std::cout << "FINITE MAX INDEED\n";
            std::cout << "MAX-LOC " << _location.transpose() << std::endl;
            number_type sq_dist((_location - pc[*_ah.begin()]).squaredNorm());
            auto r_pair = cph(cp_type(_ah.begin(), _ah.end(),
                                      std::move(sq_dist)));
            if (r_pair.first) {  // only spawn descends for new maxima
              auto max_ptr = r_pair.second;
              spawn_sub_descends(dth, _ah.size(),
                                 std::move(_location), std::move(_ah), max_ptr);
            }
            break;    // EXIT 2 - none have been dropped -> finite max
          } else {
            std::cout << "HAD TO DROP SOME\n";
            update_ray<RAY_DIR::FROM_DRIVER>(_ah, _location, lambda,
                                             driver, _ray);
          }
        } else {
          std::cout << "NOT MAX YET - NEED TO CLIMB\n";
          update_ray<RAY_DIR::FROM_DRIVER>(_ah, _location, lambda,
                                           driver, _ray);
          dropped_end = dropvec.begin();  // starting from dimension 3, this is
                                          // necessary because more than 1 point
                                          // can be dropped at once: after one
                                          // more round there is no recently
                                          // dropped point anymore
        }
      }
      vf.reset(dropped_end);  // make the vertex filter consider all points
                              // of the point cloud again on subsequent calls
    } while(true);
    std::cout << "***AT-COMPLETE***\n";
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
//    int seed = rd();
    int seed = -460063112;
    std::cout << "seed = " << seed << std::endl;
    std::mt19937 gen(seed);
    // generates numbers in [0, pc.size() - 1]
    std::uniform_int_distribution<size_type> rand_idx(0, pc.size() - 1);
    std::uniform_real_distribution<Float> rand_real;
    std::vector<size_type> sampled_indices;
    sampled_indices.reserve(pc.dim() + 1);
    auto already_sampled = [&sampled_indices] (size_type idx) {
      return sampled_indices.end() != std::find(sampled_indices.begin(),
                                                sampled_indices.end(), idx);
    };
    target.setZero();
    Float sum(0.0);
    for (size_type i = 0; i < pc.dim() + size_type(1); ++i) {
      size_type idx;
      do {
        idx = rand_idx(gen);
        assert(0 <= idx);
        assert(idx < pc.size());
      } while(already_sampled(idx));
      Float tmp = rand_real(gen);
      sum += tmp;
      sampled_indices.push_back(idx);
      target += tmp * pc[idx];
    }
    target /= sum;
  }

  affine_hull<point_cloud_type> _ah;
  eigen_vector                  _location;
  eigen_vector                  _ray;
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

