#ifndef ASCEND_TASK_HPP_
#define ASCEND_TASK_HPP_

#include <cassert>

#include <random>
#include <exception>
#include <iterator>
#include <ostream>
#include <tuple>
#include <utility>

#include <Eigen/Core>

#include "common.hpp"
#include "critical_point.hpp"
#include "descend_task.hpp"
#include "affine_hull.hpp"
#include "nn_along_ray.hpp"
#include "update_ray.hpp"
#include "vertex_filter.hpp"
#include "logger.hpp"

namespace FC {

template <typename point_cloud_t>
class descend_task;  // forward declaration

template <typename point_cloud_t>
class ascend_task {
public:
  typedef point_cloud_t                                 point_cloud_type;
  typedef typename point_cloud_type::number_type        number_type;
  typedef typename point_cloud_type::size_type          size_type;
  typedef Eigen::Matrix<number_type, Eigen::Dynamic, 1> eigen_vector;

  /**
    @brief use this constructor, to spawn an ascend_task at a random position
  */
  ascend_task(point_cloud_type const& pc)
    : _ah(pc), _location(pc.dim()), _ray(pc.dim()) {
    Logger() << "***AT-CTOR: " << this << std::endl;
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
  
  /**
    @brief spawn an ascend task that upflows from a d-1 Delaunay facet
  */
  ascend_task(affine_hull<point_cloud_type> ah,
              eigen_vector && location,
              eigen_vector && ray,
              size_type dropped_idx)
  : _ah(std::move(ah)), _location(), _ray(), _dropped_idx(dropped_idx) {
    Logger() << "***AT-CTOR: " << this << std::endl;
    assert(_ah.size() == _ah.pc().dim());
    _location.swap(location);
    _ray.swap(ray);
  }
  
  ascend_task (ascend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _ray(),
      _dropped_idx(tmp._dropped_idx)  {
    Logger() << "***AT-MOVE-CTOR: " << this << std::endl;
    _location.swap(tmp._location);
    _ray.swap(tmp._ray);
  }
  
  ~ascend_task() {
    Logger() << "***AT-DESTRUCT: " << this << std::endl;
  }
  
  ascend_task (ascend_task const&) = delete;
  ascend_task & operator=(ascend_task const&) = delete;
  ascend_task & operator=(ascend_task &&) = delete;
  
  template <class DTHandler, class ATHandler, class CPHandler>
  void execute(DTHandler & dth, ATHandler & ath, CPHandler & cph) {
    auto const& pc = _ah.pc();
    // TODO the vectors can be allocated as thread local storage
    std::vector<size_type> nnvec(pc.dim() + 1);  // TODO dynarray
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    // there's only 2 cases of ascend tasks: completely new, and those starting
    // with d points on the boundary. The first case has not dropped yet, the
    // 2nd case has always dropped before
    assert(_ah.size() == 1 or _ah.size() == pc.dim());
    auto *const dropped_begin = &_dropped_idx;
    auto * dropped_end = (_ah.size() == 1 ? dropped_begin
                                          : std::next(dropped_begin));
    // TODO vertex filter in descend task takes a single index at max too
    //      -> instead of giving a range, give a single pointer instead
    auto vf = make_vertex_filter(_ah, dropped_begin, dropped_end);
    auto nn = std::make_pair(nnvec.begin(), number_type(0));
    
    Logger() << "ASCEND-TASK-STARTS\n";
    do {
      try {
        Logger() << _ah << std::endl
                 << "LOCATION = " << _location.transpose() << std::endl
                 << "DRIVER = " << driver.transpose() << std::endl
                 << "RAY = " << _ray.transpose() << std::endl;
        auto const max_nn = pc.dim() + 1 - _ah.size();
        nn = nearest_neighbor_along_ray(_location, _ray, pc[*_ah.begin()],
                                        vf, nnvec.begin(),
                                        std::next(nnvec.begin(), max_nn));
      } catch(std::exception & e) {
        std::cerr << "error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (nn.first == nnvec.begin()) {  // no nn found -> proxy at inf
        assert(_ah.size() == pc.dim());
        Logger() << "NO STOPPER FOUND - SPAWNING SUB DESCENDS\n";
        // TODO introduce proxy store: don't descend from the same proxy twice
        // TODO find a cheaper way to get the inf-ptr
        using cp_type = critical_point<number_type, size_type>;
        auto * inf_ptr = cph(cp_type(pc.dim())).second;  // addr of cp at inf
        if (dropped_begin == dropped_end) {
          using dt = descend_task<point_cloud_type>;
          size_type const dummy_ignore = *_ah.begin();
          dth(dt(std::move(_ah), std::move(_location), inf_ptr, dummy_ignore));
        } else {  // we dropped before flowing to infinity
          Logger() << "DROPPED BEFORE FLOW TO INF\n";
          auto & pos_offsets = nnvec;  // reuse
          assert(pos_offsets.size() >= _ah.size());
          _ah.add_point(_dropped_idx);  // append the dropped point
          _ah.project(_location, lambda);
          // TODO in case we came from an ascend with multiple negative indices
          //      this ascend task has one or more twins, that might also flow to infinity
          //      and consequently will spawn the exact same descend tasks: duplicate work
          auto pos_end = get_pos_offsets(lambda, pos_offsets.begin());
          spawn_sub_descends(dth, pos_offsets.begin(), pos_end,
                             std::move(_location), std::move(_ah), inf_ptr);
        }
        break;  // EXIT 1
      } else {
        Logger() << "STOPPER FOUND\n";
        _location += nn.second * _ray;
        for (auto it = nnvec.begin(); it != nn.first; ++it)
          _ah.add_point(*it);
        // check for finite max
        if (_ah.size() == pc.dim() + 1) {
          simplex_case_upflow(std::move(_ah), std::move(_location), lambda,
                              std::move(_ray), driver,
                                        nnvec.begin(), nnvec.end(),
                              cph, ath, dth);
          break;  // EXIT 2
        } else {
          Logger() << "NOT MAX YET - NEED TO CLIMB\n";
          update_ray<RAY_DIR::FROM_DRIVER>(_ah, _location, lambda,
                                           driver, _ray);
          dropped_end = dropped_begin;
        }
      }
      assert(dropped_begin == dropped_end);
      vf.reset(dropped_end);  // make the vertex filter consider all points
                              // of the point cloud again on subsequent calls
    } while(true);
    Logger() << "***AT-COMPLETE***\n";
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
    int seed = rd();
    Logger() << "seed = " << seed << std::endl;
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
  size_type                     _dropped_idx;
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

