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

#include "affine_hull.hpp"
#include "common.hpp"
#include "critical_point.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"
#include "nn_along_ray.hpp"
#include "update_ray.hpp"
#include "vertex_filter.hpp"
#include "logger.hpp"
#include "utility.hpp"

namespace FC {

template <typename point_cloud_t>
class descend_task;  // forward declaration

template <typename point_cloud_t>
class ascend_task {
  using self_t = ascend_task<point_cloud_t>;
public:
  typedef point_cloud_t                                 point_cloud_type;
  typedef typename point_cloud_type::number_type        number_type;
  typedef typename point_cloud_type::size_type          size_type;
  typedef Eigen::Matrix<number_type, Eigen::Dynamic, 1> eigen_vector;
  typedef flow_complex<number_type, size_type>          fc_type;

  /**
    @brief use this constructor, to spawn an ascend_task at a random position
  */
  ascend_task(point_cloud_type const& pc)
    : _ah(pc), _location(pc.dim()), _ray(pc.dim()), _dropped(false, 0) {
    Logger() << "***AT-CTOR: " << this << std::endl;
    std::tuple<size_type, number_type, bool> nn;
    // reseed as long as the nearest neighbor is not unique
    do {
      gen_convex_comb(pc, _location);
      nn = pc.nearest_neighbor(_location);
    } while (std::get<2>(nn));
    // add the nearest neighbor
    _ah.append_point(std::get<0>(nn));
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
  : _ah(std::move(ah)), _location(), _ray(), _dropped(true, dropped_idx) {
    Logger() << "***AT-CTOR: " << this << std::endl;
    assert(_ah.size() == _ah.pc().dim());
    _location.swap(location);
    _ray.swap(ray);
  }
  
    /**
    @brief spawn an ascend task that upflows from a d-1 Delaunay facet
  */
  ascend_task(affine_hull<point_cloud_type> ah,
              eigen_vector && location,
              eigen_vector && ray)
  : ascend_task(std::move(ah), std::move(location), std::move(ray), 0) {
    _dropped.first = false;
  }
  
  ascend_task (ascend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _ray(), _dropped(tmp._dropped) {
    Logger() << "***AT-MOVE-CTOR: " << this << std::endl;
    _location.swap(tmp._location);
    _ray.swap(tmp._ray);
  }
  
  ~ascend_task() {
    Logger() << "***AT-DESTRUCT: " << this << std::endl;
  }

  ascend_task & operator=(ascend_task && tmp) {
    Logger() << "***AT-MOVE-ASSIGN: " << this << std::endl;
    if (this != &tmp) {
      _ah = std::move(tmp._ah);
      _location.swap(tmp._location);
      _ray.swap(tmp._ray);
      _dropped = tmp._dropped;
    }
    return *this;
  }
  
  ascend_task (ascend_task const&) = delete;
  ascend_task & operator=(ascend_task const&) = delete;
  
  template <class DTHandler, class ATHandler, class CIHandler>
  void execute(ATHandler & ath, DTHandler & dth, fc_type & fc,
               CIHandler & cih) {
    auto const& pc = _ah.pc();
    thread_local std::vector<size_type> nnvec(pc.dim() + 1);
    thread_local std::vector<size_type> idx_store(pc.size());
    thread_local eigen_vector driver(pc.dim());
    thread_local eigen_vector lambda(pc.dim() + 1);
    // there's only 2 cases of ascend tasks: completely new, and those starting
    // with d points on the boundary. The first case has not dropped yet, the
    // 2nd case has always dropped before
    assert(_ah.size() == 1 or _ah.size() == pc.dim());
    auto vf = make_at_filter(_ah, _ray,
                             (_dropped.first ? _dropped.second : *_ah.begin()),
                             idx_store.begin());
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
        std::cerr << "ASCEND-TASK error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (nn.first == nnvec.begin()) {  // no nn found -> proxy at inf
        assert(_ah.size() == pc.dim());
        Logger() << "NO STOPPER FOUND - SPAWNING SUB DESCENDS\n";
        auto * inf_ptr = fc.max_at_inf();
        if (_dropped.first) {  // we dropped before flowing to infinity
          Logger() << "DROPPED BEFORE FLOW TO INF\n";
          auto & pos_offsets = nnvec;  // reuse
          assert(pos_offsets.size() >= _ah.size());
          _ah.append_point(_dropped.second);  // append the dropped point
          using ci_type = circumsphere_ident<size_type>;
          if (cih(ci_type(_ah.begin(), _ah.end()))) { // avoid same inf descends
            _ah.project(_location, lambda);
            auto pos_end = get_pos_offsets(lambda, pos_offsets.begin());
            spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                               std::move(_location), std::move(_ah), inf_ptr);
          }
        }  // the else case covers the incidence when we ascend from a d-1 facet
           // in which case we don't need to descend back
        break;  // EXIT 1
      } else {
        Logger() << "STOPPER FOUND\n";
        _location += nn.second * _ray;
        for (auto it = nnvec.begin(); it != nn.first; ++it)
          _ah.append_point(*it);
        // check for finite max
        if (_ah.size() == pc.dim() + 1) {
          simplex_case_upflow(lambda, driver,
                              idx_store.begin(), idx_store.end(),
                              fc, ath, dth);
          break;  // EXIT 2
        } else {
          Logger() << "NOT MAX YET - NEED TO CLIMB\n";
          update_ray<RAY_DIR::FROM_DRIVER>(_ah, _location, lambda,
                                           driver, _ray);
        }
      }
      // make vf consider all ambient points again on subsequent calls
      vf.reset(_ray, _ah, idx_store.begin());
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
    std::uniform_real_distribution<Float> rand_real(0.5, 1.0);
    std::vector<size_type> sampled_indices;
    sampled_indices.reserve(pc.dim() + 1);
    auto already_sampled = [&sampled_indices] (size_type idx) {
      return sampled_indices.end() != std::find(sampled_indices.begin(),
                                                sampled_indices.end(), idx);
    };
    target.setZero();
    Float sum(0.0);
    for (size_type i = 0; i < (pc.dim() + size_type(1)); ++i) {
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
  
  /**
    @brief Handles the case when, during an upflow task, we are at the center of
           of the circumsphere of a full-dimensional-simplex.
  */
  template <class IdxIterator, class DTHandler, class ATHandler>
  void simplex_case_upflow(eigen_vector & lambda,
                           eigen_vector & driver,
                           IdxIterator begin, IdxIterator end,
                           fc_type & fc, ATHandler & ath, DTHandler & dth) {
    auto & ah = _ah;
    auto const& pc = ah.pc();
    auto & x = _location;
    auto & ray = _ray;
    assert(ah.size() == (pc.dim() + 1));
    assert(std::distance(begin, end) >= pc.dim() + 1);
    Logger() << "FINITE MAX SUSPECT\n";
    assert(lambda.size() == ah.size());
    ah.project(x, lambda);
    Logger() << "lambda = " << lambda.head(ah.size()).transpose() << std::endl;
    auto const neg_end = get_neg_offsets(lambda, begin);
    if (begin == neg_end) {  // no points negative
      Logger() << "FINITE MAX INDEED\n"
               << "MAX-LOC " << x.transpose() << std::endl;
      number_type sq_dist((x - pc[*ah.begin()]).squaredNorm());
      using cp_type = typename fc_type::cp_type;
      auto r_pair = fc.insert(cp_type(ah.begin(), ah.end(), sq_dist));
      if (r_pair.first) {  // only spawn descends for new maxima
        auto max_ptr = r_pair.second;
        auto pos_end = std::next(begin, pc.dim() + 1);
        std::iota(begin, pos_end, 0);
        spawn_sub_descends(dth, fc, begin, pos_end, std::move(x), std::move(ah),
                           max_ptr);
      }
    } else {
      Logger() << "NO FINITE MAX - SPAWN NEW ASCEND TASKS\n";
      size_type dropped_idx;
      for (auto it = begin; it != neg_end; ++it) {
        auto new_ah(ah);
        dropped_idx = *(ah.begin() + *it);
        new_ah.drop_point(new_ah.begin() + *it);
        update_ray<RAY_DIR::FROM_DRIVER>(new_ah, x, lambda, driver, ray);
        using ev = eigen_vector;
        using at = self_t;
        Logger() << "NEW TASK: " << new_ah << std::endl
                 << "DROPPED IDX = " << dropped_idx << std::endl;
        ath(at(std::move(new_ah), ev(x), ev(ray), dropped_idx));
      }
    }
  }

  affine_hull<point_cloud_type> _ah;
  eigen_vector                  _location;
  eigen_vector                  _ray;
  std::pair<bool, size_type>    _dropped;
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

