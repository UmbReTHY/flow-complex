#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cassert>

#include <algorithm>
#include <ostream>
#include <iterator>
#include <utility>

#include <Eigen/Core>

#include "affine_hull.hpp"
#include "ascend_task.hpp"
#include "critical_point.hpp"
#include "descend_task.hpp"
#include "logger.hpp"
#include "update_ray.hpp"

namespace FC {

template <class Derived, class Iterator, class Condition>
Iterator get_cond_offsets(Eigen::MatrixBase<Derived> const& lambda,
                          Iterator cond_begin, Condition cond) {
  Iterator cond_end = cond_begin;
  using lambda_t = Eigen::MatrixBase<Derived>;
  for (typename lambda_t::Index i = lambda.size(); i-- > 0;)
    if (cond(lambda[i]))
      *(cond_end++) = i;
  return cond_end;
}

template <class Derived, class Iterator>
Iterator get_pos_offsets(Eigen::MatrixBase<Derived> const& lambda,
                         Iterator pos_begin) {
  using number_t = typename Eigen::MatrixBase<Derived>::Scalar;
  return get_cond_offsets(lambda, pos_begin, [](number_t a) {return a > 0;});
}

template <class Derived, class Iterator>
Iterator get_neg_offsets(Eigen::MatrixBase<Derived> const& lambda,
                         Iterator neg_begin) {
  using number_t = typename Eigen::MatrixBase<Derived>::Scalar;
  return get_cond_offsets(lambda, neg_begin, [](number_t a) {return a < 0;});
}

/**
  @brief Spawns descend tasks for every sub-simplex of ah that is obtained
         by dropping an index from the given range of positive indices.
  @param dth Function object that spawns the given descend_task
  @param drop_pos_begin, drop_pos_end iterators over the position offsets from
                                      ah.begin() that shall be dropped
  @param x the starting position of the yet to spawn descend tasks
  @param ah the simplex to drop from
  @param succ the pointer of the successor in the Hasse-diagram for the new
              descend tasks
*/
template <class DTHandler, class PointCloud, typename number_type,
          typename size_type, class Iterator>
void spawn_sub_descends(DTHandler & dth,
                        Iterator drop_pos_begin, Iterator drop_pos_end,
                        Eigen::Matrix<number_type, Eigen::Dynamic, 1> && x,
                        affine_hull<PointCloud> ah,
                        critical_point<number_type, size_type> * succ) {
  // TODO maybe check before every descend if the cp exists already. If so
  //      the descend would be successful and does not need to be scheduled
  //      Instead, just a (possibly) new successor needs to be added.
  using dt = descend_task<PointCloud>;
  // we skip the first position to use it as the "move-case" after the loop
  assert(drop_pos_begin != drop_pos_end);
  for (Iterator it = std::next(drop_pos_begin); it != drop_pos_end; ++it) {
    auto new_ah(ah);
    new_ah.drop_point(new_ah.begin() + *it);
    auto const& dropped_idx = *(ah.begin() + *it);
    Logger() << "DT takes dropped idx = " << dropped_idx << std::endl;
    using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
    dth(dt(std::move(new_ah), eigen_vector(x), succ, dropped_idx));
  }
  assert(drop_pos_begin != drop_pos_end);
  size_type dropped_idx = *(ah.begin() + *drop_pos_begin);
  Logger() << "DT takes dropped idx = " << dropped_idx << std::endl;
  ah.drop_point(ah.begin() + *drop_pos_begin);
  dth(dt(std::move(ah), std::move(x), succ, dropped_idx));
}

/**
  @brief Handles the case when, during an upflow task, we are at the center of
         of the circumsphere of a full-dimensional-simplex.
*/
template <class PointCloud, typename number_type, class IdxIterator,
          class CPHandler, class DTHandler, class ATHandler>
void simplex_case_upflow(affine_hull<PointCloud> ah,
                         Eigen::Matrix<number_type, Eigen::Dynamic, 1> && x,
                         Eigen::Matrix<number_type, Eigen::Dynamic, 1> & lambda,
                         Eigen::Matrix<number_type, Eigen::Dynamic, 1> && ray,
                         Eigen::Matrix<number_type, Eigen::Dynamic, 1> & driver,
                         IdxIterator begin, IdxIterator end,
                         CPHandler & cph, ATHandler & ath, DTHandler & dth) {
  using size_type = typename PointCloud::size_type;
  auto const& pc = ah.pc();
  assert(ah.size() == (pc.dim() + 1));
  assert(std::distance(begin, end) >= pc.dim() + 1);
  // TODO static_assert: type of IdxIterator == size_type
  Logger() << "FINITE MAX SUSPECT\n";
  assert(lambda.size() == ah.size());
  ah.project(x, lambda);
  Logger() << "lambda = " << lambda.head(ah.size()).transpose() << std::endl;
  auto const neg_end = get_neg_offsets(lambda, begin);
  if (begin == neg_end) {  // no points negative
    Logger() << "FINITE MAX INDEED\n"
             << "MAX-LOC " << x.transpose() << std::endl;
    number_type sq_dist((x - pc[*ah.begin()]).squaredNorm());
    using cp_type = critical_point<number_type, size_type>;
    auto r_pair = cph(cp_type(ah.begin(), ah.end(), sq_dist));
    if (r_pair.first) {  // only spawn descends for new maxima
      auto max_ptr = r_pair.second;
      // we know/assume the point added last is the one in last position
      // within affine hull: so we don't make (begin + dim() + 1) the end
      // of pos coeffs since we don't want to descend back to where we came from
      // TODO because of this, we don't get all the incidences right:
      //      this max_ptr is (possibly) the succ of a cp we found from d-1 cp
      //      or a descend_task of such a d-1 facet
      auto pos_end = std::next(begin, pc.dim() + 1);
      std::iota(begin, pos_end, 0);
      spawn_sub_descends(dth, begin, pos_end, std::move(x), std::move(ah),
                         max_ptr);
    }
  } else {
    Logger() << "NO FINITE MAX - SPAWN NEW ASCEND TASKS\n";
    size_type dropped_idx;  // TODO get rid of this temporary
    // TODO move the last iteration rather than copying it
    for (auto it = begin; it != neg_end; ++it) {
      auto new_ah(ah);
      dropped_idx = *(ah.begin() + *it);
      new_ah.drop_point(new_ah.begin() + *it);
      update_ray<RAY_DIR::FROM_DRIVER>(new_ah, x, lambda, driver, ray);
      using ev = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
      using at = ascend_task<PointCloud>;
      Logger() << "NEW TASK: " << new_ah << std::endl
               << "DROPPED IDX = " << dropped_idx << std::endl;
      ath(at(std::move(new_ah), ev(x), ev(ray), dropped_idx));
    }
  }
}

}  // namespace FC

#endif  // COMMON_HPP_

