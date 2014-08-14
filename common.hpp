#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cassert>

#include <algorithm>
#include <ostream>
#include <iterator>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "affine_hull.hpp"
#include "ascend_task.hpp"
#include "critical_point.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"
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
          typename size_type, class Iterator, class Iterator2 = int *>
void spawn_sub_descends(DTHandler & dth,
                        flow_complex<number_type, size_type> const& fc,
                        Iterator drop_pos_begin, Iterator drop_pos_end,
                        Eigen::Matrix<number_type, Eigen::Dynamic, 1> && x,
                        affine_hull<PointCloud> ah,
                        critical_point<number_type, size_type> * succ,
                        Iterator2 ignore_begin = nullptr,
                        Iterator2 ignore_end = nullptr) {
  using dt_type = descend_task<PointCloud>;
  using fc_type = flow_complex<number_type, size_type>;
  using cp_type = typename fc_type::cp_type;
  // we skip the first position to use it as the "move-case" after the loop
  assert(drop_pos_begin != drop_pos_end);
  cp_type * existing_cp = nullptr;
  // stores the indices for cp creation - this avoids creating a whole
  // affine hull, in case existing_cp != nullptr
  auto const& pc = ah.pc();
  thread_local std::vector<size_type> newidx(pc.dim() + 1);
  auto const physical_end = std::copy(ah.begin(), ah.end(), newidx.begin());
  auto const logical_end = physical_end - 1;
  auto swap_dropped = [logical_end] (size_type d_idx) {
    auto it = newidx.begin();
    while (true) {
      if (*it == d_idx) {
        std::swap(*it, *logical_end);
        break;
      }
      ++it;
    }
  };
  for (Iterator it = std::next(drop_pos_begin); it != drop_pos_end; ++it) {
    auto const dropped_idx = *(ah.begin() + *it);
    swap_dropped(dropped_idx);
    if ((existing_cp = fc.find(cp_type(newidx.begin(), logical_end, 0)))) {
      Logger() << "CP ALREADY FOUND - NO DT SPAWNED,SUCCS UPDATED\n";
      existing_cp->add_successor(succ);
    } else {
      auto new_ah(ah);
      new_ah.drop_point(new_ah.begin() + *it);
      Logger() << "DT takes dropped idx = " << dropped_idx << std::endl;
      using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
      dt_type dt(std::move(new_ah), eigen_vector(x), succ,
                 ignore_begin, ignore_end);
      dt.add_ignore_idx(dropped_idx);
      dth(std::move(dt));
    }
  }
  // move-case for "first" (skipped) iteration
  // same loop body, except that we can move ah, instead of creating a copy of it
  size_type const dropped_idx = *(ah.begin() + *drop_pos_begin);
  swap_dropped(dropped_idx);
  if ((existing_cp = fc.find(cp_type(newidx.begin(), logical_end, 0)))) {
    Logger() << "CP ALREADY FOUND - NO DT SPAWNED,SUCCS UPDATED\n";
    existing_cp->add_successor(succ);
  } else {
    Logger() << "DT takes dropped idx = " << dropped_idx << std::endl;
    using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
    ah.drop_point(ah.begin() + *drop_pos_begin);
    dt_type dt(std::move(ah), eigen_vector(x), succ, ignore_begin, ignore_end);
    dt.add_ignore_idx(dropped_idx);
    dth(std::move(dt));
  }
}

}  // namespace FC

#endif  // COMMON_HPP_

