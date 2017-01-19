#ifndef NN_ALONG_RAY_HPP_
#define NN_ALONG_RAY_HPP_

#include <cassert>

#include <stdexcept>
#include <utility>

#include <Eigen/Core>
#include <glog/logging.h>

#include "vertex_filter.hpp"

namespace FC {

/**
  @tparam F callable taking a pointer to std::size_t and returning a pointer
            to an eigen_map
  @param begin, end iterators iterators for nearest neighbor storage.
                    If the number of nearest neighbors found exceeds
                    the given range, an exception is thrown.
*/
template <class Derived1, class Derived2, class Derived3, class NNIterator,
          class... Params>
std::pair<NNIterator, typename Derived1::Scalar>
nearest_neighbor_along_ray(Eigen::MatrixBase<Derived1> const& x,
                           Eigen::MatrixBase<Derived2> const& v,
                           Eigen::MatrixBase<Derived3> const& p,
                           vertex_filter<Params...> & get_next,
                           NNIterator begin, NNIterator end) {
  LOG(INFO) << "NORM OF RAY " << v.squaredNorm() << std::endl;
  using number_type = typename Derived1::Scalar;
  // compute some values that don't depend on the candidate points
  number_type const x_p = x.dot(p);
  number_type const p_p = p.dot(p);
  number_type const v_p = v.dot(p);
  // init return value
  auto r = std::make_pair(begin, number_type(0));
  typename vertex_filter<Params...>::size_type q_idx;
  auto * q_ptr = get_next(&q_idx);
  while (q_ptr) {
    auto & q = *q_ptr;
    number_type const tmp = (v_p - q.dot(v));
    if (0 == tmp)
      throw std::logic_error("division by 0");
    number_type const t = (x.dot(q) - x_p + 0.5 * (p_p - q.dot(q))) / tmp;
    LOG(INFO) << "t = " << t << " for id = " << q_idx << std::endl;
    if (t > 0 and (r.first == begin or t <= r.second)) {
      LOG(INFO) << "t-DIFF = " << (t - r.second) << std::endl;
      if (r.first != begin and Eigen::internal::isApprox(t, r.second)) {
        if (r.first == end)
          throw std::logic_error("too many nearest neighbors");
      } else {
        r.first = begin;
        r.second = t;
      }
      *r.first++ = q_idx;
    }
    q_ptr = get_next(&q_idx);
  }
  
  return r;
}

}  // namespace FC

#endif  // NN_ALONG_RAY_HPP_

