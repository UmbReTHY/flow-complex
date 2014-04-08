#ifndef NN_ALONG_RAY_HPP_
#define NN_ALONG_RAY_HPP_

#include <cassert>
#include <cstdint>
#include <cstddef>

#include <limits>
#include <tuple>

#include <Eigen/Core>

#include "point_cloud.hpp"

namespace FC {

enum class NN_ERROR : std::uint8_t {
  OK,
  DIV_BY_ZERO,  // solved by perturbing data
  MORE_THAN_ONE_NN  // solved by perturbing data
};

/**
  @tparam F callable taking a pointer to std::size_t and returning a pointer
            to an eigen_map
  @param loc point with the current location
  @param pc the point cloud
  @param p_idx index of a point of the point cloud, which is at starting
               distance from loc
  @param filter function that, given an index into pc, returns wheter or not
                this point does even need to be considered as a candidate
*/
template <typename Derived1, typename Derived2, typename Derived3,
          typename F>
std::tuple<std::size_t, typename Derived1::Scalar, NN_ERROR>
nearest_neighbor_along_ray(Eigen::MatrixBase<Derived1> const& x,
                           Eigen::MatrixBase<Derived2> const& v,
                           Eigen::MatrixBase<Derived3> const& p,
                           F & get_next) {
  using number_type = typename Derived1::Scalar;
  static_assert(std::numeric_limits<number_type>::has_infinity,
                "number_type needs to have a representation for infinity");
  // compute some values that don't depend on the candidate points
  number_type const x_p = x.dot(p);
  number_type const p_p = p.dot(p);
  number_type const v_p = v.dot(p);
  // init return value
  auto r = std::make_tuple(std::size_t(0),
                           std::numeric_limits<number_type>::infinity(),
                           NN_ERROR::OK);
  std::size_t q_idx;
  auto * q_ptr = get_next(&q_idx);
  while (q_ptr) {
    auto & q = *q_ptr;
    number_type const tmp = (v_p - q.dot(v));
    if (0 == tmp) {
      std::get<2>(r) = NN_ERROR::DIV_BY_ZERO;
      break;
    }
    number_type const t = (x.dot(q) - x_p + 0.5 * (p_p - q.dot(q))) / tmp;
    if (t > 0) {
      if (t < std::get<1>(r)) {
        std::get<0>(r) = q_idx;
        std::get<1>(r) = t;
      } else if (t == std::get<1>(r)) {
        std::get<2>(r) = NN_ERROR::MORE_THAN_ONE_NN;
        break;
      }
    }
    q_ptr = get_next(&q_idx);
  }
  
  return r;
}

}  // namespace FC

#endif  // NN_ALONG_RAY_HPP_

