#ifndef NN_ALONG_RAY_HPP_
#define NN_ALONG_RAY_HPP_

#include <cassert>
#include <cstdint>

#include <limits>
#include <tuple>

#include <Eigen/Dense>

#include "point_cloud.hpp"

namespace FC {

enum class NN_ERROR : std::uint8_t {
  OK,
  DIV_BY_ZERO,  // solved by perturbing data
  MORE_THAN_ONE_NN  // solved by perturbing data
};

/**
  @tparam F callable taking a pointer to size_type and returning a bool
  @param loc point with the current location
  @param pc the point cloud
  @param p_idx index of a point of the point cloud, which is at starting
               distance from loc
  @param filter function that, given an index into pc, returns wheter or not
                this point does even need to be considered as a candidate
*/
template <typename number_type, typename size_type, typename F>
std::tuple<size_type, number_type, NN_ERROR>
nearest_neighbor_along_ray(number_type const* loc,
                           number_type const* ray,
                           number_type const* member,
                           size_type const dim,
                           F get_next) {
  static_assert(std::numeric_limits<number_type>::has_infinity,
                "number_type needs to have a representation for infinity");
  assert(m_end - m_begin > 0);
  assert(loc != nullptr);
  assert(ray != nullptr);
  assert(p_idx < pc.size());  // needs to be part of the point cloud
  assert(not filter(p_idx));  // should not allow it's own members
  
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using eigen_cmap = Eigen::Map<eigen_vector const>;
  // map the input to eigen vectors
  eigen_cmap x(loc, dim);
  eigen_cmap v(ray, dim);
  eigen_cmap p(pc[p_idx], dim);
  // compute some values that don't depend on the candidate points
  number_type const x_p = x.dot(p);
  number_type const p_p = p.dot(p);
  number_type const v_p = v.dot(p);
  // init return value
  auto r = std::make_tuple(size_type(0),
                           std::numeric_limits<number_type>::infinity(),
                           NN_ERROR::OK);
  size_type q_idx;
  while (get_next(&q_idx)) {
    eigen_cmap q(q_idx, dim);
    number_type const tmp = (v_p - q.dot(v));
    if (number_type(0) == 0) {
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
  }
  
  return r;
}

}  // namespace FC

#endif  // NN_ALONG_RAY_HPP_

