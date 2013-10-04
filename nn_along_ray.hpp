#ifndef NN_ALONG_RAY_HPP_
#define NN_ALONG_RAY_HPP_

#include <cassert>

#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "point_cloud.hpp"

namespace fc {

template <typename number_type, typename size_type, typename dim_type,
          typename iterator>
std::pair<number_type, std::vector<size_type>>
nearest_neighbor_along_ray(number_type const* loc,
                           number_type const* ray,
                           point_cloud<number_type, size_type, dim_type> const& pc, // TODO templ templ
                           iterator members_begin, iterator members_end,
                           size_type const p_idx) {
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using eigen_cmap = Eigen::Map<const eigen_vector>;
  eigen_cmap v(ray, pc.dim());
  eigen_cmap p(pc[p_idx], pc.dim());
  eigen_cmap x(loc, pc.dim());
  number_type const x_p = x.dot(p);
  number_type const p_p = p.dot(p);
  number_type const v_p = v.dot(p);
  number_type min_t;
  std::vector<size_type> nn;
  // TODO there are possible improvements to reduce the number of iterations
  for (size_type i = 0; i < pc.size(); ++i) {
    if (members_end == std::find(members_begin, members_end, i)) {
      eigen_cmap q(pc[i], pc.dim());
      number_type const tmp = (v_p - q.dot(v));
      if (tmp.array().all() == 0)
        throw std::runtime_error("ill-conditioned data: div by 0");
      number_type const t = (2 * (x.dot(q) - x_p) + p_p - q.dot(q)) /
                            (2 * tmp);
      if (t > 0) {
        if (nn.empty() || t == min_t) {
          nn.push_back(i);
        } else if (t < min_t) {
          nn.clear();
          nn.push_back(i);
        }
      } 
    }
  }
  
  return std::make_pair(std::move(min_t), std::move(nn));
}

}  // namespace fc

#endif  // NN_ALONG_RAY_HPP_

