#include <cassert>
#include <cstdlib>
#include <cstdint>

#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <limits>
#include <iterator>
#include <utility>

#include <Eigen/Core>

#include "file_io.hpp"
#include "flow_complex.hpp"
#include "point_cloud.hpp"
#include "affine_hull.hpp"
#include "update_ray.hpp"

template <class Iterator, typename nt, typename st, bool aligned>
std::vector<nt> c_of_seb(Iterator begin, Iterator end,
                         FC::point_cloud<nt, st, aligned> const&);

int main(int argc, char ** argv) {
  try {
    using size_type = int;
    using float_t = double;
    using ps_type = FC::point_store<float_t, size_type>;
    using pc_type = FC::point_cloud<float_t, size_type, false>;
    using fc_type = FC::flow_complex<float_t, size_type>;
    // deal with program arguments
    if (6 != argc)
      throw std::invalid_argument("usage: fc2pc <from_idx> <to_idx> <in_pc> <in_fc> <out_pc>");
    int from_idx = std::atoi(argv[1]);
    int to_idx = std::atoi(argv[2]);
    if (to_idx < from_idx || from_idx < 0)
      throw std::invalid_argument("invalid index range");
    // read point cloud
    ps_type in_ps;
    {
      auto * in_pc_filename = argv[3];
      std::ifstream is(in_pc_filename);
      is >> in_ps;
    }
    // read flow complex
    fc_type fc(42, 42);  // dummy arguments
    {
      auto * in_fc_filename = argv[4];
      std::ifstream is(in_fc_filename);
      is >> fc;
    }
    // compute requested circum-centers
    ps_type out_ps;
    if (in_ps.size()) {
      const size_type dim = in_ps.dim();
      pc_type pc(in_ps.begin(), in_ps.end(), dim);
      for (auto const& cp : fc) {
        if ((!cp.is_max_at_inf()) &&
            from_idx <= cp.index() && cp.index() <= to_idx) {
          out_ps.add_point(dim) = c_of_seb(cp.idx_begin(), cp.idx_end(), pc);
          // correctness check
          using ev = Eigen::Matrix<float_t, Eigen::Dynamic, 1>;
          using em = ev::MapType;
          em new_cc(out_ps[out_ps.size() - 1].data(), pc.dim());
        }
      }
    }
    // write result pc
    auto * out_pc_filename = argv[5];
    std::ofstream os(out_pc_filename);
    os << out_ps;
  } catch (std::exception & e) {
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::exit(EXIT_SUCCESS);
}

template <class Iterator, typename nt, typename st, bool aligned>
std::vector<nt> c_of_seb(Iterator begin, Iterator end,
                         FC::point_cloud<nt, st, aligned> const& pc) {
  using ev = Eigen::Matrix<nt, Eigen::Dynamic, 1>;
  using em = typename ev::MapType;
  using pc_type = FC::point_cloud<nt, st, aligned>;
  using ah_type = FC::affine_hull<pc_type>;
  // init
  ev c = ev::Random(pc.dim());
  ah_type ah(pc);
  {
    st farthest_idx = 0; nt farthest_dist = 0;
    for (auto it = begin; it != end; ++it) {
      nt dist = (pc[*it] - c).squaredNorm();
      if (dist > farthest_dist) {
        farthest_dist = dist;
        farthest_idx = *it;
      }
    }
    ah.append_point(farthest_idx);
  }
  // helper lambda
  auto find_stopper = [begin, end] (ah_type const& ah, ev const& ray,
                                    ev const& c) {
    nt min_t = std::numeric_limits<nt>::infinity();
    st idx = 77;
    auto const& pc = ah.pc();
    auto const& p = pc[*ah.begin()];
    nt v_p = ray.dot(p);
    nt c_p = c.dot(p);
    nt p_p = p.dot(p);
    for (auto it = begin; it != end; ++it) {
      if (ah.end() == std::find(ah.begin(), ah.end(), *it)) {
        auto const& q = pc[*it];
        nt const tmp = (v_p - q.dot(ray));
        nt const t = (c.dot(q) - c_p + 0.5 * (p_p - q.dot(q))) / tmp;
        if (t < min_t) {
          min_t = t;
          idx = *it;
        }
      }
    }
    assert(end != std::find(begin, end, idx));
    return std::make_pair(idx, min_t);
  };
  // move to dual voronoi
  ev cc(pc.dim());
  ev lambda(pc.dim() + 1);
  ev ray(pc.dim());
  while (ah.size() < std::distance(begin, end)) {
    assert(ah.size() >= 1);
    FC::update_ray<FC::RAY_DIR::TO_DRIVER>(ah, c, lambda, cc, ray);
    auto rpair = find_stopper(ah, ray, c);
    c += rpair.second * ray;
    ah.append_point(rpair.first);
  }
  // return cc
  FC::update_ray<FC::RAY_DIR::TO_DRIVER>(ah, c, lambda, cc, ray);
  std::vector<nt> r(pc.dim());
  em(r.data(), r.size()) = cc;
  return r;
}
