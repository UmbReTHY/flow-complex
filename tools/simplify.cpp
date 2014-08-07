#include <cstdlib>
#include <cstdint>
#include <cmath>

#include <fstream>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <utility>
#include <tuple>
#include <limits>
#include <vector>

#include "flow_complex.hpp"

template <typename nt, typename st>
FC::flow_complex<nt, st> simplify(FC::flow_complex<nt, st> fc, nt t);

int main(int argc, char ** argv) {
try {
  using float_t = double;
  using size_type = std::uint32_t;
  if (4 != argc)
    throw std::invalid_argument("usage: simplify <t> <in_file> <out_file>");
  // read the flow complex from file
  auto fc = FC::flow_complex<float_t, size_type>(42, 42);  // dummy
  {
    auto in_filename = argv[2];
    std::ifstream in_file(in_filename);
    in_file >> fc;
  }
  // simplify
  float_t const t = std::atof(argv[1]);
  if (t < 1)
    throw std::invalid_argument("prerequisite: t >= 1");
  fc = simplify(std::move(fc), t);
  // write result flow complex to file
  char const* out_filename = argv[3];
  std::ofstream out_file(out_filename);
  out_file << fc;
} catch (std::exception & e) {
  std::cerr << e.what() << std::endl;
  std::exit(EXIT_FAILURE);
}
  std::exit(EXIT_SUCCESS);
}

template <typename nt, typename st>
std::tuple<typename FC::flow_complex<nt, st>::cp_type *,
           typename FC::flow_complex<nt, st>::cp_type *,
           nt>
get_min_incidence(FC::flow_complex<nt, st> & fc) {
  auto r = std::make_tuple(&*fc.begin(), &*fc.begin(),
                           std::numeric_limits<nt>::infinity());
  auto & min_a = std::get<0>(r);
  auto & min_b = std::get<1>(r);
  auto & min_ratio = std::get<2>(r);
  for (auto & a : fc) {
    if (0 == a.index())
      continue;  // the ratio is not defined for a_dist == 0
    nt const a_dist = std::sqrt(a.sq_dist());
    for (auto b_it = a.succ_begin(); b_it != a.succ_end(); ++b_it) {
      auto & b = **b_it;
      if ((not b.is_max_at_inf()) and (1 == (b.index() - a.index()))) {
        nt const ratio_ab = std::sqrt(b.sq_dist()) / a_dist;
        if (ratio_ab < min_ratio) {
          min_a = &a;
          min_b = &b;
          min_ratio = ratio_ab;
        }
      }
    }
  }
  return r;
}

template <typename nt, typename st>
FC::flow_complex<nt, st> simplify(FC::flow_complex<nt, st> fc,
                                  nt t) {
  // used types
  using fc_type = FC::flow_complex<nt, st>;
  using cp_type = typename fc_type::cp_type;
  // tuple: predecessor a, successor b, ratio t_ab
  auto min_incidence = get_min_incidence(fc);
  auto & a_ptr = std::get<0>(min_incidence);
  auto & b_ptr = std::get<1>(min_incidence);
  auto & min_ratio = std::get<2>(min_incidence);
  std::vector<cp_type *> in_a;
  std::vector<cp_type *> out_a;
  std::vector<cp_type *> in_b;
  while (min_ratio <= t) {
    in_a.clear(); in_b.clear(); out_a.clear(); // init
    for (auto & cp : fc) {  // determine In(a) and In(b)
      for (auto it = cp.succ_begin(); it != cp.succ_end(); ++it) {
        if (*it == a_ptr) {
          in_a.push_back(&cp); break;
        } else if (*it == b_ptr) {
          in_b.push_back(&cp); break;
        }
      }
    }
    out_a.insert(out_a.end(), a_ptr->succ_begin(), a_ptr->succ_end());  // Out(a)
    // update the incidences
    for (auto ib : in_b)
      for (auto oa : out_a)
        if (ib->succ_end() == std::find(ib->succ_begin(), ib->succ_end(), oa))
          ib->add_successor(oa);
    // remove a and be and possible incidences to them
    for (auto ia : in_a)
      ia->erase(a_ptr);
    for (auto ib : in_b)
      ib->erase(b_ptr);
    fc.erase(*a_ptr);
    fc.erase(*b_ptr);
    // prepare next iteration
    min_incidence = get_min_incidence(fc);
  }
  return fc;
}
