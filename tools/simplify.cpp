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
#include <string>
#include <ostream>

#include "flow_complex.hpp"

template <typename nt, typename st>
FC::flow_complex<nt, st> simplify(FC::flow_complex<nt, st>, nt t);
template <typename nt, typename st>
void print_histogram(std::ostream &, FC::flow_complex<nt, st>);

int main(int argc, char ** argv) {
try {
  using float_t = double;
  using size_type = std::uint32_t;
  if (4 != argc)
    throw std::invalid_argument("usage: simplify <<t> | -p> <in_file> <out_file>");
  // read the flow complex from file
  auto fc = FC::flow_complex<float_t, size_type>(42, 42);  // dummy
  {
    auto in_filename = argv[2];
    std::ifstream in_file(in_filename);
    in_file >> fc;
  }
  char const* action = argv[1];
  if (std::string(action) == "-p") {
    char const* out_filename = argv[3];
    std::ofstream out_file(out_filename);
    print_histogram(out_file, std::move(fc));
  } else {
    // simplify
    float_t const t = std::atof(argv[1]);
    if (t < 1)
      throw std::invalid_argument("prerequisite: t >= 1");
    fc = simplify(std::move(fc), t);
    // write result flow complex to file
    char const* out_filename = argv[3];
    std::ofstream out_file(out_filename);
    out_file << fc;
  }
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

/**
  @brief Performs a single reduction step.
*/
template <typename nt, typename st>
FC::flow_complex<nt, st> reduce(FC::flow_complex<nt, st> fc,
                                FC::critical_point<nt, st> * a_ptr,
                                FC::critical_point<nt, st> * b_ptr) {
  // used types
  using fc_type = FC::flow_complex<nt, st>;
  using cp_type = typename fc_type::cp_type;
  static std::vector<cp_type *> in_a;
  static std::vector<cp_type *> out_a;
  static std::vector<cp_type *> in_b;
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
  return fc;
}

/**
  @brief Performs a series of reductions until the threshold is exceeded.
*/
template <typename nt, typename st>
FC::flow_complex<nt, st> simplify(FC::flow_complex<nt, st> fc,
                                  nt t) {
  // tuple: predecessor a, successor b, ratio t_ab
  auto min_incidence = get_min_incidence(fc);
  auto & a_ptr = std::get<0>(min_incidence);
  auto & b_ptr = std::get<1>(min_incidence);
  auto & min_ratio = std::get<2>(min_incidence);
  while (min_ratio <= t) {
    fc = reduce(std::move(fc), a_ptr, b_ptr);
    // prepare next iteration
    min_incidence = get_min_incidence(fc);
  }
  return fc;
}

template <typename nt, typename st>
void print_histogram(std::ostream & os, FC::flow_complex<nt, st> fc) {
  auto const INF = std::numeric_limits<nt>::infinity();
  auto min_incidence = get_min_incidence(fc);
  auto & a_ptr = std::get<0>(min_incidence);
  auto & b_ptr = std::get<1>(min_incidence);
  auto & min_ratio = std::get<2>(min_incidence);
  using fc_type = FC::flow_complex<nt, st>;
  auto print_hist = [&os] (nt ratio, fc_type const& fc) {
    auto hist = compute_hist(fc);
    os << ratio;
    for (auto & p : hist)
      os << " " << p.second;
    os << std::endl;
  };
  print_hist(0, fc);  // print initial fc
  while (INF != min_ratio) {
    fc = reduce(std::move(fc), a_ptr, b_ptr);
    print_hist(min_ratio, fc);
    // prepare next iteration
    min_incidence = get_min_incidence(fc);
  }
}

