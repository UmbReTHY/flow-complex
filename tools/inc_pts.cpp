#include <cstdlib>

#include <exception>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <ostream>
#include <fstream>
#include <utility>

#include "critical_point.hpp"
#include "flow_complex.hpp"

template <typename nt, typename st>
using inc_map = std::unordered_map<FC::critical_point<nt, st> const*,
                                   std::unordered_set<st>>;
template <typename nt, typename st>
inc_map<nt, st> incident_points(FC::flow_complex<nt, st> const&);
template <typename nt, typename st>
std::ostream & operator<<(std::ostream &, inc_map<nt, st> const&);

int main(int argc, char ** argv) {
try {
  using float_t = double;
  using size_type = std::uint32_t;
  if (3 != argc)
    throw std::invalid_argument("usage: inc-pts <in_file> <out_file>");
  // read the flow complex from file
  auto fc = FC::flow_complex<float_t, size_type>(42, 42);  // dummy
  {
    auto in_filename = argv[1];
    std::ifstream in_file(in_filename);
    in_file >> fc;
  }
  // compute incident points (index-0 critical points)
  auto incidences = incident_points(fc);
  // write incidence list to file
  char const* out_filename = argv[2];
  std::ofstream out_file(out_filename);
  out_file << incidences;
} catch (std::exception & e) {
  std::cerr << e.what() << std::endl;
  std::exit(EXIT_FAILURE);
}
  std::exit(EXIT_SUCCESS);
}

template <typename nt, typename st>
inc_map<nt, st> augment_inc_map(FC::critical_point<nt, st> const* cp,
                                st idx,
                                inc_map<nt, st> im) {
  using set_t = std::unordered_set<st>;
  if (cp->index() > 0) {
    auto rpair = im.emplace(cp, set_t());
    rpair.first->second.insert(idx);
  }   
  for (auto it = cp->succ_begin(); it != cp->succ_end(); ++it)
    if (!(*it)->is_max_at_inf())
      im = augment_inc_map(*it, idx, std::move(im));
  return im;
}

template <typename nt, typename st>
inc_map<nt, st> incident_points(FC::flow_complex<nt, st> const& fc) {
  inc_map<nt, st> im;
  for (st i = 0; i < fc.num_minima(); ++i)
    im = augment_inc_map(fc.minimum(i), i, std::move(im));
  return im;
}

template <typename nt, typename st>
std::ostream & operator<<(std::ostream & os, inc_map<nt, st> const& im) {
  for (auto const& p : im) {
    os << *p.first << "|";
    for (auto idx : p.second)
      os << " " << idx;
    os << std::endl;
  }
  return os;
}

