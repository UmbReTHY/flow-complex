#include <cstddef>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <stdexcept>
#include <exception>
#include <tuple>
#include <fstream>
#include <string>

#include <Eigen/Core>

#include "file_io.hpp"
#include "point_cloud.hpp"

int main(int argc, char ** argv) {
try {
  if (3 != argc)
    throw std::invalid_argument("usage: perturb <in_file> <out_file>");
  // read in_file
  using float_t = double;
  char const* in_filename = argv[1];
  std::ifstream in_file(in_filename);
  if (not in_file)
    throw std::runtime_error(std::string("could not open ") + in_filename);
  FC::point_store<float_t> ps;
  in_file >> ps;
  // perform pertubation
  if (ps.size()) {
    using size_type = std::size_t;
    auto pc = FC::point_cloud<float_t, size_type, false>(ps.cbegin(), ps.cend(),
                                                         ps[0].size());
    float_t const FACTOR = 0.1;
    std::vector<size_type> nn(2);
    std::vector<float_t> sq_dists(2);
    for (auto & pv : ps) {
      using EigenVector = Eigen::Matrix<float_t, Eigen::Dynamic, 1>;
      auto p = typename EigenVector::MapType(pv.data(), pc.dim());
      pc.k_nearest_neighors(p, 2, nn.data(), sq_dists.data());
      auto dist = sq_dists[1];
      p += (FACTOR * dist) * EigenVector::Random(pc.dim()).normalized();
    }
  }
  // write out_file
  char const* out_filename = argv[2];
  std::ofstream out_file(out_filename);
  if (not out_filename)
    throw std::runtime_error(std::string("could not open ") + out_filename);
  out_file << ps;
} catch (std::exception & e) {
  std::cout << e.what() << std::endl;
  std::exit(EXIT_FAILURE);
}
  std::exit(EXIT_SUCCESS);
}
