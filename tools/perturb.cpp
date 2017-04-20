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

int main(int argc, char ** argv) {
try {
  if (4 != argc)
    throw std::invalid_argument("usage: perturb <eps> <in_file> <out_file>");
  // read in_file
  using float_t = double;
  using size_type = int;
  char const* in_filename = argv[2];
  std::ifstream in_file(in_filename);
  if (!in_file)
    throw std::runtime_error(std::string("could not open ") + in_filename);
  FC::point_store<float_t, size_type> ps;
  in_file >> ps;
  // perform pertubation
  float_t const FACTOR = std::atof(argv[1]);
  if (ps.size()) {
    std::vector<size_type> nn(2);
    std::vector<float_t> sq_dists(2);
    for (auto & pv : ps) {
      using EigenVector = Eigen::Matrix<float_t, Eigen::Dynamic, 1>;
      auto p = typename EigenVector::MapType(pv.data(), pv.size());
      p += FACTOR * EigenVector::Random(pv.size()).normalized();
    }
  }
  // write out_file
  char const* out_filename = argv[3];
  std::ofstream out_file(out_filename);
  if (!out_filename)
    throw std::runtime_error(std::string("could not open ") + out_filename);
  out_file << ps;
} catch (std::exception & e) {
  std::cout << e.what() << std::endl;
  std::exit(EXIT_FAILURE);
}
  std::exit(EXIT_SUCCESS);
}
