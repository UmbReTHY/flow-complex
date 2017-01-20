#include <cstdint>

#include <fstream>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>

#include "file_io.hpp"
#include "compute.hpp"
#include "clean.hpp"

int main(int argc, char ** argv) {
  try {
    using float_t = double;
    if (argc != 2)
      throw std::invalid_argument("usage: fc <point_cloud_file>");
    auto * filename = argv[1];
    std::ifstream in_file(filename);
    if (not in_file)
      throw std::runtime_error(std::string("could not open ") + filename);
    FC::point_store<float_t> ps;
    in_file >> ps;
    if (not ps.size())
      throw std::invalid_argument("empty data sets");

    using size_type = std::int64_t;
    size_type const dim = ps[0].size();
    auto fc = FC::compute_flow_complex<size_type>(ps.begin(), ps.end(), dim);
    if (not FC::validate(fc))
      std::cout << "warning: the computed flow complex is not valid. "
                   "This is probably the result of numerical inaccuracies or "
                   "degenerate input\n";
    // cleansing
    fc = clean_incidences(std::move(fc));
    // printing
    auto fc_filename = std::string(filename) + ".fc";
    std::ofstream f(fc_filename);
    if (not f)
      throw std::runtime_error("could not write to file " + fc_filename);
    f << fc;
  } catch (std::exception & e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
  return 0;
}

