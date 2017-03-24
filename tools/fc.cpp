#include <cstdint>

#include <fstream>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>

#include <gflags/gflags.h>

#include "clean.hpp"
#include "compute.hpp"
#include "file_io.hpp"
#include "flow_complex.hpp"

DEFINE_string(point_cloud, "", "path to file containing a point cloud");
DEFINE_bool(hist, false, "flag to toggle a print of the histogram of"
                         " computed critical points");
DEFINE_int32(num_threads, -1, "number of threads to use");

int main(int argc, char ** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_point_cloud.empty()) << "point cloud file missing";
  try {
    using float_t = long double;
    std::ifstream in_file(FLAGS_point_cloud);
    if (!in_file)
      throw std::runtime_error("could not open " + FLAGS_point_cloud);
    FC::point_store<float_t> ps;
    in_file >> ps;
    if (0 == ps.size()) throw std::invalid_argument("empty data sets");
    using size_type = std::int64_t;
    size_type const dim = ps[0].size();
    auto fc = FC::compute_flow_complex<size_type>(ps.begin(), ps.end(), dim,
                                                  FLAGS_num_threads);
    if (not FC::validate(fc))
      std::cout << "warning: the computed flow complex is not valid. "
                   "This is probably the result of numerical inaccuracies or "
                   "degenerate input\n";
    if (FLAGS_hist) {
      std::cout << "** printing histogram **" << std::endl;
      std::cout << "index\tcount" << std::endl;
      for (const auto& cp_pair : FC::compute_hist(fc))
        std::cout << cp_pair.first << "\t" << cp_pair.second << std::endl;
    }
    // cleansing
    fc = clean_incidences(std::move(fc));
    // printing
    auto fc_filename = FLAGS_point_cloud + ".fc";
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

