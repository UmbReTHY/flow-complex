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
DEFINE_bool(bench, false, "only compute, but don't write to disk");

int main(int argc, char ** argv) {
  google::InitGoogleLogging(argv[0]);
  gflags::SetUsageMessage("call with --helpshort parameter for available flags");
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  CHECK(!FLAGS_point_cloud.empty()) << "point cloud file missing";
  try {
    using float_t = long double;
    using size_type = int;
    std::ifstream in_file(FLAGS_point_cloud);
    if (!in_file)
      throw std::runtime_error("could not open " + FLAGS_point_cloud);
    FC::point_store<float_t, size_type> ps;
    in_file >> ps;
    CHECK (ps.size() > 0) << "empty data sets";
    auto fc = FC::compute_flow_complex<size_type>(ps.begin(), ps.end(),
                                                  ps.dim(), FLAGS_num_threads);
    if (!FC::validate(fc))
      std::cout << "warning: the computed flow complex is not valid. "
                   "This is probably the result of numerical inaccuracies or "
                   "degenerate input\n";
    if (FLAGS_hist) {
      std::cout << "** printing histogram **" << std::endl;
      std::cout << "index\tcount" << std::endl;
      for (const auto& cp_pair : FC::compute_hist(fc))
        std::cout << cp_pair.first << "\t" << cp_pair.second << std::endl;
    }
    std::cout.flush();
    // cleansing
    fc = clean_incidences(std::move(fc));
    // printing
    if (!FLAGS_bench) {
      auto fc_filename = FLAGS_point_cloud + ".fc";
      std::ofstream f(fc_filename);
      if (!f)
        throw std::runtime_error("could not write to file " + fc_filename);
      f << fc;
    }
  } catch (std::exception & e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
  return 0;
}

