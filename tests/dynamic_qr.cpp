#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstdint>

#include "dynamic_qr.hpp"

using namespace fc;

using number_type = float;
using size_type = std::uint8_t;

size_type const NUM_ROWS = 7;
size_type const NUM_COLS = 5;

int main() {
  using dyn_qr_t = dynamic_qr<number_type, size_type>;
  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  
  eigen_matrix const static_mat = eigen_matrix::Random(NUM_ROWS, NUM_COLS);
    
  dyn_qr_t dyn_qr(NUM_ROWS);
  std::vector<number_type> x(num_rows);
  Eigen::
  std::vector<number_type> b(num_rows);
  Eigen::Map<eigen_vector> b_map(b.data(), b.size());
  b_map = eigen_vector::Random(b.size());
  // append test
  for (size_type i = 0; i < static_mat.cols(); ++i) {
    dyn_qr.append_column(static_mat.col(i).data());
    dyn_qr.solve(b.data(), x.data());
    auto const static_sol = static_mat._A.householderQr().solve(b_map);
    assert( );
  }
  
  std::exit(EXIT_SUCCESS);
}
