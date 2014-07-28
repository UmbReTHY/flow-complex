#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <ctime>

#include <iostream>
#include <vector>
#include <random>

#include <Eigen/Core>
#include <Eigen/QR>

#include "dynamic_qr.hpp"

using namespace FC;

using number_type = float;
using size_type = std::uint8_t;

size_type const NUM_ROWS = 4;
size_type const NUM_COLS = 3;

// TODO convert to CATCH

int main() {
  using dyn_qr_t = dynamic_qr<number_type>;
  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  
  std::srand(std::time(nullptr));
  
  eigen_matrix static_mat = eigen_matrix::Random(NUM_ROWS, NUM_COLS);
    
  dyn_qr_t dyn_qr(NUM_ROWS);
  eigen_vector x(NUM_COLS);
  eigen_vector b = eigen_vector::Random(NUM_ROWS);

  // append test
  for (size_type i = 0; i < static_mat.cols(); ++i) {
    dyn_qr.append_column(static_mat.col(i));
    assert(dyn_qr.num_cols() == (i + 1));
    dyn_qr.solve(b, x.head(i + 1));
    eigen_vector const static_sol = static_mat.leftCols(i + 1).householderQr().solve(b);
    assert(static_sol.isApprox(x.head(i + 1)));
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> unitRand;
  // delete test
  while (dyn_qr.num_cols() > 0) {
    // delete columns in random order
    auto const pos = std::floor(dyn_qr.num_cols() * unitRand(gen));
    {  // update the static matrix
      auto const num_cols_after = static_mat.cols() - pos - 1;
      static_mat.block(0, pos, NUM_ROWS, num_cols_after) = static_mat.rightCols(num_cols_after);
      static_mat.conservativeResize(Eigen::NoChange_t(), static_mat.cols() - 1);
    }
    // get the dynamic solution
    dyn_qr.delete_column(pos);
    assert(eigen_matrix::Index(dyn_qr.num_cols()) == static_mat.cols());
    if (dyn_qr.num_cols() > 0) {
      dyn_qr.solve(b, x.head(dyn_qr.num_cols()));
      eigen_vector const static_sol = static_mat.householderQr().solve(b);
      assert(static_sol.isApprox(x.head(dyn_qr.num_cols())));
    }
  }
  
  // generate a new random matrix
  static_mat = eigen_matrix::Random(NUM_ROWS, NUM_COLS);
  
  // second append test - to check whether this still works after deletion
  for (size_type i = 0; i < static_mat.cols(); ++i) {
    dyn_qr.append_column(static_mat.col(i));
    dyn_qr.solve(b, x.head(i + 1));
    eigen_vector const static_sol = static_mat.leftCols(i + 1).householderQr().solve(b);
    assert(static_sol.isApprox(x.head(i + 1)));
  }
  
  // size-getters test
  assert(dyn_qr.num_cols() == NUM_COLS);
  assert(dyn_qr.num_rows() == NUM_ROWS);
  
  // rank-one update
  eigen_vector const u = eigen_vector::Random(NUM_ROWS);
  eigen_vector const v = eigen_vector::Random(NUM_COLS);
  dyn_qr.rank_one_update(u, v);
  static_mat += u * v.transpose();
  dyn_qr.solve(b, x);
  eigen_vector const static_sol = static_mat.householderQr().solve(b);
  assert(static_sol.isApprox(x));
  
  std::exit(EXIT_SUCCESS);
}
