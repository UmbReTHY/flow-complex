#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <ctime>

#include <iostream>
#include <vector>

#include "dynamic_qr.hpp"

using namespace fc;

using number_type = float;
using size_type = std::uint8_t;

size_type const NUM_ROWS = 7;
size_type const NUM_COLS = 5;

// TODO convert to CATCH

int main() {
  using dyn_qr_t = dynamic_qr<number_type, size_type>;
  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  
  std::srand(std::time(nullptr));
  
  eigen_matrix static_mat = eigen_matrix::Random(NUM_ROWS, NUM_COLS);
    
  dyn_qr_t dyn_qr(NUM_ROWS);
  std::vector<number_type> x(NUM_COLS);
  Eigen::Map<eigen_vector> x_map(x.data(), x.size());
  std::vector<number_type> b(NUM_ROWS);
  Eigen::Map<eigen_vector> b_map(b.data(), b.size());
  b_map = eigen_vector::Random(b.size());

  // append test
  for (size_type i = 0; i < static_mat.cols(); ++i) {
    dyn_qr.append_column(static_mat.col(i));
    dyn_qr.solve(b.data(), x.data());
    eigen_vector const static_sol = static_mat.leftCols(i + 1).householderQr().solve(b_map);
    assert(static_sol == x_map.head(i + 1));
  }

  // delete test
  while (dyn_qr.num_cols() > 0) {
    // delete columns in random order
    auto const pos = std::floor(dyn_qr.num_cols() *
                                (rand() / static_cast<float>(RAND_MAX)));
    {  // update the static matrix
      eigen_matrix new_static(NUM_ROWS, static_mat.cols() - 1);
      for (size_type i = 0; i < static_mat.cols(); ++i)
        if (i < pos)
          new_static.col(i) = static_mat.col(i);
        else if (i > pos)
         new_static.col(i - 1) = static_mat.col(i);
      static_mat = new_static;
    }
    // get the dynamic solution
    dyn_qr.delete_column(pos);
    if (dyn_qr.num_cols() > 0) {
      dyn_qr.solve(b.data(), x.data());
      eigen_vector const static_sol = static_mat.householderQr().solve(b_map);
      assert(static_sol == x_map.head(dyn_qr.num_cols()));
    }
  }
  
  // generate a new random matrix
  static_mat = eigen_matrix::Random(NUM_ROWS, NUM_COLS);
  
  // second append test
  for (size_type i = 0; i < static_mat.cols(); ++i) {
    dyn_qr.append_column(static_mat.col(i));
    dyn_qr.solve(b.data(), x.data());
    eigen_vector const static_sol = static_mat.leftCols(i + 1).householderQr().solve(b_map);
    assert(static_sol == x_map.head(dyn_qr.num_cols()));
  }
  
  // size-getters test
  assert(dyn_qr.num_cols() == NUM_COLS);
  assert(dyn_qr.num_rows() == NUM_ROWS);
  
  // rank-one update
  eigen_vector const u = eigen_vector::Random(NUM_ROWS);
  dyn_qr.rank_one_update(u, eigen_vector::Ones(NUM_COLS));
  static_mat += u.replicate(1, NUM_COLS);
  dyn_qr.solve(b.data(), x.data());
  eigen_vector const static_sol = static_mat.householderQr().solve(b_map);
  assert(static_sol == x_map);
  
  std::exit(EXIT_SUCCESS);
}
