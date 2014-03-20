#ifndef DYNAMIC_QR_HPP_
#define DYNAMIC_QR_HPP_

#include <cassert>

#include <Eigen/Dense>

namespace FC {

template <typename _number_type, typename _size_type>
class dynamic_qr {
  using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
  using eigen_matrix = Eigen::Matrix<_number_type, Eigen::Dynamic, Eigen::Dynamic>;

  public:
    typedef _number_type number_type;
    typedef _size_type size_type;
    
    dynamic_qr(size_type num_rows)
      : _A(num_rows, 0) {
      assert(num_rows > 0);
    }
    
    size_type num_rows() const noexcept {
      return _A.rows();
    }
    
    size_type num_cols() const noexcept {
      return _A.cols();
    }
    
    template <typename DerivedVector>
    void append_column(Eigen::MatrixBase<DerivedVector> const& col) {
      assert(num_cols() < num_rows());
      // TODO seems buggy, check that
//      _A.conservativeResize(num_rows(), num_cols() + 1);
//      eigen_matrix const tmp = _A;
//      _A.resize(num_rows(), num_cols() + 1);
//      _A.leftCols(num_cols() - 1) = tmp;
//      _A.rightCols(1) = col;
      eigen_matrix tmp(num_rows(), num_cols() + 1);
      tmp.leftCols(num_cols()) = _A;
      tmp.rightCols(1) = col;
      _A.swap(tmp);
    }
    
    void delete_column(size_type const pos) {
      assert(pos >= 0);
      assert(pos < num_cols());
      if (pos != num_cols() - 1) {
        auto const num_cols_to_move = num_cols() - pos - 1;
        _A.block(0, pos, num_rows(), num_cols_to_move) =
        _A.block(0, pos + 1, num_rows(), num_cols_to_move);
      }
      _A.conservativeResize(num_rows(), num_cols() - 1);
    }
    
    /**
      @brief updates the QR decomposition for A + u * v^T
    */
    template <class DerivedVector1,
              class DerivedVector2>
    void rank_one_update (Eigen::MatrixBase<DerivedVector1> const& u,
                         Eigen::MatrixBase<DerivedVector2> const& v) {
      static_assert(DerivedVector1::ColsAtCompileTime == 1,
                    "u needs to be a column vector");
      static_assert(DerivedVector2::ColsAtCompileTime == 1,
                    "v needs to be a column vector");
      assert(u.rows() == num_rows());
      assert(v.rows() == num_cols());
      _A += u * v.transpose();
    }
    
     /**
      @brief finds the least-sqaures solution to QR * x = b
    */
    void solve(number_type const* b, number_type * x) const {
      assert(num_cols() > 0);
      Eigen::Map<eigen_vector> x_map(x, num_cols());
      Eigen::Map<eigen_vector const> b_map(b, num_rows());
      x_map = _A.householderQr().solve(b_map);
    }
    
  private:
    // TODO not efficient yet
    Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic> _A;
};

}  // namespace FC

#endif  // DYNAMIC_QR_HPP_

