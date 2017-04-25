#ifndef DYNAMIC_QR_HPP_
#define DYNAMIC_QR_HPP_

#include <cassert>
#include <cstddef>
#include <cmath>

#include <algorithm>
#include <ostream>
#include <vector>
#include <numeric>

#include <Eigen/Core>
#include <glog/logging.h>

namespace FC {

template <typename _number_type>
class dynamic_qr;

template <typename number_type>
std::ostream & operator<<(std::ostream &, dynamic_qr<number_type> const&);

template <typename _number_type>
class dynamic_qr {
  using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
  using eigen_map = Eigen::Map<eigen_vector>;
  using eigen_mat = Eigen::Matrix<_number_type, Eigen::Dynamic, Eigen::Dynamic>;
public:
  typedef _number_type number_type;
  typedef std::size_t  size_type;

  /**
    @brief will create a dynamic qr for at most num_rows columns
  */    
  dynamic_qr(size_type num_rows)
    : _q(eigen_mat::Identity(num_rows, num_rows)),
      _r_raw(new number_type[num_rows * num_rows]),
      _r_begin(new number_type * [num_rows]), _r_end(_r_begin) {
    DLOG(INFO) << "**QR-CTOR " << this << std::endl;
    assert(num_rows > 0);
    std::vector<size_type> offsets(num_rows);
    std::iota(offsets.begin(), offsets.end(), 0);
    init_row_ptr(offsets.begin());
  }
  
  // copy-constructor
  dynamic_qr(dynamic_qr const& orig)
    // note: orig.num_rows() == maximum number of columns possible
    : _q(orig._q), _r_raw(nullptr),
      _r_begin(new number_type * [orig.num_rows()]),
      _r_end(_r_begin + orig.num_cols()) {
    DLOG(INFO) << "**QR-COPY-CTOR " << this << std::endl;
    size_type const num_elements = orig.num_rows() * orig.num_rows();
    _r_raw = new number_type[num_elements];

    // since some column pointers might have been swapped,
    // we need to copy the exact same offsets
    std::vector<size_type> offsets(num_rows());
    for (size_type i = 0; i < offsets.size(); ++i)
      offsets[i] = (orig._r_begin[i] - orig._r_raw) / orig.num_rows();
    init_row_ptr(offsets.begin());
    
    std::copy(orig._r_raw, orig._r_raw + num_elements, _r_raw);
    
    CHECK(num_rows() == orig.num_rows());
    CHECK(num_cols() == orig.num_cols());
    CHECK(_q == orig._q);
    DCHECK(std::equal(_r_raw, _r_raw + num_elements, orig._r_raw));
  }

  // move-constructor
  dynamic_qr(dynamic_qr && tmp)
    : _q(), _r_raw(tmp._r_raw), _r_begin(tmp._r_begin), _r_end(tmp._r_end) {
    DLOG(INFO) << "**QR-MOVE-CTOR " << this << std::endl;
    _q.swap(tmp._q);
    // ownership has been transferred
    tmp._r_raw = nullptr;
    tmp._r_begin = nullptr;
  }
  
  // move-assignment
  dynamic_qr & operator=(dynamic_qr && tmp) {
    DLOG(INFO) << "**QR-MOVE-ASSIGN " << this << std::endl;
    if (this != &tmp) {
      _q.swap(tmp._q);
      _r_raw = tmp._r_raw;
      _r_begin = tmp._r_begin;
      _r_end = tmp._r_end;
      // ownership has been transferred
      tmp._r_raw = nullptr;
      tmp._r_begin = nullptr;
    }
    return *this;
  }

  dynamic_qr & operator=(dynamic_qr const&) = delete;
  
  size_type num_rows() const noexcept {
    assert(_q.rows() >= 0);
    return _q.rows();
  }
  
  size_type num_cols() const noexcept {
    assert(_r_end >= _r_begin);
    return _r_end - _r_begin;
  }
  
  template <typename DerivedVector>
  void append_column(Eigen::MatrixBase<DerivedVector> const& col) {
    assert(num_cols() < num_rows());
    eigen_map new_col(*_r_end, num_rows());
    new_col = _q.transpose() * col;
    number_type c, s;
    for (size_type i = num_rows() - 1; i-- > num_cols();) {
      // i is the index of 'a' in the givens routine
      givens(new_col[i], new_col[i + 1], &c, &s);  // c and s
      // Note: we don't use Eigen here for the rotation matrix multiplication
      //       because due to aliasing, we would have to evaluate the rotation
      //       result into a temporary and assign afterwards.
      update_q(c, s, i);
      apply_half_givens(c, s, *_r_end, i);  // update R
      // new_col[i + 1] is implicitely set to 0
    }
    ++_r_end;
  }

  void delete_column(size_type pos) {
    assert(pos >= 0);
    assert(pos < num_cols());
    if (pos != num_cols() - 1) {
      // move all columns right of pos one to the left
      number_type * tmp = _r_begin[pos];
      for (size_type i = pos, max = num_cols() - 1; i < max; ++i)
        _r_begin[i] = _r_begin[i + 1];
      --_r_end;  // num_cols() is now up-to-date again
      *_r_end = tmp;  // restored the address of the deleted column
      // R is now upper Hessenberg starting at column pos
      hessenberg_update(pos, num_cols());
    } else {
      --_r_end;
    }
  }
  
  /**
    @brief updates the QR decomposition for A + u * v^T
  */
  template <class DerivedVector1,
            class DerivedVector2>
  void rank_one_update (Eigen::MatrixBase<DerivedVector1> const& u,
                       Eigen::MatrixBase<DerivedVector2> const& v) {
    assert(size_type(u.rows()) == num_rows());
    assert(size_type(v.rows()) == num_cols());
    
    // we now need to actually zero the subdiagonal entries, since we
    // access those during the updates of R
    for (size_type i = 0; i < (num_cols() - 1); ++i)
      _r_begin[i][i + 1] = 0;
    bool const is_not_square = num_rows() > num_cols();
    if (is_not_square)  // only the the last column has a sub-diagonal element
      _r_begin[num_cols() - 1][num_cols()] = 0;
    
    eigen_vector w = _q.transpose() * u;
    number_type c ,s;
    for (size_type i = num_rows() - 1; i-- > 0;) {
      givens(w[i], w[i + 1], &c, &s);
      apply_half_givens(c, s, w.data(), i);  // clear element i + 1 of w
      if (i < num_cols())
        apply_givens(c, s, i, i);
      update_q(c, s, i);
    }
    // compute H = R + w * v^T
    // Since w = +-||w|| * e_1, this amounts to adding w[0] * v^t to R.row(0)
    for (size_type i = 0; i < num_cols(); ++i)
      _r_begin[i][0] += w[0] * v[i];
    // R is now upper Hessenberg
    hessenberg_update(0, (is_not_square ? num_cols() : num_cols() - 1));
  }
  
   /**
    @brief finds the least-squares solution to QR * x = b
  */
  template <typename Derived1, typename Derived2>
  void solve(Eigen::MatrixBase<Derived1> const& b,
             Eigen::MatrixBase<Derived2> const& x) const {
    assert(num_cols() > 0);
    auto & x_w = const_cast<Eigen::MatrixBase<Derived2> &>(x);
    assert(size_type(x_w.rows()) == num_cols());
    eigen_vector u = _q.transpose() * b;
    // now back substitution: R * x = u
    for (size_type i = num_cols(); i-- > 0;) {
      for (size_type j = i + 1; j < num_cols(); ++j)
        u[i] -= x_w[j] * _r_begin[j][i];
      x_w[i] = u[i] / _r_begin[i][i];
    }
  }
  
  ~dynamic_qr() {
    DLOG(INFO) << "**QR-DESTRUCT " << this << std::endl;
    if (_r_begin) {
      delete[] _r_begin;
      _r_begin = nullptr;
    }
    if (_r_raw) {      
      delete[] _r_raw;
      _r_raw = nullptr;
    }
  }

private:
  /**
    @brief helper function to avoid code duplication in constructors
    @param begin iterator to the offsets of the columns into the raw pointer
  */
  template <class Iterator>
  void init_row_ptr(Iterator begin) {
    size_type const max_cols = _q.cols();
    auto const& num_rows = max_cols;
    assert(nullptr != _r_raw);  // to make sure constructors initialized
                                // the storage
    for (size_type i = 0; i < max_cols; ++i)
      _r_begin[i] = _r_raw + *(begin++) * num_rows;
  }

  void givens(number_type const& a, number_type const& b,
              number_type * c, number_type * s) const {
    if (0.0 == b) {
      *c = 1.0;
      *s = 0.0;
    } else {
      if (abs(b) > abs(a)) {
        number_type const tau = -a / b;
        *s = 1 / sqrt(1 + (tau * tau));
        *c = *s * tau;
      } else {
        number_type const tau = -b / a;
        *c = 1 / sqrt(1 + (tau * tau));
        *s = *c * tau;
      }
    }
    DCHECK(abs(*c) <= 1);
    DCHECK(abs(*s) <= 1);
  }
  
  void update_q(number_type const& c, number_type const& s, size_type k) {
    assert(0 <= k);
    assert(typename eigen_mat::Index(k + 1) < _q.cols());

    number_type tmp;    
    for (typename eigen_mat::Index i = 0; i < _q.cols(); ++i) {
      tmp = c * _q(i, k) - s * _q(i, k + 1);
      _q(i, k + 1) = s * _q(i, k) + c * _q(i, k + 1);
      _q(i, k) = tmp;
    }
  }
  
  void apply_half_givens(number_type const& c, number_type const& s,
                         number_type * col, size_type a_idx) {
    assert(col);
    assert(0 <= a_idx);
    assert(a_idx < num_rows() - 1);

    col[a_idx] = c * col[a_idx] - s * col[a_idx + 1];
  }
  
  void apply_givens(number_type const& c, number_type const& s,
                    size_type r_col_begin, size_type a_idx) {
    assert(0 <= r_col_begin);
    assert(r_col_begin <= num_cols());
    assert(0 <= a_idx);
    assert(a_idx < num_rows() - 1);

    number_type tmp;
    for (size_type j = r_col_begin; j < num_cols(); ++j) {
      tmp = c * _r_begin[j][a_idx] - s * _r_begin[j][a_idx + 1];
      _r_begin[j][a_idx + 1] = s * _r_begin[j][a_idx] +
                               c * _r_begin[j][a_idx + 1];
      _r_begin[j][a_idx] = tmp;
    };
  }
  
  /**
    @param col_end the index of the column past the last column to update
  */
  void hessenberg_update(size_type k, size_type col_end) {
    number_type c, s;
    for (size_type i = k; i < col_end; ++i) {
      givens(_r_begin[i][i], _r_begin[i][i + 1], &c, &s);
      apply_half_givens(c, s, _r_begin[i], i);
      apply_givens(c, s, i + 1, i);
      update_q(c, s, i);
    }
  }

  template <typename _nr_type>
  friend std::ostream & operator<<(std::ostream &, dynamic_qr<_nr_type> const&);

  eigen_mat   _q;
  number_type *  _r_raw;
  number_type ** _r_begin;
  number_type ** _r_end;
};

template <typename number_type>
std::ostream & operator<<(std::ostream & ostr,
                          dynamic_qr<number_type> const& dqr) {
  using size_type = typename dynamic_qr<number_type>::size_type;
  ostr << "Q = \n" << dqr._q << '\n';
  ostr << "R = \n";
  for (size_type i = 0; i < dqr.num_rows(); ++i) {
    for (size_type j = 0; j < dqr.num_cols(); ++j) {
      if (i > j)
        ostr << 0 << '\t';
      else
        ostr << dqr._r_begin[j][i] << '\t';
    }
    ostr << '\n';
  }
  return ostr;
}

}  // namespace FC

#endif  // DYNAMIC_QR_HPP_

