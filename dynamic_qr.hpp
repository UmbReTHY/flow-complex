#ifndef DYNAMIC_QR_HPP_
#define DYNAMIC_QR_HPP_

#include <cstddef>

#include <Eigen/Dense>

namespace FC {

template <typename _number_type>
class dynamic_qr {
  public:
    typedef _number_type number_type;
    
    template <typename dim_type>
    dynamic_qr(dim_type dim)
      : _A(dim, 0) {
    }
    
    std::size_t num_rows() const noexcept {
      return _A.rows();
    }
    
    std::size_t num_cols() const noexcept {
      return _A.cols();
    }
    
    void append_column(number_type const* col) {
    }
    
    template <typename idx_type>
    void delete_column(idx_type pos) {
    }
    
    void rank_one_update() {
    }
    
  private:
    Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic> _A;
    
};

}  // namespace FC

#endif  // DYNAMIC_QR_HPP_

