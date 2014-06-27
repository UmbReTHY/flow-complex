#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <algorithm>

#include <type_traits>

#include <Eigen/Core>

#include "affine_hull.hpp"
#include "descend_task.hpp"

namespace FC {

template <typename T>
struct base_t {
  private:
    using base1 = typename std::remove_reference<T>::type;
    using base2 = typename std::remove_cv<base1>::type;
  public:
  typedef base2 type;
};

template <bool Aligned> struct eigen_align {
  using eigen_align_t = decltype(Eigen::Aligned);
  static eigen_align_t const value;
};
template <> typename eigen_align<true>::eigen_align_t
const eigen_align<true>::value = Eigen::Aligned;
template <> typename eigen_align<false>::eigen_align_t
const eigen_align<false>::value = Eigen::Unaligned;

/**
  @return the postitions returned are in descending order,
          to assist drop_neg_coeffs
*/
template <class Iterator, class Derived1, class Derived2, class PointCloud>
Iterator get_neg_pos(Eigen::MatrixBase<Derived1> const& x,
                     Eigen::MatrixBase<Derived2> const& lambda_const,
                     affine_hull<PointCloud> const& ah,
                     Iterator neg_begin) {
  using lambda_type = Eigen::MatrixBase<Derived2>;
  assert(lambda_const.size() == ah.size());
  auto & lambda = const_cast<lambda_type &>(lambda_const);
  Iterator neg_end = neg_begin;
  ah.project(x, lambda);
  using size_type = typename affine_hull<PointCloud>::size_type;
  for (size_type i = static_cast<size_type>(lambda_const.size()); i-- > 0;) {
    if (lambda[i] < 0) {  // TODO put stability threshold here
      std::cout << "NEG COEFF " << lambda[i] << " for index = " << *(ah.begin() + i) << std::endl;
      *(neg_end++) = i;
    } else {
      std::cout << "POS COEFF " << lambda[i] << " for index = " << *(ah.begin() + i) << std::endl;
    }
  }
  return neg_end;
}

/**
  @param it_begin OutputIterator to store at least ah.pc().dim() iterators
                  of type affine_hull<PointCloud>::iterator
  @param coeff_begin OutputIterator to store at least ah.pc().dim() coefficients
  @return true, if points were dropped, false otherwise
*/
template <typename PointCloud, typename Derived1, typename Derived2,
          typename Iterator>
Iterator drop_neg_indices(Eigen::MatrixBase<Derived1> const& x,
                          Eigen::MatrixBase<Derived2> const& lambda_const,
                          affine_hull<PointCloud> & ah,
                          Iterator dropped_begin) {
  using lambda_type = Eigen::MatrixBase<Derived2>;
  lambda_type & lambda = const_cast<lambda_type &>(lambda_const);
  // resuse the dropped range for the negative positions
  auto dropped_end = get_neg_pos(x, lambda, ah, dropped_begin);
  using size_type = typename affine_hull<PointCloud>::size_type;
  size_type tmp_idx;
  for (auto it = dropped_begin; it != dropped_end; ++it) {
    // assert descending order of positions
    assert((it + 1) == dropped_end or (*it > *(it + 1)));
    tmp_idx = *(ah.begin() + *it);
    ah.drop_point(ah.begin() + *it);
    *it = tmp_idx;  // now it turned into the index, not the position
  }
  return dropped_end;
}

/**
  @brief schedules descend tasks for every point but one that is dropped in
         the range. One point is dropped on ah directly and is intended to be
         reused in other tasks instead.
  @param num_dropped the number of points for which to schedule a descend task
                     after being dropped, starting at ah.begin()
*/
template <class DTHandler, class PointCloud, typename number_type,
          typename size_type, class Iterator>
void spawn_sub_descends(DTHandler & dth,
                        Iterator drop_pos_begin, Iterator drop_pos_end,
                        Eigen::Matrix<number_type, Eigen::Dynamic, 1> && x,
                        affine_hull<PointCloud> ah,
                        critical_point<number_type, size_type> * succ) {
  using dt = descend_task<PointCloud>;
  // we skip the first position to use it as the "move-case" after the loop
  assert(drop_pos_begin != drop_pos_end);
  for (Iterator it = std::next(drop_pos_begin); it != drop_pos_end; ++it) {
    auto new_ah(ah);
    new_ah.drop_point(new_ah.begin() + *it);
    auto const& dropped_idx = *(ah.begin() + *it);
    std::cout << "DT takes dropped idx = " << dropped_idx << std::endl;
    using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
    dth(dt(std::move(new_ah), eigen_vector(x), succ,
           &dropped_idx, std::next(&dropped_idx)));
  }
  assert(drop_pos_begin != drop_pos_end);
  size_type dropped_idx = *(ah.begin() + *drop_pos_begin);
  std::cout << "DT takes dropped idx = " << dropped_idx << std::endl;
  ah.drop_point(ah.begin() + *drop_pos_begin);
  dth(dt(std::move(ah), std::move(x), succ,
         &dropped_idx, std::next(&dropped_idx)));
}

template <typename _size_type>
class circumsphere_ident {
  public:
    typedef _size_type size_type;
    
    template <class Iterator>
    circumsphere_ident(Iterator begin, Iterator end) : _support(begin, end) {
      std::sort(_support.begin(), _support.end());
    }
    
    bool operator==(circumsphere_ident const& rhs) const {
      return (_support.size() == rhs._support.size()) and
             std::equal(_support.begin(), _support.end(), rhs._support.begin());
    }
    
    bool operator!=(circumsphere_ident const& rhs) const {
      return not operator==(rhs);
    }
    
  private:
    std::vector<_size_type> _support;
};

}  // namespace FC

#endif  // UTILITY_HPP_

