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
  @param it_begin OutputIterator to store at least ah.pc().dim() iterators
                  of type affine_hull<PointCloud>::iterator
  @param coeff_begin OutputIterator to store at least ah.pc().dim() coefficients
  @return true, if points were dropped, false otherwise
*/
template <typename PointCloud, typename Derived1, typename Derived2,
          typename Iterator>
Iterator drop_neg_coeffs(Eigen::MatrixBase<Derived1> const& x,
                         Eigen::MatrixBase<Derived2> const& lambda,
                         affine_hull<PointCloud> & ah, Iterator dropped_begin) {
  using size_type = typename affine_hull<PointCloud>::size_type;
  assert(static_cast<size_type>(lambda.size()) == ah.size());
  ah.project(x, const_cast<Eigen::MatrixBase<Derived2> &>(lambda));
  auto m_it = ah.end();
  auto dropped_end = dropped_begin;
  for (size_type i = static_cast<size_type>(lambda.size()); i-- > 0;) {
    --m_it;
    if (lambda[i] < 0) {  // TODO put stability threshold here
      std::cout << "NEG COEFF" << lambda[i] << " for index = " << *(ah.begin() + i) << std::endl;
      assert(ah.begin() <= m_it);
      assert(m_it < ah.end());
      *(dropped_end++) = *(ah.begin() + i);  // before drop, or iterator inval.
      ah.drop_point(m_it);
    } else {
      std::cout << "POS COEFF" << lambda[i] << " for index = " << *(ah.begin() + i) << std::endl;
    }
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
          typename size_type>
void spawn_sub_descends(DTHandler & dth,
                        size_type num_dropped,
                        Eigen::Matrix<number_type, Eigen::Dynamic, 1> && x,
                        affine_hull<PointCloud> ah,
                        critical_point<number_type, size_type> * succ) {
  size_type dropped_idx;
  using dt = descend_task<PointCloud>;
  for (size_type offset = 1; offset < num_dropped; ++offset) {
    auto new_ah(ah);
    dropped_idx = *(new_ah.begin() + offset);
    new_ah.drop_point(new_ah.begin() + offset);
    dth(dt(std::move(new_ah), x, succ, &dropped_idx, &dropped_idx + 1));
  }
  assert(num_dropped >= 1);
  dropped_idx = *ah.begin();
  ah.drop_point(ah.begin());
  dth(dt(std::move(ah), std::move(x), succ, &dropped_idx, &dropped_idx + 1));
}

}  // namespace FC

#endif  // UTILITY_HPP_

