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
template <typename PointCloud, typename Derived1, typename Derived2>
bool drop_neg_coeffs(Eigen::MatrixBase<Derived1> const& x,
                     Eigen::MatrixBase<Derived2> const& lambda,
                     affine_hull<PointCloud> & ah) {
  using size_type = typename affine_hull<PointCloud>::size_type;
  assert(static_cast<size_type>(lambda.size()) == ah.size());
  ah.project(x, const_cast<Eigen::MatrixBase<Derived2> &>(lambda));
  auto m_it = ah.end();
  for (size_type i = static_cast<size_type>(lambda.size()); i-- > 0;) {
    --m_it;
    if (lambda[i] < 0) {
      std::cout << "NEG COEFF" << lambda[i] << std::endl;
      assert(ah.begin() <= m_it);
      assert(m_it < ah.end());
      ah.drop_point(m_it);
    }
  }
  return static_cast<size_type>(lambda.size()) != ah.size();
}

/**
  TODO this is not useful for descend tasks: one id should not be dropped
  @brief schedules descend tasks for every point but one that is dropped in
         the range. One point is dropped on ah directly and is intended to be
         reused in other tasks instead.
  @param num_dropped the number of points for which to schedule a descend task
                     after being dropped, starting at ah.begin()
*/
template <class DTHandler, class PointCloud, class Derived,
          typename number_type, typename size_type>
void spawn_sub_descends(DTHandler const& dth,
                        Eigen::MatrixBase<Derived> const& x,
                        size_type num_dropped,
                        critical_point<number_type, size_type> * succ,
                        affine_hull<PointCloud> & ah) {
  // the dt corresponding to offset = 0 is continued by *this
  for (size_type offset = 1; offset < num_dropped; ++offset) {
    auto new_ah(ah);
    new_ah.drop_point(new_ah.begin() + offset);
    using dt = descend_task<PointCloud>;
    using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
    dth(dt(std::move(new_ah), eigen_vector(x), succ));
  }
  assert(num_dropped >= 1);
  ah.drop_point(ah.begin());
}

}  // namespace FC

#endif  // UTILITY_HPP_

