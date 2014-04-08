#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <type_traits>

#include <Eigen/Core>

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

}  // namespace FC

#endif  // UTILITY_HPP_

