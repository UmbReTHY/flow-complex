#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <algorithm>
#include <type_traits>
#include <vector>

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

// TODO make succ for descend tasks part of the ident. This avoids duplicate
//      descend tasks from infinity, after "spawned more-than-one-neg" ascend tasks
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

