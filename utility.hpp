#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <algorithm>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "makros.h"

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

template <typename _size_type>
class circumsphere_ident {
using container_type = std::vector<_size_type>;
public:
  typedef _size_type                          size_type;
  typedef typename container_type::const_iterator const_iterator;
  
  template <class Iterator>
  circumsphere_ident(Iterator begin, Iterator end) : _support(begin, end) {
    std::sort(_support.begin(), _support.end());
  }
  
  bool operator==(circumsphere_ident const& rhs) const {
    assert(_support.size() == rhs._support.size());
    return (_support.size() == rhs._support.size()) and
           std::equal(_support.begin(), _support.end(), rhs._support.begin());
  }
  
  bool operator!=(circumsphere_ident const& rhs) const {
    return not operator==(rhs);
  }
  
  const_iterator cbegin() const {
    return _support.cbegin();
  }
  
  const_iterator cend() const {
    return _support.cend();
  }
  
private:
  std::vector<_size_type> _support;
};

struct RangeHash {
  template <class Iterator>
  __pure std::size_t operator()(Iterator begin, Iterator end) const {
    assert(std::is_sorted(begin, end));
    // the code below is boost's hash_range - thx boost
    std::size_t seed = 0;
    using size_type = typename base_t<decltype(*Iterator())>::type;
    std::hash<size_type> _idx_hash;
    while (begin != end)
      seed ^= _idx_hash(*begin++) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct CIHash {
  template <typename size_type>
  std::size_t operator()(circumsphere_ident<size_type> const& ci) const {
    RangeHash range_hash;
    return range_hash(ci.cbegin(), ci.cend());
  }
};

}  // namespace FC

#endif  // UTILITY_HPP_

