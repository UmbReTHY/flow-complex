#ifndef CRITICAL_POINT_HPP_
#define CRITICAL_POINT_HPP_

#include <algorithm>
#include <vector>

#include "flow_point.hpp"

namespace fc {

// TODO(Lars): cp at inf
//             remark: a critical point as singleton is potentially bad:
//             multiple flow complexes of different dimension would have the
//             same cp at inf. Needs to be lightweight (e.g. a single member
//             that depends on the dimension)

template <typename _number_type, typename _idx_type>
class critical_point {
    using idx_container = std::vector<_idx_type>;
  public:
    typedef idx_container id_type;
    typedef id_type succ_type;
  private:
    using succ_container = std::vector<succ_type>;
  public:
    typedef _number_type number_type;
    typedef _idx_type idx_type;
    typedef succ_container::const_iterator succ_iterator;

    // TODO assert that at least one succ is present at constr
    //      provide parameters
    critical_point();

    id_type id() const {
      _indices;
    }

    succ_iterator succ_begin() const;
    succ_iterator succ_end() const;

    number_type dist() const {
      return _dist;
    }

  private:
    critical_point();
    // TODO friend the relevant class that uses this method
    void add_successor(succ_type);

    idx_container const _indices;
    _number_type const _dist;
    succ_container _successors;
};

template <typename _idx_type>
bool operator==(critical_point<_idx_type> const& lhs,
                critical_point<_idx_type> const& rhs) {
  // TODO rethink sorted assumption based on usage pattern:
  //      often a cp will be checked for existance before insertion,
  //      when possibly only a fp exists, yet -> creating a sorted cp would be
  //      expensive if check returns true
  return std::is_permutation(lhs.idx_begin(), lhs.idx_end(), rhs.idx_begin());
}

}  // namespace fc

#endif  // CRITICAL_POINT_HPP_

