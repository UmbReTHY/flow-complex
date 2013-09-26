#ifndef CRITICAL_POINT_HPP_
#define CRITICAL_POINT_HPP_

#include <algorithm>
#include <vector>

#include "crit_point_base.hpp"
#include "flow_point.hpp"

namespace FC {

// TODO(Lars): cp at inf
//             remark: a critical point as singleton is potentially bad:
//             multiple flow complexes of different dimension would have the
//             same cp at inf. Needs to be lightweight (e.g. a single member
//             that depends on the dimension)

template <typename _idx_type>
class critical_point<_idx_type> : public crit_point_base<_idx_type> {
  public:
  private:
    using succ_container = std::vector<succ_type>;
  public:
    typedef _idx_type                           idx_type;
    typedef succ_container::const_iterator succ_iterator;

    critical_point(flow_point);  // TODO(Lars): template parameters still missing

    succ_iterator succ_begin() const nothrow;
    succ_iterator succ_end() const nothrow;

    void add_successor(succ_type);

  private:
    critical_point();

    succ_container _succsessors;
};

template <typename _idx_type>
bool operator==(critical_point<_idx_type> const& lhs,
                critical_point<_idx_type> const& rhs) {
  return std::is_permutation(lhs.idx_begin(), lhs.idx_end(), rhs.idx_begin());
}

}  // namespace FC



#endif  // CRITICAL_POINT_HPP_

