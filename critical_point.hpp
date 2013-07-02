#ifndef CRITICAL_POINT_HPP_
#define CRITICAL_POINT_HPP_

#include <algorithm>
#include <vector>

namespace FC {

class no_distance_policy {
};

template <typename number_type>
class distance_policy {
  public:
    number_type distance() const {
      return _distance;
    }

  private:
    number_type _distance;
};

// TODO(Lars): cp at inf
template <typename _idx_type, class _distance_policy>
class critical_point<_idx_type, _distance_policy> : public _distance_policy {
  public:
    typedef critical_point * succ_type;
  private:
    using idx_container  = std::vector<_idx_type>;
    using succ_container = std::vector<succ_type>;
  public:
    typedef _idx_type                     idx_type;
    typedef idx_container::iterator   idx_iterator;
    typedef succ_container::iterator succ_iterator;

    /** the indices are guaranteed to be sorted in ascending order */
    idx_iterator idx_begin() const nothrow;
    idx_iterator idx_end() const nothrow;

    succ_iterator succ_begin() const nothrow;
    succ_iterator succ_end() const nothrow;

    void add_successor(succ_type);

  private:
    idx_container _indices;
    succ_container _succsessors;
};

template <typename _idx_type, class _distance_policy>
bool operator==(critical_point<_idx_type,_distance_policy> const& lhs,
                critical_point<_idx_type,_distance_policy> const& rhs) {
  // Note: this requires the index ranges to be sorted
  return std::equal(lhs.idx_begin(), lhs.idx_end(),
                    rhs.idx_begin(), rhs.idx_end());
}

}  // namespace FC

namespace std {
  template <typename _idx_type, class _distance_policy>
  struct hash<critical_point<_idx_type,_distance_policy>> {
    using cp_type = critical_point<_idx_type,_distance_policy>;
    size_t operator()(cp_type const& cp) const {
      // the code below is boost's hash_range - thx boost
      size_t seed = 0;
      hash<typename cp_type::idx_type> hash_func;
      for (auto it = cp.idx_begin(); it != cp.idx_end(); ++it)
        seed ^= hash_func(*it) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
}  // namespace std

#endif  // CRITICAL_POINT_HPP_

