#ifndef CRITICAL_POINT_HPP_
#define CRITICAL_POINT_HPP_

#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>
#include <mutex>
#include <vector>

namespace fc {

// forward declarations
template <typename _number_type, typename _size_type>
class critical_point;
template <typename number_type, typename size_type>
size_type index(critical_point<number_type, size_type> const&);
template <typename number_type, typename size_type>
bool is_cp_at_inf(critical_point<number_type, size_type> const& cp);

template <typename _number_type, typename _size_type>
class critical_point {
    using self_t = critical_point<_number_type, _size_type>;
    using cp_cref_t = std::reference_wrapper<self_t const>;
    using succ_container = std::vector<cp_cref_t>;
    using idx_container = std::vector<_size_type>;
  public:
    typedef _number_type                           number_type;
    typedef _size_type                               size_type;
    typedef idx_container::const_iterator   idx_const_iterator;
    typedef succ_container::const_iterator succ_const_iterator;
 
    number_type dist() const noexcept {
      return _dist;
    }
    
    idx_const_iterator m_cbegin() const noexcept {
      return _indices.cbegin();
    }
    
    idx_const_iterator m_cend() const noexcept {
      return _indices.cend();
    }

    succ_const_iterator s_cbegin() const noexcept {
      return _succs.cbegin();
    }
    
    succ_const_iterator s_cend() const noexcept {
      return _succs.cend();
    }
  private:
    // TODO friend the factory
    template <typename member_iterator, typename succ_iterator>
    critical_point(member_iterator m_begin, member_iterator m_end,
                   succ_iterator s_begin, succ_iterator s_end,
                   number_type const distance)
    : _indices(m_begin, m_end), _succs(s_begin, s_end), _dist(distance) {
      if (_dist == std::numeric_limits<_number_type>::infinity()) {
        assert(_indices.empty() && "cps at inf shall not have members");
        assert(_succs.empty() && "cps at inf don't have successors");
      } else {
        assert(_dist >= 0);
        assert(!_indices.empty() && "only cps at inf can have no members");
        if (_dist == 0)
          assert(_indices.size() == 1u &&
                 "only index-0 cps have only 1 member");
        for (auto const& succ : _succs)
          assert((is_cp_at_inf(succ) || index(succ) > index(*this)) &&
                 _dist < succ._dist);
      }
    }

    // TODO friend the relevant class that uses this method
    void add_successor(self_t const& succ) {    
      _succs_mutex.lock();
      if (_succs.end() == std::find(_succs.begin(), _succs.end(), succ))
        _succs.push_back(succ);
      _succs_mutex.unlock();
    }

    idx_container const _indices;
    _number_type const _dist;
    succ_container _succs;
    std::mutex _succs_mutex;
};

template <typename number_type, typename size_type>
bool is_cp_at_inf(critical_point<number_type, size_type> const& cp) {
  return cp.dist() == std::numeric_limits<number_type>::infinity();
}

template <typename number_type, typename size_type>
size_type index(critical_point<number_type, size_type> const& cp) {
  // this implicates, that if the cp at inf is ecoded by an empty member list
  // this function cannot be called on a cp at inf
  assert(std::distance(cp.m_cbegin(), cp.m_cend()) > 0);
  return std::distance(cp.m_cbegin(), cp.m_cend()) - 1;
}

template <typename number_type, typename size_type>
bool operator==(critical_point<number_type, size_type> const& lhs,
                critical_point<number_type, size_type> const& rhs) {
  assert(std::is_sorted(lhs.m_cbegin(), lhs.m_cend()));
  assert(std::is_sorted(rhs.m_cbegin(), rhs.m_cend()));
  // we want to allow cps of different indices to be compared
  // (e.g. when stored in the same hash container, == is possibly called on
  // two instances with different index), however this should evaluate to
  // false then
  return (std::distance(lhs.m_cbegin(), lhs.m_cend()) ==
          std::distance(rhs.m_cbegin(), rhs.m_cend())) &&
         std::equal(lhs.m_cbegin(), lhs.m_cend(), rhs.m_cbegin());
}

template <typename number_type, typename size_type>
bool operator!=(critical_point<number_type, size_type> const& lhs,
                critical_point<number_type, size_type> const& rhs) {
  return !(lhs == rhs);
}

}  // namespace fc

#endif  // CRITICAL_POINT_HPP_
