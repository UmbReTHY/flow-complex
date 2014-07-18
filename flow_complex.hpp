#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <algorithm>
#include <utility>
#include <map>

#include <tbb/concurrent_unordered_set.h>

#include "critical_point.hpp"
#include "logger.hpp"

namespace FC {

template <typename _number_type, typename _size_type>
class flow_complex {
  public:
    typedef _number_type                             number_type;
    typedef _size_type                                 size_type;
    typedef critical_point<_number_type, _size_type>     cp_type;

  private:
    using cp_container = tbb::concurrent_unordered_set<cp_type, CPHash>;
    
  public:
    typedef typename cp_container::iterator             iterator;
    typedef typename cp_container::const_iterator const_iterator;
  
    /**
      @brief initializes the flow complex with a maximum at infinity
    */
    flow_complex(size_type dim, size_type num_pts)
    : _minima(num_pts, nullptr) {
      _max_at_inf = &*_cps.insert(cp_type(dim)).first;
      // init fc with id-0 critical points
      for (size_type i = 0; i < num_pts; ++i)
        _minima[i] = insert(cp_type(&i, &i + 1, 0)).second;
    }
    // copy- and move-constructor
    flow_complex(flow_complex const&) = default;
    flow_complex(flow_complex &&) = default;
    // copy- and move-constructor
    flow_complex & operator=(flow_complex const&) = default;
    flow_complex & operator=(flow_complex &&) = default;
  
    // iterators
    iterator begin() {
      return _cps.begin();
    }
    
    iterator end() {
      return _cps.end();
    }

    // const-iterators
    const_iterator begin() const {
      return _cps.begin();
    }

    const_iterator end() const {
      return _cps.end();
    }
    
    // modifiers
    void erase(cp_type const&);
    void erase(iterator);
    void erase(const_iterator);
    
    /**
      @return the first member of the pair is true, if the critical point was
              inserted by this call, and false if it already existed
    */
    std::pair<bool, cp_type *> insert(cp_type && cp) {
      std::pair<bool, cp_type *> r;
      auto ret_pair = _cps.insert(cp);
      r.first =   ret_pair.second;
      r.second = &*ret_pair.first;
      if (r.first)
        Logger() << "CRITICAL POINT IS NEW!\n";
      else
        Logger() << "CRITICAL POINT WAS FOUND ALREADY\n";
      return r;
    }
    
    cp_type const* max_at_inf() const {
      return _max_at_inf;
    }
    
    cp_type * minimum(size_type pos) const {
      return _minima[pos];
    }
    
    cp_type * find(cp_type const& cp) const {
      auto it = _cps.find(cp);
      return (it == _cps.end() ? nullptr : &*it);
    }
    
  private:
    cp_type *              _max_at_inf;
    cp_container                  _cps;
    std::vector<cp_type *>     _minima;
};

template <typename number_type, typename size_type>
bool validate(flow_complex<number_type, size_type> const& fc) {
  std::map<size_type, int> hist;
  for (auto const& cp : fc)
    if (not cp.is_max_at_inf()) {
      auto r_pair = hist.emplace(cp.index(), 1);
      if (not r_pair.second)
        ++r_pair.first->second;
    }
  int sum = 0;
  for (auto const& el : hist)
    sum += (0 == (el.first % 2) ? el.second : -el.second);
  return 1 == sum;
}

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

