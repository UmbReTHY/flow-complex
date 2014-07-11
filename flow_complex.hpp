#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <algorithm>
#include <utility>

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
    iterator begin() noexcept {
      return _cps.begin();
    }
    
    iterator end() noexcept {
      return _cps.end();
    }
    
    // const-iterators
    const_iterator cbegin() const {
      return _cps.cbegin();
    }
    
    const_iterator cend() const {
      return _cps.cend();
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

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

