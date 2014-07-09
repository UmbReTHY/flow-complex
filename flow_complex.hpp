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
    flow_complex(size_type dim) {
      _max_at_inf = &*_cps.insert(cp_type(dim)).first;
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
    
  private:
    cp_type *    _max_at_inf;
    cp_container        _cps;
};

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

