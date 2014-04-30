#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <algorithm>
#include <vector>
#include <utility>

#include "critical_point.hpp"

namespace FC {

template <typename _number_type, typename _size_type>
class flow_complex {
  public:
    typedef _number_type                             number_type;
    typedef _size_type                                 size_type;
    typedef critical_point<_number_type, _size_type>     cp_type;

  private:
    using cp_container = std::vector<cp_type>;
    
  public:
    typedef typename cp_container::iterator             iterator;
    typedef typename cp_container::const_iterator const_iterator;
  
    /**
      @brief initializes the flow complex with a maximum at infinity
    */
    flow_complex(size_type dim) {
      _cps.emplace_back(dim);
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
      auto it = std::find(_cps.begin(), _cps.end(), cp);
      if (_cps.end() == it) {
        // insert the new cp
        _cps.push_back(std::move(cp));
        r.first  = true;
        r.second = &_cps.back();
      } else {
        r.first = false;
        r.second = &*it;
      }
      return r;
    }
    
  private:  
    cp_container _cps;
};

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

