#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <vector>

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
  
    flow_complex() = default;
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
    
  private:  
    cp_container _cps;
};

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

