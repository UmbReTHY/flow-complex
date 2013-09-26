#ifndef AFFINE_HULL_HPP_
#define AFFINE_HULL_HPP_

#include <cassert>

#include <algorithm>
#include <vector>

#include "point_cloud.hpp"
#include "dynamic_qr.hpp"

namespace fc {

template <typename _number_type, typename _size_type, typename _dim_type>
class affine_hull {
  using member_container = std::vector<_size_type>;

  public:  
    typedef _size_type size_type;
    typedef member_container::const_iterator const_iterator;
    typedef point_cloud<number_type, dim_type, size_type> point_cloud_t;
    
    affine_hull(point_cloud_t const& pc) : _pc(pc) {
    }
    
    void add_point(size_type const idx) {
      assert(!is_member(idx));
      _members.push_back(idx);
    }
    
    void drop_point(size_type const idx) {
      assert(is_member(idx));
      auto it = std::find(_members.begin(), _members.end(), idx);
      _members.erase(it);
    }
    
    /**
      @return the number of points 
    */
    size_type size() const {
      return _members.size();
    }
    
    const_iterator begin() const {
      return _members.begin();
    }
    
    const_iterator end() const {
      return _members.end();
    }
  
  private:
    bool is_member(_size_type const idx) const {
      return _members.end() != std::find(_members.begin(), _members.end(), idx);
    }

    point_cloud_t const& _pc;  
    std::vector<_size_type> _members;
    dynamic_qr<_number_type> _dyn_qr;  // TODO: here
};

}  // namespace fc

#endif  // AFFINE_HULL_HPP_

