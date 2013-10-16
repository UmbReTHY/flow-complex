#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <utility>
#include <vector>

#include <tbb/concurrent_unordered_set.h>

#include "point_cloud.hpp"
#include "critical_point.hpp"
#include "flow_point.hpp"

namespace fc {

template <typename _number_type, typename _size_type>
class flow_complex {
  public:
    typedef _number_type number_type;
    typedef _size_type     size_type;
  
  private:
    using cp_t = critical_point<_number_type, _size_type>;
    using cp_container = tbb::concurrent_unordered_set<cp_t>;
    using flow_point_t = flow_point<number_type, size_type>;

  public:
    template <class point_iterator>
    flow_complex(dim_type dim, point_iterator begin, point_iterator end)
      : _point_cloud(dim, begin, end), _critical_points(/*TODO insert max at inf*/),
        _cp_at_inf(/*TODO make this the ref of the first inserted cp*/) {
      auto const first_max = find_first_max();
      // TODO insert this first max as cp if not at inf into container
    }

  private:
    flow_point_t find_first_max() {
      std::vector<number_type> loc(_point_cloud.dim());
      gen_seed_point(_point_cloud.dim(),
                     _point_cloud.begin(), _point_cloud.end(),
                     loc.data());  // TODO
      flow_point_t fp(_point_cloud, loc.data());
      
      while (!fp.is_proxy_at_inf() && !cp.is_finite_max())
        ascend();
        
      return std::move(fp);
    }
    
    cp_t const& cp_at_inf() const {
      return _cp_at_inf;
    }
  
    point_cloud<number_type, dim_type, size_type> const _point_cloud;
    cp_container _critical_points;
    cp_t const _cp_at_inf;  // TODO make this a ref
};

}  // namespace fc

#endif  // FLOW_COMPLEX_HPP_

