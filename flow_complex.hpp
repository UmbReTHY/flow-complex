#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <vector>

#include <tbb/concurrent_unordered_set.h>

#include "point_cloud.hpp"
#include "critical_point.hpp"
#include "flow_point.hpp"
#include "utility.hpp"

namespace FC {

template <typename _number_type, typename _size_type, typename _dim_type>
class flow_complex {
  public:
    typedef _number_type number_type;
    typedef _size_type     size_type;
    typedef _dim_type       dim_type;
  
  private:
    using cp_container = tbb::concurrent_unordered_set<critical_point>;  // TODO template parameters for cp
    using flow_point_t = flow_point<number_type, size_type, dim_type>;

  public:
    
    template <class point_iterator>
    flow_complex(dim_type dim, point_iterator begin, point_iterator end)
      : (dim, begin, end) {
    }



  private:
    flow_point_t find_first_max() {
      std::vector<number_type> loc(_point_cloud.dim());
      gen_seed_point(_point_cloud.dim(),
                     _point_cloud.begin(), _point_cloud.end(),
                     loc.data());
      flow_point_t fp(_point_cloud, loc.data());
      
      while () { // TODO: complete loop
        ascend();
      }
    }
  
    point_cloud<number_type, dim_type, size_type> const _point_cloud;
    cp_container _critical_points;
};

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

