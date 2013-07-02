#ifndef FLOW_POINT_HPP_
#define FLOW_POINT_HPP_

#include <utility>

namespace FC {

template <typename number_type, typename dim_type>
class flow_point {
  public:
    flow_point(dim_type const dim,
               point_cloud const& P,
               std:scoped_ptr<numbertype[]> location);

  private:
    std:scoped_ptr<numbertype[]> _location;
};

}  // namespace FC

#endif  // FLOW_POINT_HPP_

