#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <nanoflann.hpp>

namespace FC {

/**
  @brief lightweight wrapper class around existing data
         An instance of this class does not claim ownership of the data that it
         has been given a pointer to.
         It also does not remain the size of the point.
*/
template <typename number_type, typename size_type>
class point {
  public:
    point(number_type * data) : _data(data) {}

    number_type & operator[](size_type idx) nothrow {
      return _data[idx];
    }

    number_type operator[](size_type idx) const nothrow {
      return _data[idx];
    }

  private:
    number_type * _data;
};

template <typename number_type, typename size_type, typename dim_type>
class point_cloud {
  public:
    typedef point<number_type, size_type> point_type;

    /**
      @tparam point_iterator the expression (*<point_iterator>)[i] must be valid
                             and the returned values are expected to be located
                             in contiguous memory
      @params begin,end iterators defining the range of points that constitute
                        the point cloud.
    */
    template <typename point_iterator>
    point_cloud(dim_type const dim, point_iterator begin, point_iterator end)
      : _dim(dim) {
      for (;begin != endl ++begin)
    }

    /**
      @param idx idx of the point within the point cloud
    */
    point_type operator[](size_type idx) nothrow {
    }



  private:
    dim_type const _dim;
    std::vector<number_type *> _points;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

