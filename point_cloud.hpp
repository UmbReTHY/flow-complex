#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <cassert>

#include <iterator>
#include <vector>

namespace FC {

/**
  @brief lightweight wrapper class around existing data
         An instance of this class does not claim ownership of the data that it
         has been given a pointer to.
         It also does not remain the size of the point. Hence,
         out-of-bounce-indexing is not checked.
*/
template <typename _number_type>
class point {
  public:
    typedef _number_type number_type;

    point(number_type const* data) : _data(data) {
      assert(data && "POINTER 'data' IS NULL-POINTER");
    }

    number_type const* data() const {
      return _data;
    }

    template <typename dim_type>
    number_type operator[](dim_type idx) const noexcept {
      return _data[idx];
    }

  private:
    number_type const* _data;
};

template <typename _number_type, typename _dim_type, typename _size_type>
class point_cloud {
  public:
    typedef _number_type number_type;
    typedef _dim_type dim_type;
    typedef _size_type size_type;
    typedef point<number_type> point_type;

  private:
    using point_container = std::vector<point_type>;    

  public:
    typedef typename point_container::const_iterator const_iterator;

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
      _points.reserve(std::distance(begin, end));
      for (;begin != end; ++begin)
        _points.push_back(point_type(&(*begin)[0]));
    }

    point_cloud() = delete;
    point_cloud(point_cloud const&) = delete;
    point_cloud & operator=(point_cloud const&) = delete;
    point_cloud(point_cloud const&&) = delete;
    point_cloud & operator=(point_cloud const&&) = delete;

    const_iterator begin() const {
      return _points.cbegin();
    }

    const_iterator end() const {
      return _points.cend();
    }

    /**
      @param idx idx of the point within the point cloud
    */
    point_type const& operator[](size_type const idx) const noexcept {
      return _points[idx];
    }

    dim_type dim() const noexcept {
      return _dim;
    }

    size_type size() const {
      return _points.size();
    }

  private:

    dim_type const _dim;
    point_container _points;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

