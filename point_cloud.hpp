#ifndef POINT_CLOUD_HPP_
#define POINT_CLOUD_HPP_

#include <iterator>
#include <vector>

#include <Eigen/Dense>
#include <nanoflann.hpp>

namespace FC {

/**
  @brief lightweight wrapper class around existing data
         An instance of this class does not claim ownership of the data that it
         has been given a pointer to.
         It also does not remain the size of the point. Hence,
         out-of-bounce-indexing is not checked.
*/
template <typename number_type, typename dim_type>
class point {
  public:
    point(number_type const* data) : _data(data) {}

    number_type operator[](dim_type idx) const noexcept {
      return _data[idx];
    }

  private:
    number_type const* _data;
};

template <typename _number_type, typename _dim_type>
class point_cloud {
  public:
    typedef _number_type number_type;
    typedef _dim_type dim_type;
    typedef point<number_type, dim_type> point_type;

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
      : _dim(dim), _pc_adaptor(*this), _kd_tree(dim, _pc_adaptor) {  // TODO(Lars): consider using params
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
    template <typename size_type>
    point_type const& operator[](size_type idx) const noexcept {
      return _points[idx];
    }

    dim_type dim() const noexcept {
      return _dim;
    }

    std::size_t size() const {
      return _points.size();
    }

  private:
    /**
      @brief this class
    */
    class point_cloud_adaptor {
      public:
        point_cloud_adaptor(point_cloud const& pc) : _pc(pc) {
        }

        std::size_t kdtree_get_point_count() const {
          return _pc.size();
        }

        number_type kdtree_distance(number_type const* p1,
                                    std::size_t const idx_p2,
                                    std::size_t const dim) const {
          using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
          return (Eigen::Map<eigen_vector>(p1, dim) -
                  Eigen::Map<eigen_vector>(&_pc[idx_p2][0], dim)).squaredNorm();
        }

        number_type kdtree_get_pt(std::size_t const idx, int dim) const {
          return _pc[idx][dim];
        }

      private:
        point_cloud const& _pc;
    };

    dim_type const _dim;
    point_container _points;
    // members for nearest neighbor queries
    point_cloud_adaptor _pc_adaptor;
    nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<number_type, point_cloud_adaptor>,
      point_cloud_adaptor> _kd_tree;
};

}  // namespace FC

#endif  // POINT_CLOUD_HPP_

