//+++HDR+++
//======================================================================
//  This file is part of the SL software library.
//
//  Copyright (C) 1993-2014 by Enrico Gobbetti (gobbetti@crs4.it)
//  Copyright (C) 1996-2014 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://www.crs4.it/vvr/.
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//
//======================================================================
//---HDR---//
#ifndef SL_AXIS_ALIGNED_BOX_HPP
#define SL_AXIS_ALIGNED_BOX_HPP

#include <sl/bounding_volume.hpp>
#include <sl/linear_map.hpp>
#include <sl/fixed_size_plane.hpp>
#include <sl/interval.hpp>

namespace sl {

  namespace detail {
    template <size_t N_dimension, class T> 
    struct axis_aligned_box_data {
      typedef fixed_size_point<N_dimension,T> point_t;
      
      point_t limits_[2];

      inline axis_aligned_box_data() {};
      inline axis_aligned_box_data(const point_t& p0,
                                   const point_t& p1) {
        limits_[0] = p0;
        limits_[1] = p1;
      }
    public: // Serialization
    
      void store_to(output_serializer& s) const {
        s << limits_[0] << limits_[1];
      }
    
      void retrieve_from(input_serializer& s) {
        s >> limits_[0] >> limits_[1];
      }

    };
  }
      
  /**
   * An axis-aligned N-dimensional box
   */
  template <size_t N_dimension, class T> 
  class axis_aligned_box: 
    public 
      bounding_volume_base<
        N_dimension,
        T,
        axis_aligned_box<N_dimension,T>,
        detail::axis_aligned_box_data<N_dimension,T> 
      > 
  {
  public:
    enum { dimension = N_dimension };

    typedef bounding_volume_base<
        N_dimension,
        T,
        axis_aligned_box<N_dimension,T>,
        detail::axis_aligned_box_data<N_dimension,T> 
    > super_t;

    typedef typename super_t::derived_t         this_t;
    typedef typename super_t::value_t           value_t;
    typedef typename super_t::data_t            data_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::line_t            line_t;
    typedef typename super_t::ray_t             ray_t;
    typedef typename super_t::line_segment_t    line_segment_t;

    typedef fixed_size_vector<column_orientation,N_dimension,T> vector_t;
    typedef fixed_size_plane<N_dimension,T>                     plane_t;

  public: // Creation & Destruction

    /// Default init (positioned at the origin, zero size)
    inline axis_aligned_box() { 
    }

    /// Init from params
    inline axis_aligned_box(const point_t& pmin,
			    const point_t& pmax): super_t(data_t(pmin,pmax)) {
    }

    /// The box containing the single point p
    inline explicit axis_aligned_box(const point_t& p): super_t(data_t(p,p)) {
    }

  public: // Access

    /// the i-th extremum corner (0: minimum, 1: maximum)
    inline point_t& operator[](int i) { 
      return (this->data_).limits_[i]; 
    }
	
    /// the i-th extremum corner (0: minimum, 1: maximum)
    inline const point_t& operator[](int i) const { 
      return (this->data_).limits_[i]; 
    }

  public: // Queries
		
    /// True iff the box is empty (i.e. one of the lengths is < 0)
    inline bool is_empty() const {
	bool result = true;
	for (size_t i=0; i<dimension && result; i++) {
	  result = (*this)[0][i] > (*this)[1][i];
	}
	return result;
    }

    /// The main diagonal of the box
    inline vector_t diagonal() const {
      return is_empty() ? vector_t() : (*this)[1] - (*this)[0];
    }
	
    /// The half side lengths of the box
    inline vector_t half_side_lengths() const {
      const value_t one_half = reciprocal(scalar_math<value_t>::two());
      return is_empty() ? vector_t() : one_half * ((*this)[1] - (*this)[0]);
    }

    /// The center of the box
    inline point_t center() const {
      const value_t one_half = reciprocal(scalar_math<value_t>::two());
      return is_empty() ? point_t() : (*this)[0].lerp((*this)[1], one_half);
    }
      
    /// The number of vertices of the box
    inline size_t corner_count() const {
      return scalar_math<size_t>::one()<<dimension;
    }

    /// The vertex number i
    inline point_t corner(size_t i) const {
      SL_REQUIRE("Good index", i < corner_count());
      
      point_t result = tags::not_initialized();

      for (size_t d=0; d<dimension; d++) {
	if (i & (scalar_math<size_t>::one()<<d)) {
	  result[d] = (*this)[1][d];
	} else {
	  result[d] = (*this)[0][d];
	}
      }

      return result;
    }
      
    /// The i-th coordinate axis
    inline vector_t axis(size_t i) const {
      return vector_t::unit(i);
    }
    
    inline this_t inf_box(std::size_t k,
                          float d) const {
      point_t plo = (*this)[0];
      point_t phi = (*this)[1];
      phi[k] = d;
      return this_t(plo, phi);
    }

    inline this_t sup_box(std::size_t k,
                          float d) const {
      point_t plo = (*this)[0];
      point_t phi = (*this)[1];
      plo[k] = d;
      return this_t(plo, phi);
    }

    /** 
     *  The plane corresponding to the i-th face of the box.
     *  The numbering of the faces is -x, x, -y, y, ...
     *  The normal is oriented towards the outside of the box.
     */
    inline plane_t plane(size_t i) const {
      SL_REQUIRE("Good face index", i<2*dimension);
      size_t coord  = i/2;
      bool positive = (i%2) == 1;
      vector_t n = axis(coord);
      if (!positive) n = -n;
      point_t p = center() + half_side_lengths()[coord] * n;
      return plane_t(as_dual(n), p);
    }
 
    /// the volume of the box
    value_t volume() const {
      value_t result = scalar_math<value_t>::zero();
      if (!is_empty()) {
	result = scalar_math<value_t>::one();
	vector_t sz = diagonal();
	for (size_t i=0; i<dimension; i++) {
	  result *= sz[i];
	}
      }
      return result;
    }

    /// The side area of the face normal to the coordinate axis number i
    value_t side_area(size_t i) const {
      SL_REQUIRE("Good dimension", i<dimension);
      value_t result = scalar_math<value_t>::one(); 
      vector_t d = diagonal();
      for (size_t j=0; j<dimension; j++) {
	if (j!=i) {
	  result *= d[j];
	}
      }
      return result;
    }

    /// The projected area of the box in direction d
    value_t projected_area(const vector_t& d) const {
      value_t result = scalar_math<value_t>::one(); 
      for (size_t j=0; j<dimension; j++) {
	result += sl::abs(d[j]) * side_area(j);
      }
      return result;
    }

    /// The square of the minimum Euclidean distance to the point p
    value_t square_distance_to(const point_t& p) const {
      value_t result = scalar_math<value_t>::zero();
      for (size_t j=0; j<dimension; ++j) {
	if (p[j] < (*this)[0][j] ) {
	  result += sl::sqr(p[j] - (*this)[0][j]);
	} else if (p[j] > (*this)[1][j]) {
	  result += sl::sqr(p[j] - (*this)[1][j]);
	}
      }
      return result;
    }

    /// The minimum Euclidean distance to the point p
    inline value_t distance_to(const point_t& p) const {
      return std::sqrt(square_distance_to(p));
    }

    /// The range of signed distances to the plane hp
    inline interval<value_t> signed_distance_range(const plane_t& hp) const {
      point_t p_hi = tags::not_initialized();
      point_t p_lo = tags::not_initialized();

      for (size_t i=0; i<dimension; ++i) {
	if (sl::is_positive(hp[i])) {
	  p_lo[i] = (*this)[0][i];
	  p_hi[i] = (*this)[1][i];
	} else {
	  p_lo[i] = (*this)[1][i];
	  p_hi[i] = (*this)[0][i];
	}
      }

      return interval<value_t>(hp.value(p_lo), hp.value(p_hi));
    }

    /// Is the box fully on the side of the plane indicated by its normal? 
    inline bool is_fully_above(const plane_t& hp) const {
      point_t p_near = tags::not_initialized();
      for (size_t i=0; i<dimension; ++i) {
	p_near[i] = (*this)[sl::is_positive(hp[i]) ? 0 : 1][i];
      }
      return sl::is_positive(hp.value(p_near));
    }

    /// Is the box fully on the opposite side of the plane indicated by its normal? 
    inline bool is_fully_below(const plane_t& hp) const {
      point_t p_far = tags::not_initialized();
      for (size_t i=0; i<dimension; ++i) {
	p_far[i] = (*this)[sl::is_positive(hp[i]) ? 1 : 0][i];
      }
      return sl::is_non_positive(hp.value(p_far));
    }

  public: // Containment 

    /// True iff p (in global coordinates) is inside this
    bool contains(const point_t& p) const {
      bool result = true;
      for (size_t i=0; i<dimension && result; i++) {
	result = ((*this)[0][i]<=p[i]) && (p[i]<=(*this)[1][i]);
      }
      return result;
    }
	
  public: // Ray-tracing

    /**
     * Clip a linear component with this axis aligned box.
     * on entry, tmin and tmax determine the extent of
     * the linear component. On exit, off is true if
     * the linear component is outside of the box, and
     * tmin, tmax contain the portion of the component
     * inside the box.
     */
    void clip(const point_t& org,
	      const vector_t& dir,
	      bool&    off,
	      value_t& tmin,
	      value_t& tmax) const {
      off = (tmin > tmax);
      for (size_t i=0; (i< dimension) && (!off); ++i) {
	const value_t min_i = (*this)[0][i];
	const value_t max_i = (*this)[1][i];
	const value_t org_i = org[i];
	const value_t dir_i = dir[i];
	
#if 1
	if (is_zero(dir_i)) {
	  off = (org_i < min_i || org_i > max_i);
	} else {
	  const value_t inv_dir_i = reciprocal(dir_i);
	  const value_t t1 = (min_i - org_i) * inv_dir_i;
	  const value_t t2 = (max_i - org_i) * inv_dir_i;
	  if (t1 < t2) {
	    if (t1 > tmin) tmin = t1;
	    if (t2 < tmax) tmax = t2;
	  } else {
	    if (t2 > tmin) tmin = t2;
	    if (t1 < tmax) tmax = t1;
	  }
	  off = (tmin > tmax);
	}
#else
	if (is_negative(dir_i)) {
	  const value_t inv_dir_i = reciprocal(dir_i);
	  value_t t = (min_i - org_i) * inv_dir_i;
	  off = (t < tmin);
	  if (!off) {
	    if (t < tmax) {
	      tmax = t;
	    }
	    t = (max_i - org_i) * inv_dir_i;
	    if (t >= tmin) {
	      off = (t > tmax);
	      if (!off) {
		tmin = t;
	      }
	    }
	  }
	} else if (is_positive(dir_i)) {
	  const value_t inv_dir_i = reciprocal(dir_i);
	  value_t t = (max_i - org_i) * inv_dir_i;
	  off = (t < tmin);
	  if (!off) {
	    if (t < tmax) {
	      tmax = t;
	    }
	    t = (min_i - org_i) * inv_dir_i;
	    if (t >= tmin) {
	      off = (t > tmax);
	      if (!off) {
		tmin = t;
	      }
	    }
	  }	    
	} else { // is_zero(dir_i)
	  off = (org_i < min_i || org_i > max_i);
	}
#endif
      }
    }

    /// True iff the ray intersects the volume
    inline bool intersection_exists(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }


    /// True iff the line intersects the volume
    inline bool intersection_exists(const line_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }

    /// True iff the line segment intersects the volume
    inline bool intersection_exists(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }

    /**
     * ray - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest ray parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      if (off) {
	return std::pair<bool,value_t>(false,tmin);
      } else if (is_zero(tmin)) {
	  // No` near intersection farther than epsilon away, look at far intersection
	return std::pair<bool,value_t>(true,tmax);
      } else {
	return std::pair<bool,value_t>(true,tmin);
      }
    }

    /**
     * line - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest line parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const line_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

    /**
     * segment - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      if (off) {
	return std::pair<bool,value_t>(false,tmin);
      } else if (is_zero(tmin)) {
	return std::pair<bool,value_t>(!is_one(tmax),tmax);
      } else {
	return std::pair<bool,value_t>(true,tmin);
      }
    }

    /**
     * ray - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest ray parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

    /**
     * line - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest line parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const line_t& r) const {
      return intersection(r);
    }

    /**
     * segment - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

  public: // Commands

    /// Set this to the zero-sized box containing p
    inline void to(const point_t& p) {
      (*this)[0] = p;
      (*this)[1] = p;
    }

    /// Set this to the empty box
    inline void to_empty() {
      (*this)[0].fill(scalar_math<value_t>::upper_bound());
      (*this)[1].fill(scalar_math<value_t>::lower_bound());
    }

    /// Set this to the zero-sized box containing the origin
    inline void to_zero() {
      to(point_t());
    }
	
    /// Set this to the largest representable box
    inline void to_huge() {
      (*this)[0].fill(scalar_math<value_t>::lower_bound());
      (*this)[1].fill(scalar_math<value_t>::upper_bound());
    }

    /// Add point p to this, eventually adjusting the size of this
    void merge(const point_t& p) {
      for (size_t i =0; i< dimension; i++) {
	(*this)[0][i]=min((*this)[0][i],p[i]);
	(*this)[1][i]=max((*this)[1][i],p[i]);
      }
    }
      
    /// Merge box bb with this, eventually adjusting the size of this
    void merge(const this_t& bb) {
      for (size_t i =0; i< dimension; i++) {
	(*this)[0][i]=min((*this)[0][i],bb[0][i]);
	(*this)[1][i]=max((*this)[1][i],bb[1][i]);
      }
    }

  public: // Transformations

    /// the axis-aligned box containing this transformed by map 
    template <class T_tag>
    this_t transformed_by(const sl::linear_map_base<T_tag, dimension, value_t>& m) const {
      this_t result(m * corner(0));
      for (size_t i=1; i<corner_count(); i++) {
	result.merge(m * corner(i));
      }

      return result;
    }
    
  }; // axis_aligned_box

  template <std::size_t DIM, typename OUT_ET>
  class conv_to< axis_aligned_box<DIM, OUT_ET> > {
  public:
    typedef axis_aligned_box<DIM, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const axis_aligned_box<DIM, IN_ET>& in) {
      return result_t(conv_to< fixed_size_point<DIM,OUT_ET> >::from(in[0]),
		      conv_to< fixed_size_point<DIM,OUT_ET> >::from(in[1]));
    }
  };

}; // namespace sl

/// the axis_aligned box b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::axis_aligned_box<N_dimension,T> transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
						      const sl::axis_aligned_box<N_dimension,T>& b) {
  return b.transformed_by(map);
}

/// the axis_aligned box b transformed by the inverse of map
template <class T_tag, size_t N_dimension, class T>
inline sl::axis_aligned_box<N_dimension,T> inverse_transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
							      const sl::axis_aligned_box<N_dimension,T>& b) {
  return b.transformed_by(map.inverse());
}

/// the axis-aligned box b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::axis_aligned_box<N_dimension,T> operator * (const sl::linear_map_base<T_tag, N_dimension, T>& map, 
						     const sl::axis_aligned_box<N_dimension,T>& b) {
  return b.transformed_by(map);
}

// I/O

template <size_t N_dimension, class T>
std::ostream& operator <<(std::ostream& s, const sl::axis_aligned_box<N_dimension,T>& b) {
  s << b[0] << std::endl;;
  s << b[1] << std::endl;;
  return s;
}
    
template <size_t N_dimension, class T>
std::istream& operator >>(std::istream& s, sl::axis_aligned_box<N_dimension,T>& b) {
  typename sl::axis_aligned_box<N_dimension,T>::point_t p0;
  typename sl::axis_aligned_box<N_dimension,T>::point_t p1;

  s >> p0;
  s >> p1;

  b.to(p0);
  b.merge(p1);

  return s;
}

namespace sl {


  /// An object for building axis-aligned boxes from point clouds
  template <size_t N_dimension, class T> 
  class axis_aligned_box_builder: 
    public bounding_volume_builder<axis_aligned_box<N_dimension,T> > {
  public:

    typedef bounding_volume_builder<axis_aligned_box<N_dimension,T> > super_t;
    typedef axis_aligned_box_builder<N_dimension, T> this_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::value_t           value_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::bounding_volume_t bounding_volume_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  public: // Creation & Destruction

    /// Default init
    inline axis_aligned_box_builder() { }

    inline virtual ~axis_aligned_box_builder() { }

  public: // Construction

    /// Start point cloud
    virtual inline void begin_model() {
      super_t::begin_model();
      (this->last_bounding_volume_).to_empty();
    }

    /// Add point to current cloud
    virtual inline void put_point(const point_t& p) {
      (this->last_bounding_volume_).merge(p);
    }

    /// Finish point cloud and build box
    virtual inline void end_model() {
      super_t::end_model();
    }

  }; // class axis_aligned_box_builder

}; // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  // A 2D axis-aligned box with single precision floating point components
  typedef axis_aligned_box<2,float> aabox2f;
  // A 3D axis-aligned box with single precision floating point components
  typedef axis_aligned_box<3,float> aabox3f;
  // A 4D axis-aligned box with single precision floating point components
  typedef axis_aligned_box<4,float> aabox4f;

  // A 2D axis-aligned box with double precision floating point components
  typedef axis_aligned_box<2,double> aabox2d;
  // A 3D axis-aligned box with double precision floating point components
  typedef axis_aligned_box<3,double> aabox3d;
  // A 4D axis-aligned box with double precision floating point components
 typedef axis_aligned_box<4,double> aabox4d;
}

#endif
