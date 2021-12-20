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
#ifndef SL_ORIENTED_BOX_HPP
#define SL_ORIENTED_BOX_HPP

#include <sl/bounding_volume.hpp>
#include <sl/rigid_body_map.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/linear_map_factory.hpp>
#include <sl/convex_hull.hpp>
#include <vector>

#define BUGGY_ROTATING_CALIPER 1

namespace sl {

  namespace detail {
    template <size_t N_dimension, class T> 
    struct oriented_box_data {
      typedef rigid_body_map<N_dimension,T>                       rigid_body_map_t;
      typedef fixed_size_vector<column_orientation,N_dimension,T> vector_t;
      
      rigid_body_map_t                    from_box_space_map_;
      vector_t                            half_side_length_;

      inline oriented_box_data() {}
      inline oriented_box_data(const rigid_body_map_t& rbm,
			       const vector_t& hsl) 
	: from_box_space_map_(rbm), half_side_length_(hsl) {
      }
    public: // Serialization
    
      void store_to(output_serializer& s) const {
        s << from_box_space_map_ << half_side_length_;
      }
      
      void retrieve_from(input_serializer& s) {
        s >> from_box_space_map_ >> half_side_length_;
      }
    };
  }

  /**
   * An oriented N-dimensional box
   */
  template <size_t N_dimension, class T> 
  class oriented_box :
    public 
      bounding_volume_base<
        N_dimension,
        T,
        oriented_box<N_dimension,T>,
        detail::oriented_box_data<N_dimension,T> 
      > 
  {
  public:
    enum { dimension = N_dimension };

    typedef bounding_volume_base<
        N_dimension,
        T,
        oriented_box<N_dimension,T>,
        detail::oriented_box_data<N_dimension,T> 
    > super_t;

    typedef typename super_t::derived_t         this_t;
    typedef typename super_t::value_t           value_t;
    typedef typename super_t::data_t            data_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::line_t            line_t;
    typedef typename super_t::ray_t             ray_t;
    typedef typename super_t::line_segment_t    line_segment_t;

    typedef fixed_size_vector<column_orientation,N_dimension,T> vector_t;
    typedef fixed_size_vector<row_orientation,dimension,T>      dual_vector_t;
    typedef fixed_size_plane<N_dimension,T>                     plane_t;
    typedef axis_aligned_box<dimension,T>                       axis_aligned_box_t;
    typedef rigid_body_map<dimension,T>                         rigid_body_map_t;
    typedef linear_map_factory<dimension,value_t>               linear_map_factory_t;

  public: // Creation & Destruction

    /// Default creation (identity position, zero-size)
    inline oriented_box() { 
    }

    /// Explicit creation 
    inline oriented_box(const rigid_body_map_t& from_box_space_map,
			const vector_t& half_side_length):
      super_t(data_t(from_box_space_map,half_side_length)) {
    }

    /// Init from axis-aligned box
    inline oriented_box(const axis_aligned_box_t& aab)
      :
      super_t(data_t(linear_map_factory_t::translation(as_vector(aab.center())), aab.half_side_lengths())) {
    }

  public: // Access

    /// The map transforming local to global coordinates
    inline const rigid_body_map_t& from_box_space_map() const { 
      return (this->data_).from_box_space_map_; 
    }
	
    /// The map transforming global to local coordinates
    inline rigid_body_map_t to_box_space_map() const { 
      return (this->data_).from_box_space_map_.inverse(); 
    }

    /// The half side lengths of the box
    inline const vector_t& half_side_lengths() const { 
      return (this->data_).half_side_length_;
    }
      
  public: // Transformations of points/vectors

    /// Point p transformed to local coordinates
    inline point_t to_box_space(const point_t& p) const {
      return inverse_transformation(from_box_space_map(), p);
    }
      
    /// Vector v transformed to local coordinates
    inline vector_t to_box_space(const vector_t& v) const {
      return inverse_transformation(from_box_space_map(), v);
    }

  public: // Queries
		
    /// True iff the box is empty (i.e. one of the lengths is < 0)
    bool is_empty() const {
	bool result = true;
	for (size_t i=0; i<dimension && result; ++i) {
	  result = is_negative((this->data_).half_side_length_[i]);
	}
	return result;
    }

    /// The number of vertices of the box
    inline size_t corner_count() const {
      return 1<<dimension;
    }


    /// The vertex number i, in local coordinates
    inline point_t box_space_corner(size_t i) const {
      SL_REQUIRE("Good index", i < corner_count());
      
      point_t result = tags::not_initialized();
      for (int d=0; d<dimension; ++d) {
	if (i & (1<<d)) {
	  result[d] = half_side_lengths()[d];
	} else {
	  result[d] = -half_side_lengths()[d];
	}
      }
      return result;
    }

    /// The vertex number i, in global coordinates.
    inline point_t corner(size_t i) const {
      return from_box_space_map() * box_space_corner(i);
    }

    /// The main diagonal of the box, in local coordinates.
    inline vector_t box_space_diagonal() const {
      return box_space_corner(corner_count()-1) - box_space_corner(0);
    }

    /// The main diagonal of the box, in global coordinates.
    inline vector_t diagonal() const {
      return from_box_space_map() * box_space_diagonal();
    }
	
    /// The center of the box, in global coordinates.
    inline point_t center() const {
      return from_box_space_map() * point_t();
    }
      
    /// The volume of the box ("area" in 2D)
    value_t volume() const {
      value_t result = scalar_math<value_t>::zero();
      if (!is_empty()) {
	result = scalar_math<value_t>::one();
	vector_t sz = diagonal();
	for (size_t i=0; i<dimension; ++i) {
	  result *= sz[i];
	}
      }
      return result;
    }

    /// The side area of the face normal to the coordinate axis number i
    value_t side_area(size_t i) const {
      SL_REQUIRE("Good dimension", i<dimension);
      value_t result = scalar_math<value_t>::one(); 
      for (size_t j=0; j<dimension; ++j) {
	if (j!=i) {
	  result *= (this->data_).half_side_length_[j] + (this->data_).half_side_length_[j];
	}
      }
      return result;
    }

    /// The projected area of the box in direction d
    value_t projected_area(const vector_t& d) const {
      vector_t local_d = to_box_space(d);

      value_t result = scalar_math<value_t>::zero(); 
      for (size_t j=0; j<dimension; ++j) {
	result += sl::abs(local_d[j]) * side_area(j);
      }
      return result;
    }

    /// The i-th coordinate axis
    vector_t axis(size_t i) const {
      SL_REQUIRE("Good dimension", i<dimension);
      vector_t result = tags::not_initialized();
      for (size_t j=0; j<dimension; ++j) {
	result[j] = (this->data_).from_box_space_map_(j,i);
      }
      return result;
    }

    /** 
     *  The plane corresponding to the i-th face of the box.
     *  The numbering of the faces is -x, x, -y, y, ...
     *  The normal is oriented towards the outside of the box.
     */
    plane_t plane(size_t i) const {
      SL_REQUIRE("Good face index", i<2*dimension);
      size_t coord  = i/2;
      bool positive = (i%2) == 1;
      vector_t n = axis(coord);
      if (!positive) n = -n;
      point_t p = center() + half_side_lengths()[coord] * n;
      return plane_t(as_dual(n), p);
    }
      

  public: // Containment 

    /// True iff p (in global coordinates) is inside this
    bool contains(const point_t& p) const {
      bool result = true;
      point_t localp = to_box_space(p);
      for (size_t i=0; i<dimension && result; ++i) {
	result = (-half_side_lengths()[i]<=localp[i]) && (localp[i]<=half_side_lengths()[i]);
      }
      return result;
    }
	
    /// True iff this is not touching other
    bool is_disjoint(const this_t& other) const {
      return obb_disjoint(to_box_space_map() * other.from_box_space_map(),
			  half_side_lengths(),
			  other.half_side_lengths());
    }

    /// True iff this does totally contain other
    bool contains(const this_t& other) const {
      bool result = true;
      for (size_t i=0; i< corner_count() && result; ++i) {
	result = contains(other.corner(i));
      }
      return result;
    }

    /// The square of the minimum Euclidean distance to the point p
    value_t square_distance_to(const point_t& p) const {
      point_t local_p =  inverse_transformation(from_box_space_map(), p);
      value_t result = sl::zero(value_t());
      for (size_t j=0; j<dimension; ++j) {
	if ( local_p[j] < -(this->data_).half_side_length_[j] ) {
	  result += sl::sqr(local_p[j] + (this->data_).half_side_length_[j]);
	} else if (local_p[j] > (this->data_).half_side_length_[j]) {
	  result += sl::sqr(local_p[j] - (this->data_).half_side_length_[j]);
	}
      }
      return result;
    }

    /// The minimum Euclidean distance to the point p
    inline value_t distance_to(const point_t& p) const {
      return std::sqrt(square_distance_to(p));
    }

    /// The range of signed distances to the plane hp
    interval<value_t> signed_distance_range(const plane_t& hp) const {
      plane_t local_hp = inverse_transformation(from_box_space_map(), hp);

      value_t v_lo = local_hp[dimension];
      value_t v_hi = local_hp[dimension];
      for (size_t j=0; j<dimension; ++j) {
	value_t local_hp_j_times_dj = local_hp[j] * (this->data_).half_side_length_[j];
	if (sl::is_positive(local_hp_j_times_dj)) {
	  v_lo -= local_hp_j_times_dj;
	  v_hi += local_hp_j_times_dj;
	} else {
	  v_lo += local_hp_j_times_dj;
	  v_hi -= local_hp_j_times_dj;
	}
      }

      return interval<value_t>(v_lo, v_hi);
    }
      
    /// The range of signed distances to the plane passing to the origin and orthogonal to the coordinate axis i
    interval<value_t> signed_distance_range(size_t i) const {
      SL_REQUIRE("Good coordinate", i < dimension);

      value_t v_lo = (this->data_).from_box_space_map_(i,dimension);
      value_t v_hi = (this->data_).from_box_space_map_(i,dimension);
      for (size_t j=0; j<dimension; ++j) {
	value_t local_hp_j_times_dj = (this->data_).from_box_space_map_(i,j) * (this->data_).half_side_length_[j];
	if (sl::is_positive(local_hp_j_times_dj)) {
	  v_lo -= local_hp_j_times_dj;
	  v_hi += local_hp_j_times_dj;
	} else {
	  v_lo += local_hp_j_times_dj;
	  v_hi -= local_hp_j_times_dj;
	}
      }
      return interval<value_t>(v_lo, v_hi);
    }

    /// The axis-aligned-box containing this
    axis_aligned_box_t bounding_axis_aligned_box() const {
      axis_aligned_box_t result;
      for (size_t i=0; i<dimension; ++i) {
	value_t v_lo = (this->data_).from_box_space_map_(i,dimension);
	value_t v_hi = (this->data_).from_box_space_map_(i,dimension);
        for (size_t j=0; j<dimension; ++j) {
	  value_t local_hp_j_times_dj = (this->data_).from_box_space_map_(i,j) * (this->data_).half_side_length_[j];
	  if (sl::is_positive(local_hp_j_times_dj)) {
	    v_lo -= local_hp_j_times_dj;
	    v_hi += local_hp_j_times_dj;
	  } else {
	    v_lo += local_hp_j_times_dj;
	    v_hi -= local_hp_j_times_dj;
	  }
	}
	result[0][i] = v_lo;
	result[1][i] = v_hi;
      }
      return result;
    }

    /// Is the box fully on the side of the plane indicated by its normal? 
    bool is_fully_above(const plane_t& hp) const {
      plane_t local_hp = inverse_transformation(from_box_space_map(), hp);
      point_t p_near = tags::not_initialized();
      for (size_t i=0; i<dimension; ++i) {
	p_near[i] = 
	  (sl::is_positive(local_hp[i]) ? 
	   -(this->data_).half_side_length_[i] : 
	   (this->data_).half_side_length_[i]);
      }
      return sl::is_positive(local_hp.value(p_near));
    }

    /// Is the box fully on the side of the plane indicated by its normal? 
    bool is_fully_below(const plane_t& hp) const {
      plane_t local_hp = inverse_transformation(from_box_space_map(), hp);
      point_t p_far = tags::not_initialized();
      for (size_t i=0; i<dimension; ++i) {
	p_far[i] =
	  (sl::is_positive(local_hp[i]) ? 
	   (this->data_).half_side_length_[i] : 
	   -(this->data_).half_side_length_[i]);
      }
      return sl::is_non_positive(local_hp.value(p_far));
    }
      
  public: // Ray-tracing

   /**
    * Clip a linear component with this oriented box.
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
      const typename rigid_body_map_t::matrix_t& m = (this->data_).from_box_space_map_.as_matrix();

      off = (tmin > tmax);
      for (size_t i=0; (i<dimension) && (!off); ++i) {
	const value_t min_i = -(this->data_).half_side_length_[i];
	const value_t max_i =  (this->data_).half_side_length_[i];

	// transform i-th component of origin/direction to local coordinates
	value_t org_i = m(0,i) * (org[0] - m(0,dimension));
	value_t dir_i = m(0,i) * (dir[0]);
	for (size_t j=1; j<dimension; ++j) {
	  org_i+= m(j,i) * (org[j] - m(j,dimension));
	  dir_i+= m(j,i) * (dir[j]);
	}

	// Clip 
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
    bool intersection_exists(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }


    /// True iff the line intersects the volume
    bool intersection_exists(const line_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }

    /// True iff the line segment intersects the volume
    bool intersection_exists(const line_segment_t& r) const {
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
    std::pair<bool,value_t> intersection(const ray_t& r) const {
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
    std::pair<bool,value_t> intersection(const line_t& r) const {
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
    std::pair<bool,value_t> intersection(const line_segment_t& r) const {
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
    std::pair<bool,value_t> bounds_intersection(const ray_t& r) const {
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
    std::pair<bool,value_t> bounds_intersection(const line_t& r) const {
      return intersection(r);
    }

    /**
     * segment - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    std::pair<bool,value_t> bounds_intersection(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

  public: // Transformations

    /// This transformed by map
    template <class T_tag>
    this_t transformed_by(const sl::linear_map_base<T_tag, dimension, value_t>& m) const {
      this_t result;

      if (T_tag::is_rigid_body) {
	result.data_.half_side_length_   = (this->data_).half_side_length_;
	result.data_.from_box_space_map_ = m * result.data_.from_box_space_map_;
      } else {
	// Non rigid body xform - we actually build an axis_aligned box!
	axis_aligned_box_t aab(m * corner(0));
	for (size_t i=1; i<corner_count(); ++i) {
	  aab.merge(m * corner(i));
	}
	result = aab;
      }
	
      return result;
    }

  public: // Helpers

    /**
     * Overlap test for 3D boxes
     *
     * from S. Gottschalk, M. Lin and D. Manocha (1996) 
     * "OBB-Tree: A Hierarchical Structure for Rapid Interference Detection", Proc. SIGGRAPH.
     * This is a test between two boxes, box A and box B.  It is assumed that
     * the coordinate system is aligned and centered on box A.  
     * ``BtoA'' is the rigid_body map that specifies the position of B in the reference frame
     * of A, The dimensions of the boxes are specified in ``a'' and ``b''.
     * dimensions of box A are given in a.
     */
    static int obb_disjoint(const rigid_body_map_t& BtoA,
			    const vector_t& a, 
			    const vector_t& b) {
      SL_REQUIRE("Implementation limit", dimension == 3);
 
      const value_t reps = std::sqrt(scalar_math<value_t>::epsilon());
      int r;

      // T = translation
      value_t Tr[3];
      Tr[0] = BtoA(0,3);
      Tr[1] = BtoA(1,3);
      Tr[2] = BtoA(2,3);

      // Bf = fabs(BtoA) + eps
      value_t Bf[3][3];
      Bf[0][0] = sl::abs(BtoA(0, 0)) + reps;
      Bf[0][1] = sl::abs(BtoA(0, 1)) + reps;
      Bf[0][2] = sl::abs(BtoA(0, 2)) + reps;
      Bf[1][0] = sl::abs(BtoA(1, 0)) + reps;
      Bf[1][1] = sl::abs(BtoA(1, 1)) + reps;
      Bf[1][2] = sl::abs(BtoA(1, 2)) + reps;
      Bf[2][0] = sl::abs(BtoA(2, 0)) + reps;
      Bf[2][1] = sl::abs(BtoA(2, 1)) + reps;
      Bf[2][2] = sl::abs(BtoA(2, 2)) + reps;
  
      value_t t, s;
  
      // if any of these tests are one-sided, then the polyhedra are disjoint
      r = 1;
      
      // A1 x A2 = A0
      t = sl::abs(Tr[0]);
      
      r &= (t <= 
	    (a[0] + b[0] * Bf[0][0] + b[1] * Bf[1][0] + b[2] * Bf[2][0]));
      if (!r) return 1;
      
      // B1 x B2 = B0
      s = Tr[0]*BtoA(0,0) + Tr[1]*BtoA(0,1) + Tr[2]*BtoA(0,2);
      t = sl::abs(s);
      
      r &= ( t <=
	     (b[0] + a[0] * Bf[0][0] + a[1] * Bf[0][1] + a[2] * Bf[0][2]));
      if (!r) return 2;
      
      // A2 x A0 = A1
      t = sl::abs(Tr[1]);
      
      r &= ( t <= 
	     (a[1] + b[0] * Bf[0][1] + b[1] * Bf[1][1] + b[2] * Bf[2][1]));
      if (!r) return 3;
      
      // A0 x A1 = A2
      t = sl::abs(Tr[2]);
      
      r &= ( t <= 
	     (a[2] + b[0] * Bf[0][2] + b[1] * Bf[1][2] + b[2] * Bf[2][2]));
      if (!r) return 4;
      
      // B2 x B0 = B1
      s = Tr[0]*BtoA(1,0) + Tr[1]*BtoA(1,1) + Tr[2]*BtoA(1,2);
      t = sl::abs(s);
      
      r &= ( t <=
	     (b[1] + a[0] * Bf[1][0] + a[1] * Bf[1][1] + a[2] * Bf[1][2]));
      if (!r) return 5;
      
      // B0 x B1 = B2
      s = Tr[0]*BtoA(2,0) + Tr[1]*BtoA(2,1) + Tr[2]*BtoA(2,2);
      t = sl::abs(s);
      
      r &= ( t <=
	     (b[2] + a[0] * Bf[2][0] + a[1] * Bf[2][1] + a[2] * Bf[2][2]));
      if (!r) return 6;

      // A0 x B0
      s = Tr[2] * BtoA(0,1) - Tr[1] * BtoA(0,2);
      t = sl::abs(s);
      
      r &= ( t <= 
	     (a[1] * Bf[0][2] + a[2] * Bf[0][1] +
	      b[1] * Bf[2][0] + b[2] * Bf[1][0]));
      if (!r) return 7;
      
      // A0 x B1
      s = Tr[2] * BtoA(1,1) - Tr[1] * BtoA(1,2);
      t = sl::abs(s);

      r &= ( t <=
	     (a[1] * Bf[1][2] + a[2] * Bf[1][1] +
	      b[0] * Bf[2][0] + b[2] * Bf[0][0]));
      if (!r) return 8;
      
      // A0 x B2
      s = Tr[2] * BtoA(2,1) - Tr[1] * BtoA(2,2);
      t = sl::abs(s);
      
      r &= ( t <=
	     (a[1] * Bf[2][2] + a[2] * Bf[2][1] +
	      b[0] * Bf[1][0] + b[1] * Bf[0][0]));
      if (!r) return 9;

      // A1 x B0
      s = Tr[0] * BtoA(0,2) - Tr[2] * BtoA(0,0);
      t = sl::abs(s);

      r &= ( t <=
	     (a[0] * Bf[0][2] + a[2] * Bf[0][0] +
	      b[1] * Bf[2][1] + b[2] * Bf[1][1]));
      if (!r) return 10;

      // A1 x B1
      s = Tr[0] * BtoA(1,2) - Tr[2] * BtoA(1,0);
      t = sl::abs(s);

      r &= ( t <=
	     (a[0] * Bf[1][2] + a[2] * Bf[1][0] +
	      b[0] * Bf[2][1] + b[2] * Bf[0][1]));
      if (!r) return 11;

      // A1 x B2
      s = Tr[0] * BtoA(2,2) - Tr[2] * BtoA(2,0);
      t = sl::abs(s);

      r &= (t <=
	    (a[0] * Bf[2][2] + a[2] * Bf[2][0] +
	     b[0] * Bf[1][1] + b[1] * Bf[1][0]));
      if (!r) return 12;

      // A2 x B0
      s = Tr[1] * BtoA(0,0) - Tr[0] * BtoA(0,1);
      t = sl::abs(s);

      r &= (t <=
	    (a[0] * Bf[0][1] + a[1] * Bf[0][0] +
	     b[1] * Bf[2][2] + b[2] * Bf[1][2]));
      if (!r) return 13;

      // A2 x B1
      s = Tr[1] * BtoA(1,0) - Tr[0] * BtoA(1,1);
      t = sl::abs(s);

      r &= ( t <=
	     (a[0] * Bf[1][1] + a[1] * Bf[1][0] +
	      b[0] * Bf[2][2] + b[2] * Bf[0][2]));
      if (!r) return 14;

      // A2 x B2
      s = Tr[1] * BtoA(2,0) - Tr[0] * BtoA(2,1);
      t = sl::abs(s);

      r &= ( t <=
	     (a[0] * Bf[2][1] + a[1] * Bf[2][0] +
	      b[0] * Bf[1][2] + b[1] * Bf[0][2]));
      if (!r) return 15;
      
      return 0;  // should equal 0
    }

  }; // oriented_box

  template <std::size_t DIM, typename OUT_ET>
  class conv_to< oriented_box<DIM, OUT_ET> > {
  public:
    typedef oriented_box<DIM, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const oriented_box<DIM, IN_ET>& in) {
      return result_t(conv_to< typename oriented_box<DIM, IN_ET>::rigid_body_map_t >::from(in.from_box_space_map()),
		      conv_to< typename oriented_box<DIM, IN_ET>::vector_t >::from(in.half_side_lengths()));
    }
  };

} // namespace sl


/// the oriented box b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::oriented_box<N_dimension,T> transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
						      const sl::oriented_box<N_dimension,T>& b) {
  return b.transformed_by(map);
}

/// the oriented box b transformed by the inverse of map
template <class T_tag, size_t N_dimension, class T>
inline sl::oriented_box<N_dimension,T> inverse_transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
							      const sl::oriented_box<N_dimension,T>& b) {
  return b.transformed_by(map.inverse());
}

/// the oriented box b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::oriented_box<N_dimension,T> operator * (const sl::linear_map_base<T_tag, N_dimension, T>& map, 
						 const sl::oriented_box<N_dimension,T>& b) {
  return transformation(map,b);
}

// I/O

template <size_t N_dimension, class T>
std::ostream& operator <<(std::ostream& s, const sl::oriented_box<N_dimension,T>& b) {
  s << b.from_box_space_map() << std::endl;;
  s << b.half_side_lengths() << std::endl;;
  return s;
}
    
template <size_t N_dimension, class T>
std::istream& operator >>(std::istream& s, sl::oriented_box<N_dimension,T>& b) {
  typename sl::oriented_box<N_dimension,T>::rigid_body_map_t m;
  typename sl::oriented_box<N_dimension,T>::vector_t     sz;

  s >> m;
  s >> sz;

  b = sl::oriented_box<N_dimension,T>(m,sz);

  return s;
}

// --------------------------------------------------------------
// sl::oriented_box_builder_base<N_dimension,T>
// --------------------------------------------------------------

namespace sl {

  /// Base class for objects constructing oriented bounding boxes from point clouds.
  template <size_t N_dimension, class T> 
  class oriented_box_builder_base: 
    public bounding_volume_builder<oriented_box<N_dimension,T> > {
  public:
    typedef bounding_volume_builder<oriented_box<N_dimension,T> > super_t;
    typedef oriented_box_builder_base<N_dimension, T> this_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::value_t                    value_t;
    typedef typename super_t::point_t                    point_t;
    typedef typename super_t::bounding_volume_t          bounding_volume_t;

    typedef typename bounding_volume_t::vector_t         vector_t;
    typedef typename bounding_volume_t::dual_vector_t    dual_vector_t;
    typedef typename bounding_volume_t::rigid_body_map_t rigid_body_map_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  protected: // State

    std::vector<point_t> points_;

    value_t hsl_eps_;

  public: // Creation & Destruction

    oriented_box_builder_base() : hsl_eps_(scalar_math<value_t>::epsilon()) { }

    virtual ~oriented_box_builder_base() {}

  public: // Construction

    /// Begin model
    virtual void begin_model() {
      super_t::begin_model();
      points_.clear();
    }
    
    /// Add point to current cloud
    virtual void put_point(const point_t& p) {
      points_.push_back(p);
    }

    /// End model, construct box
    virtual void end_model() {
      if (points_.size() == 0) {
	// Empty - build default box
	(this->last_bounding_volume_) = bounding_volume_t();
      } else if (points_.size() == 1) {
	// Single point - build zero-sized box centered on the point
	(this->last_bounding_volume_) = bounding_volume_t(linear_map_factory<dimension,value_t>::translation(points_[0].as_vector()),
						  vector_t());
      } else {
	// At least two points, construct box using generic algorithm
	construct();
      }
      expand();
      points_.clear();
      super_t::end_model();
    }

    void set_half_side_length_epsilon(value_t x) {
      hsl_eps_ = x;
    }

    value_t half_side_length_epsilon() const {
      return hsl_eps_;
    }

  protected: // Construct

    typedef fixed_size_square_matrix<dimension,value_t>             covariance_matrix_t;
    typedef fixed_size_square_matrix<dimension,value_t>             eigen_matrix_t;
    typedef fixed_size_vector<column_orientation,dimension,value_t> eigen_vector_t;

    /// Expand half side lengths of current box by half_side_length_epsilon()
    virtual void expand() {
      vector_t hsl_new = (this->last_bounding_volume_).half_side_lengths();
      for (size_t i=0; i<dimension; ++i) {
	hsl_new[i] += hsl_eps_;
      }
      (this->last_bounding_volume_) = bounding_volume_t((this->last_bounding_volume_).from_box_space_map(),
						hsl_new);
    }

    /// Compute eigen basis from current point cloud
    virtual rigid_body_map_t eigen_basis() const {
      SL_REQUIRE("At least one point", points_.size() > 0);
      // Compute mean
      const value_t       s = sl::reciprocal(value_t(points_.size()));
      point_t             mean = points_[0];
      for (size_t i=1; i<points_.size(); ++i) {
	mean += points_[i].as_vector();
      }
      mean.as_vector() *= s;

      // Compute covariance matrix;
      covariance_matrix_t cov;
      for (size_t i=0; i<points_.size(); ++i) {
	vector_t xp = points_[i] - mean;
	for (size_t d1=0; d1<dimension; ++d1) {
	  for (size_t d2=0; d2<dimension; ++d2) {
	    cov(d1,d2) += s * (xp[d1] * xp[d2]);
	  }
	}
      }

      // Compute eigenvectors 
      eigen_vector_t eigen_values;
      eigen_matrix_t eigen_vectors;

      bool ok = false;
      cov.symmetric_sorted_eigen_in(eigen_values,
				    eigen_vectors,
				    &ok);
      if (!ok) {
	// Choose random basis
	eigen_vectors.to_identity();
      }
      eigen_vectors.make_right_handed();

#if 0
      std::cerr << "EIGEN_BASIS: " << std::endl;
      for (std::size_t i=0; i<dimension; ++i) {
	std::cerr << "axis(" << i << ") = " << as_dual(eigen_vectors.axis(i)) << " val = " << eigen_values(i) << std::endl;
      } 
#endif

      // Build map
      typename rigid_body_map_t::matrix_t m;
      { 
	for (size_t i=0; i< dimension; ++i) {
	  for (size_t j=0; j< dimension; ++j) {
	    m(i,j) = eigen_vectors(i,j);
	  }
	  m(i,dimension) = mean(i);
	}
	m(dimension,dimension) = scalar_math<value_t>::one();
      }

      SL_TRACE_OUT(1) << "m = " << m << std::endl;

      SL_CHECK("Right handed", m.is_right_handed());

      return rigid_body_map_t(m);
    }

    /// Construct using given basis
    virtual void construct_from_basis(const rigid_body_map_t& from_box_space_map) {
      if (points_.empty()) {
	(this->last_bounding_volume_) = bounding_volume_t(from_box_space_map, vector_t());
      } else {
	static const value_t one_half = reciprocal(scalar_math<value_t>::two());

	// Construct a local AAB
	point_t p_local_inf = inverse_transformation(from_box_space_map,points_[0]);
	point_t p_local_sup = p_local_inf;
	for (size_t i=0; i<points_.size(); ++i) {
	  point_t p_local = inverse_transformation(from_box_space_map,points_[i]);
	  for (size_t d=0; d<dimension; ++d) {
	    const value_t p_local_d = p_local[d];
	    if (p_local_d > p_local_sup[d]) { p_local_sup[d] = p_local_d; }
	    if (p_local_d < p_local_inf[d]) { p_local_inf[d] = p_local_d; }
	  }
	}
	// Compute hsl and new center
	vector_t hsl = one_half * (p_local_sup-p_local_inf);
	point_t  p_center = transformation(from_box_space_map, p_local_inf.lerp(p_local_sup, one_half));

	// Recenter 
	typename rigid_body_map_t::matrix_t m = from_box_space_map.as_matrix();
	for (size_t i=0; i<dimension; ++i) {
	  m(i,dimension) = p_center[i];
	}	
	
	// Update volume
	(this->last_bounding_volume_) = bounding_volume_t(rigid_body_map_t(m), hsl);
      }
    }

    /// Build box from current cloud using specific algorithm
    virtual void construct();

  }; // class oriented_box_builder_base

  template <size_t N_dimension, class T> 
  void oriented_box_builder_base<N_dimension,T>::construct() {
    // By default, construct from eigenvectors
    this->construct_from_basis(eigen_basis());
  }

  /**
   *  Objects constructing oriented bounding boxes from point clouds in
   *  N_dimension. The method used is approximate (orientation given
   *  by the eigen-vectors of the covariance matrix)
   */
  template <size_t N_dimension, class T> 
  class oriented_box_builder: public oriented_box_builder_base<N_dimension,T> {
  public:
    enum { dimension = N_dimension };

    typedef T                                             value_t; 
    typedef oriented_box_builder<N_dimension, T>          this_t;
    typedef oriented_box_builder_base<N_dimension, T>     super_t;

    typedef oriented_box<N_dimension,T>                   bounding_volume_t;
    typedef typename bounding_volume_t::vector_t          vector_t;
    typedef typename bounding_volume_t::point_t           point_t;
    typedef typename bounding_volume_t::dual_vector_t     dual_vector_t;
    typedef typename bounding_volume_t::rigid_body_map_t  rigid_body_map_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  public: // Creation & Destruction

    oriented_box_builder() { }

    virtual ~oriented_box_builder() {}

  }; // oriented_box_builder<N,T>

  /**
   *  Objects constructing oriented bounding boxes from point clouds in
   *  two dimensions. The method used (rotating caliper on the
   *  point cloud convex hull) is optimal and runs in O(N*log(N)).
   */
  template <class T> 
  class oriented_box_builder<2,T>: public oriented_box_builder_base<2,T> {
  public:
    enum { dimension = 2 };

    typedef T                                       value_t; 
    typedef oriented_box_builder<dimension, T>      this_t;
    typedef oriented_box_builder_base<dimension, T> super_t;

    typedef oriented_box<dimension,T>               bounding_volume_t;
    typedef typename bounding_volume_t::vector_t       vector_t;
    typedef typename bounding_volume_t::point_t        point_t;
    typedef typename bounding_volume_t::dual_vector_t  dual_vector_t;
    typedef typename bounding_volume_t::rigid_body_map_t  rigid_body_map_t;
    typedef fixed_size_square_matrix<dimension+1,T> matrix_t;

    typedef convex_hull_builder<dimension,value_t, std::vector<point_t> > hull_builder_t;
    typedef typename hull_builder_t::hull_t hull_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  public: // Creation & Destruction

    oriented_box_builder() { }

    virtual ~oriented_box_builder() {}

  protected:

    /**
     *  Build minimum area box using rotating caliper algorithm O(N*log(N))
     *  TODO: reduce memory footprint by avoiding duplication of
     *        data between convex hull builder and box builder
     */
    virtual void construct();
  
    void get_aab_extrema(const hull_t& hull,
			 size_t& i_left,
			 size_t& i_right,
			 size_t& i_bottom,
			 size_t& i_top) const;

  }; // oriented_box_builder<2,T>


  template <class T>
  void oriented_box_builder<2,T>::get_aab_extrema(const hull_t& hull,
						  size_t& i_left,
						  size_t& i_right,
						  size_t& i_bottom,
						  size_t& i_top) const {
    const size_t Nedges = hull.size();
    value_t xmin = (this->points_)[hull[0][0]][0];
    value_t ymin = (this->points_)[hull[0][0]][1];
    value_t xmax = (this->points_)[hull[0][0]][0];
    value_t ymax = (this->points_)[hull[0][0]][1];
    i_left = 0;
    i_right = 0;
    i_bottom = 0;
    i_top = 0;
    for (size_t i=1; i<Nedges; ++i) {
      value_t x = (this->points_)[hull[i][0]][0];
      if (x <= xmin) {
	xmin = x;
	i_left = i;
      } else if (x >= xmax) {
	xmax = x;
	i_right = i;
      } 
      value_t y = (this->points_)[hull[i][0]][1];
      if (y <= ymin) {
	ymin = y;
	i_bottom = i;
      } else if (y >= ymax) {
	ymax = y;
	i_top = i;
      }
    }
    {
      // Deal with axis-aligned edges
      size_t i = i_left;
      while ((this->points_)[hull[i][0]][0] == 
	     (this->points_)[hull[(i+1) % Nedges][0]][0]) {
	i = (i+1)%Nedges; if (i == i_left) break;
      }
      i_left = i;
      
      i = i_right;
      while ((this->points_)[hull[i][0]][0] == 
	     (this->points_)[hull[(i+1) % Nedges][0]][0]) {
	i = (i+1)%Nedges; if (i == i_right) break;
      }
      i_right = i;
      
      i = i_bottom;
      while ((this->points_)[hull[i][0]][1] == 
	     (this->points_)[hull[(i+1) % Nedges][0]][1]) {
	i = (i+1)%Nedges; if (i == i_bottom) break;
      }
      i_bottom = i;
      
      i = i_top;
      while ((this->points_)[hull[i][0]][1] == 
	     (this->points_)[hull[(i+1) % Nedges][0]][1]) {
	i = (i+1)%Nedges; if (i == i_top) break;
      }
      i_top = i;
    }

    SL_TRACE_OUT(1) << "INDICES ARE" << i_left << " " << i_right << " " << i_bottom << " " << i_top << std::endl;
    
    SL_ENSURE("Edge check", (((this->points_)[hull[i_left][0]][0] == (this->points_)[hull[i_right][0]][0]) ||
			     ((this->points_)[hull[(i_left+1)%Nedges][0]][0] > (this->points_)[hull[i_left][0]][0])));
    SL_ENSURE("Edge check", (((this->points_)[hull[i_left][0]][0] == (this->points_)[hull[i_right][0]][0]) ||
			     ((this->points_)[hull[(i_right+1)%Nedges][0]][0] < (this->points_)[hull[i_right][0]][0])));
    SL_ENSURE("Edge check", (((this->points_)[hull[i_bottom][0]][1] == (this->points_)[hull[i_top][0]][1]) ||
			     ((this->points_)[hull[(i_bottom+1)%Nedges][0]][1] > (this->points_)[hull[i_bottom][0]][1])));
    SL_ENSURE("Edge check", (((this->points_)[hull[i_bottom][0]][1] == (this->points_)[hull[i_top][0]][1]) ||
			     ((this->points_)[hull[(i_top+1)%Nedges][0]][1] < (this->points_)[hull[i_top][0]][1])));
  }   

#if BUGGY_ROTATING_CALIPER 
 
  template <class T>
  void oriented_box_builder<2,T>::construct() {
    this->construct_from_basis(this->eigen_basis());
  }

#else
  template <class T>
  void oriented_box_builder<2,T>::construct() {
      // 1. Build the convex hull of the point set.
      SL_TRACE_OUT(1) << "BUILDING CONVEX HULL" << std::endl;

      hull_builder_t hb;
      hb.build((this->points_));
      const hull_t& hull = hb.last_hull();
      const size_t Nedges = hull.size();

      SL_TRACE_OUT(1) << "CONVEX HULL HAS " << Nedges << " EDGES" << std::endl;

      //////////////////////////////// HACK
      if (Nedges <3) {
	SL_TRACE_OUT(1) << "DEGENERATE HULL - QUICK FIX!" << std::endl;
	this->construct_from_basis(this->eigen_basis());
	return;	
      }
      //////////////////////////////// END HACK

      // 2. Compute unit-length edge directions of convex polygon
      SL_TRACE_OUT(1) << "BUILDING EDGE DIRECTIONS" << std::endl;
      std::vector<vector_t> edge(Nedges);
      std::vector<bool>     visited(Nedges);
      for (size_t i=0; i<Nedges; ++i) {
	edge[i] = ((this->points_)[hull[i][1]] - (this->points_)[hull[i][0]]).ok_normalized();
	visited[i] = false;
      }

      // 3. Find the smalles AAB containing the points. Keep track of the
      //    extremum indices i_left, i_right, B, i_top.

      SL_TRACE_OUT(1) << "FINDING AAB" << std::endl;

      size_t i_left, i_right, i_bottom, i_top;
      get_aab_extrema(hull,
		      i_left,i_right,i_bottom,i_top);

      // 4. Initialize OOB to the AAB found above
      value_t hsx = ((this->points_)[hull[i_right][0]][0]-(this->points_)[hull[i_left][0]][0])/two(value_t());
      value_t hsy = ((this->points_)[hull[i_top][0]][1]-(this->points_)[hull[i_bottom][0]][1])/two(value_t());
      value_t area_div_4 = hsx*hsy;
      (this->last_bounding_volume_) = bounding_volume_t(rigid_body_map_t(),
						vector_t(hsx, hsy));
      
      value_t min_area_div_4 = area_div_4;

      SL_TRACE_OUT(1) << "AREA IS " <<  min_area_div_4 << std::endl;

      // 5. Improve box using rotating calipers algorithm
      enum { F_NONE, F_LEFT, F_RIGHT, F_BOTTOM, F_TOP };
      vector_t u = vector_t::unit(0);
      vector_t v = vector_t::unit(1);

      SL_TRACE_OUT(1) << "BEGINNING ROTATING CALIPER " <<  std::endl;

      bool done = false;
      while ( !done ) {
	// determine edge that forms smallest angle with current box edges
	int flag = F_NONE;
	value_t max_dot = zero(value_t());
	
	value_t d = u.dot(edge[i_bottom]);
	if (d > max_dot) {
	  max_dot = d;
	  flag = F_BOTTOM;
	}
	d = v.dot(edge[i_right]);
	if (d > max_dot) {
	  max_dot = d;
	  flag = F_RIGHT;
	}
	d = -u.dot(edge[i_top]);
	if (d > max_dot) {
	  max_dot = d;
	  flag = F_TOP;
	}
	d = -v.dot(edge[i_left]);
	if (d > max_dot) {
	  max_dot = d;
	  flag = F_LEFT;
	}
	
	switch ( flag ) {
	case F_BOTTOM:
	  {
	    done = visited[i_bottom];
	    if (!done) {
	      u = edge[i_bottom];
	      v = -u[1], u[0];
	      visited[i_bottom] = true;
	      i_bottom = (i_bottom+1)%Nedges;
	    }
	  }
	  break;
	case F_RIGHT:
	  { 
	    done = visited[i_right];
	    if (!done) {
	      v = edge[i_right];
	      u = v[1], -v[0];
	      visited[i_right] = true;
	      i_right = (i_right+1)%Nedges;
	    }
	  }
	  break;
	case F_TOP:
	  {
	    done = visited[i_top];
	    if (!done) {
	      u = -edge[i_top];
	      v = -u[1], u[0];
	      visited[i_top] = true;
	      i_top = (i_top+1)%Nedges;
	    }
	  }
	  break;
	case F_LEFT:
	  {
	    done = visited[i_left];
	    if (!done) {
	      v = -edge[i_left];
	      u = v[1], -v[0];
	      visited[i_left] = true;
	      i_left = (i_left+1)%Nedges;
	    }
	  }
	  break;
	case F_NONE:
	  // polygon is a rectangle
	  {
	    done = true;
	  }
	  break;
	default:
	  SL_CHECK("Wrong tag!", false);
	}
	
	if (!done) {
	  hsx = u.dot((this->points_)[hull[i_right][0]] - (this->points_)[hull[i_left][0]])/two(value_t());
	  hsy = v.dot((this->points_)[hull[i_top][0]] - (this->points_)[hull[i_bottom][0]])/two(value_t());
	  area_div_4 = hsx*hsy;
	  if (area_div_4 < min_area_div_4) {
	    // Update current box
	    const vector_t tmp = 
	      ((this->points_)[hull[i_top][0]].as_vector()+(this->points_)[hull[i_bottom][0]].as_vector())/two(value_t()) -
	      (this->points_)[hull[i_left][0]].as_vector();
	    const point_t o = (this->points_)[hull[i_left][0]] + hsx * u + (v.dot(tmp)) * v;
	    
	    matrix_t m = tags::not_initialized();
	    m = 
	      u[0], v[0], o[0],
	      u[1], v[1], o[1],
	      zero(value_t()), zero(value_t()), one(value_t());

	    SL_CHECK("Right handed", m.is_right_handed());

	    (this->last_bounding_volume_) = bounding_volume_t(rigid_body_map_t(m),
						      vector_t(hsx, hsy));
	    
	    min_area_div_4 = area_div_4;
	  }
	}

	SL_TRACE_OUT(1) << "ROTATING CALIPER: AREA/4 = " << min_area_div_4 << " DONE = " << done << std::endl;
      } // while !done

      //////////////////////////////// HUGE HACK
      {
	this->construct_from_basis((this->last_bounding_volume_).from_box_space_map());
	if ((sl::abs(hsx-(this->last_bounding_volume_).half_side_lengths()[0]) > 100.0 * sl::epsilon(value_t()) + 0.01 * sl::max(sl::abs(hsx), sl::abs((this->last_bounding_volume_).half_side_lengths()[0]))) ||
	    (sl::abs(hsy-(this->last_bounding_volume_).half_side_lengths()[1]) > 100.0 * sl::epsilon(value_t()) + 0.01 * sl::max(sl::abs(hsy), sl::abs((this->last_bounding_volume_).half_side_lengths()[1])))) {
	  SL_TRACE_OUT(-1) << "WRONG ROTATING CALIPER HULL - QUICK FIX!" << std::endl;
	  SL_TRACE_OUT(-1) << "hsx = " << hsx << " vs. " << (this->last_bounding_volume_).half_side_lengths()[0] << std::endl;
	  SL_TRACE_OUT(-1) << "hsy = " << hsy << " vs. " << (this->last_bounding_volume_).half_side_lengths()[1] << std::endl;
	}
      }
      //////////////////////////////// END HACK

      SL_CHECK("Good caliper orientation", sl::is_non_negative(std::sqrt(sl::epsilon(min_area_div_4))+min_area_div_4));

      SL_TRACE_OUT(1) << "END ROTATING CALIPER: AREA/4 = " << min_area_div_4 << " DONE = " << done << std::endl;
    }
#endif  


  /**
   *  Objects constructing oriented bounding boxes from point clouds in
   *  three dimensions. If the user specifies a "preferred" normal vector,
   *  the points are projected onto the plane orthogonal to that vector,
   *  and the box is constructed from the optimal 2D bounding rectangle. 
   *  If no "preferred" normal vector is given, it is first computed using
   *  principal components analysis. 
   */
  template <class T> 
  class oriented_box_builder<3,T>: public oriented_box_builder_base<3,T> {
  public:
    enum { dimension = 3 };

    typedef T                                         value_t; 
    typedef oriented_box_builder<dimension, T>        this_t;
    typedef oriented_box_builder_base<dimension, T>   super_t;

    typedef oriented_box<dimension,T>                 bounding_volume_t;
    typedef typename bounding_volume_t::vector_t      vector_t;
    typedef typename bounding_volume_t::point_t       point_t;
    typedef typename bounding_volume_t::dual_vector_t dual_vector_t;
    typedef typename bounding_volume_t::rigid_body_map_t  rigid_body_map_t;
    typedef fixed_size_square_matrix<dimension+1,T>   matrix_t;

    typedef oriented_box<2,T>                         oriented_box2_t;
    typedef typename oriented_box2_t::vector_t        vector2_t;
    typedef typename oriented_box2_t::point_t         point2_t;
    typedef typename oriented_box2_t::dual_vector_t   dual_vector2_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

    dual_vector_t preferred_normal_;
    bool     has_normal_;

  public: // Creation & Destruction

    inline oriented_box_builder() { 
      has_normal_ = false; 
    }

    virtual ~oriented_box_builder() {}

    virtual void begin_model() {
      super_t::begin_model();
      has_normal_ = false;
      SL_ENSURE("No normal", !has_normal());
    }

    inline bool has_normal() const {
      return has_normal_;
    }
    
    /// The current preferred normal
    inline const dual_vector_t& normal() {
      SL_REQUIRE("Has normal", has_normal());
      return preferred_normal_;
    }

    /// Set the current preferred normal to n
    virtual void set_normal(const dual_vector_t& n) {
      SL_REQUIRE("Is building", this->is_building());
      SL_REQUIRE("No normal given", !has_normal());
      // SL_REQUIRE("Unit n", ...);

      preferred_normal_ = n;
      has_normal_ = true;
      SL_ENSURE("Has normal", has_normal());
    }

  protected:
      
    virtual void construct() {
      if (!has_normal()) {
#if BUGGY_ROTATING_CALIPER
	// By default, construct from eigenvectors
	this->construct_from_basis(this->eigen_basis());
	return;
#else 
	this->make_normal();
#endif
      }
      this->construct_from_points_and_normal();
    }

    /// Compute a preferred normal from the point cloud's eigenbasis
    virtual void make_normal() {
      SL_REQUIRE("No normal given", !has_normal());
      const rigid_body_map_t eb = this->eigen_basis();
      set_normal(as_dual(eb.axis(2))); // Get the axis with the smallest eigenvalue
      SL_ENSURE("Has normal", has_normal());
    }

    /** 
     *  Build the box given the point cloud and a preferred normal
     *  TODO: reduce memory footprint by avoiding duplication of
     *        data between 2D and 3D box builders
     */
    virtual void construct_from_points_and_normal();
  
  }; // oriented_box_builder<N,T>

  template <class T>
  void oriented_box_builder<3,T>::construct_from_points_and_normal() {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();

      SL_TRACE_OUT(1) << "n = " << normal() << std::endl;

      // Build a map that transforms n to the z-axis
      const rigid_body_map_t n_to_z = linear_map_factory<dimension,value_t>::rotation(as_dual(normal()),
										      vector_t::unit(2));
      const rigid_body_map_t z_to_n = ~n_to_z;

      SL_CHECK("Right handed", z_to_n.as_matrix().is_right_handed());

      // Project points along normal and build minimum area 2D rectangle
      value_t z_min = scalar_math<value_t>::finite_upper_bound();
      value_t z_max = scalar_math<value_t>::finite_lower_bound();
      oriented_box_builder<2,value_t> obuilder2;
      obuilder2.set_half_side_length_epsilon(sl::zero(value_t()));
      obuilder2.begin_model();
      for (size_t i=0; i<(this->points_).size(); ++i) {
	const point_t  p_prime = transformation(n_to_z, (this->points_)[i]);
	if (p_prime[2] < z_min) z_min = p_prime[2];
	if (p_prime[2] > z_max) z_max = p_prime[2];
	obuilder2.put_point(point2_t(p_prime[0], p_prime[1]));
      }
      obuilder2.end_model();
      const oriented_box2_t& box2 = obuilder2.last_bounding_volume();
      const value_t z_mid = (z_max + z_min)/two;
      const value_t z_hsl = z_max - z_mid;

      // Transform to 3D by setting z range and applying inverse transform
      const rigid_body_map<2,value_t>& map2 = box2.from_box_space_map();
      matrix_t m = tags::not_initialized();
      m = 
	map2(0,0), map2(0,1), zero, map2(0,2),
	map2(1,0), map2(1,1), zero, map2(1,2),
	zero,      zero,      one,  z_mid,
        zero,      zero,      zero, one;
      const vector2_t& xy_hsl = box2.half_side_lengths();

      SL_CHECK("Right handed", m.is_right_handed());

      SL_TRACE_OUT(1) << "m = " << m << std::endl;
      SL_TRACE_OUT(1) << "hsl = " << vector_t(xy_hsl[0], xy_hsl[1], z_hsl) << std::endl;

      (this->last_bounding_volume_) =  bounding_volume_t(z_to_n * rigid_body_map_t(m),
						 vector_t(xy_hsl[0], xy_hsl[1], z_hsl));


#if 0
      this->construct_from_basis(z_to_n * rigid_body_map_t(m));
      if ((sl::abs(xy_hsl[0]-(this->last_bounding_volume_).half_side_lengths()[0]) > 100.0 * sl::epsilon(value_t()) + 0.01 * sl::max(sl::abs(xy_hsl[0]), sl::abs((this->last_bounding_volume_).half_side_lengths()[0]))) ||
	  (sl::abs(xy_hsl[1]-(this->last_bounding_volume_).half_side_lengths()[1]) > 100.0 * sl::epsilon(value_t()) + 0.01 * sl::max(sl::abs(xy_hsl[1]), sl::abs((this->last_bounding_volume_).half_side_lengths()[1]))) ||
	  (sl::abs(z_hsl    -(this->last_bounding_volume_).half_side_lengths()[2]) > 100.0 * sl::epsilon(value_t()) + 0.01 * sl::max(sl::abs(z_hsl),     sl::abs((this->last_bounding_volume_).half_side_lengths()[2])))) {
	SL_TRACE_OUT(-1) << "WRONG PROJECTED BOX HSL - QUICK FIX!" << std::endl;
	SL_TRACE_OUT(-1) << "hsx = " << xy_hsl[0] << " vs. " << (this->last_bounding_volume_).half_side_lengths()[0] << std::endl;
	SL_TRACE_OUT(-1) << "hsy = " << xy_hsl[1] << " vs. " << (this->last_bounding_volume_).half_side_lengths()[1] << std::endl;
	SL_TRACE_OUT(-1) << "hsz = " << z_hsl     << " vs. " << (this->last_bounding_volume_).half_side_lengths()[2] << std::endl;
      }
#endif
  }

} // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  /// A 2D oriented bounding box with single precision floating point components
  typedef oriented_box<2,float> obox2f;
  /// A 3D oriented bounding box with single precision floating point components
  typedef oriented_box<3,float> obox3f;
  /// A 4D oriented bounding box with single precision floating point components
  typedef oriented_box<4,float> obox4f;

  /// A 2D oriented bounding box with double precision floating point components
  typedef oriented_box<2,double> obox2d;
  /// A 3D oriented bounding box with double precision floating point components
  typedef oriented_box<3,double> obox3d;
  /// A 4D oriented bounding box with double precision floating point components
  typedef oriented_box<4,double> obox4d;
}

#endif




