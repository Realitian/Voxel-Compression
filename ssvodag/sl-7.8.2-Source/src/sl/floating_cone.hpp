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
#ifndef SL_FLOATING_CONE_HPP
#define SL_FLOATING_CONE_HPP

#include <sl/ball.hpp>
#include <sl/quaternion.hpp>
#include <sl/linear_map_factory.hpp>

namespace sl {
  
  /**
   * A N-dimensional floating_cone, i.e. a cone of directions
   */
  template <size_t N_dimension, class T> 
  class floating_cone {
  public:
    enum { dimension = N_dimension };

    typedef floating_cone<dimension,T> this_t;
    typedef T                          value_t;
    
    typedef fixed_size_vector<column_orientation,N_dimension,T> vector_t;
    typedef fixed_size_vector<row_orientation,N_dimension,T>    dual_vector_t;
    typedef fixed_size_point<N_dimension,T>                     point_t;
    typedef ball<N_dimension,T>                                 ball_t;
    
  protected:

    vector_t axis_;
    value_t  sin_angle_;
    value_t  cos_angle_;
    
  public: // Serialization
    
    inline void store_to(output_serializer& s) const {
      s << axis_ << sin_angle_ << cos_angle_;
    }
    
    inline void retrieve_from(input_serializer& s) {
      s >> axis_ >> sin_angle_ >> cos_angle_;
    }

  public: // Creation & Destruction

    /// Default init (empty cone)
    inline floating_cone() { 
      axis_ = vector_t::unit(0); // Random
      sin_angle_ = sl::scalar_math<value_t>::zero();
      cos_angle_ = sl::scalar_math<value_t>::zero();
    }

    /// Init from axis only, angle assumed to be zero
    inline floating_cone(const vector_t& n): axis_(n), sin_angle_(sl::scalar_math<value_t>::zero()), cos_angle_(sl::scalar_math<value_t>::one()) {
    }

    /// Init from axis, angle. 
    inline floating_cone(const vector_t& n,
                         const value_t& alpha): axis_(n), sin_angle_(std::sin(alpha)), cos_angle_(std::cos(alpha)) {
    }

    /// The from axis, sin angle, cos angle
    inline explicit floating_cone(const vector_t& n,
                                  const value_t& sin_alpha,
                                  const value_t& cos_alpha ): axis_(n), sin_angle_(sin_alpha), cos_angle_(cos_alpha) {
    }

    /// Init from eye point and looked up ball
    inline floating_cone(const point_t& e,
                         const ball_t& b) {
      axis_ = b.center() - e;
      value_t a = axis_.two_norm();
      value_t R  = b.radius();
      if (a <= R) {
        // Full cone -- angle = 180
        axis_ = vector_t::unit(0); // Random
        cos_angle_ = -sl::scalar_math<value_t>::one();
        sin_angle_ =  sl::scalar_math<value_t>::zero();
      } else {
	axis_ /= a;
	sin_angle_ = R / a;
	cos_angle_ = std::sqrt(1.0f-sin_angle_*sin_angle_);
      }
    }
      
  public: // Access

    /// The floating_cone center
    inline const vector_t& axis() const {
      return axis_;
    }
	
    /// The floating_cone sin angle
    inline const value_t& sin_angle() const { 
      return sin_angle_;
    }

    /// The floating_cone cos angle
    inline const value_t& cos_angle() const { 
      return cos_angle_;
    }

    /// The floating cone angle
    inline value_t angle() const {
      return std::acos(sl::median(cos_angle(),
                                  -sl::scalar_math<value_t>::one(),
                                  sl::scalar_math<value_t>::one()));
    }
    
  public: // Setting

    inline void set_axis(const vector_t& x) {
      axis_ = x;
    }

    inline void set_angle(const value_t& x) {
      cos_angle_ = std::cos(x);
      sin_angle_ = std::sin(x);
    }
    
  public: // Queries
		
    /// True iff the floating_cone is empty (i.e. radius is negative)
    inline bool is_empty() const {
      return sl::is_zero(cos_angle()) && sl::is_zero(sin_angle());
    }

    /// Is angle < Pi/2?
    inline bool is_acute() const {
      return sl::is_positive(cos_angle());
    }
    
  public: // Assignment

    /// Set this to the zero-sized floating_cone containing p
    inline void to(const vector_t& n) {
      axis_ = n;
      cos_angle_ = sl::scalar_math<value_t>::one();
      sin_angle_ = sl::scalar_math<value_t>::zero();
    }

    /// Set this to the empty floating_cone
    inline void to_empty() {
      axis_.to_zero();
      cos_angle_ = sl::scalar_math<value_t>::zero();
      sin_angle_ = sl::scalar_math<value_t>::zero();
    }
	
    /// Set this to the largest representable floating_cone
    inline void to_huge() {
      axis_ = vector_t::unit(0); // Random
      cos_angle_ = -sl::scalar_math<value_t>::one();
      sin_angle_ =  sl::scalar_math<value_t>::zero();
    }

    /// find rotation axis to move  other.axis to axis
    vector_t rotation_axis(const vector_t& other ) {
      static value_t const epsilon = value_t(0.0001);
      vector_t rot_axis = axis().cross( other );
      value_t rot_norm = rot_axis.two_norm();
      if ( rot_norm > epsilon ) 
	return (rot_axis / rot_norm);

      // parallel axes, try to find an orthogonal axis with cross Ex
      vector_t e_x;
      e_x[ 0 ] = sl::scalar_math<value_t>::one();
      rot_axis = axis().cross( e_x );
      rot_norm = rot_axis.two_norm();
      if ( rot_norm > epsilon ) 
	return (rot_axis / rot_norm);

      // axis parallel to Ex, find an orthogonal axis with cross Ey
      vector_t e_y;
      e_y[ 1 ] = sl::scalar_math<value_t>::one(); 
      rot_axis = axis().cross( e_y ).ok_normalized();
      return rot_axis;
    }
 
    /// Merge floating_cone other with this, eventually adjusting axis and angle
    void merge(const this_t& other) {
      value_t delta_angle = std::acos(sl::median( axis().dot( other.axis() ),
						  -sl::scalar_math<value_t>::one(),
						  sl::scalar_math<value_t>::one())); // O..Pi
      if (other.is_empty()) {
        // nothing to do
      } else if (is_empty()) {
        (*this) = other;
      } else if ( angle() >= delta_angle + other.angle() ) {
	// this contain other, nothing to do
      } else if ( other.angle() >= delta_angle + angle() ) {
	//  other contain this, set this
	(*this) = other;
      } else {
	// merge cones
	vector_t rot_axis = rotation_axis( other.axis() );
	value_t result_angle = std::min( ( delta_angle + angle() + other.angle() ) / 2.0f,
					 sl::scalar_math<value_t>::Pi());
	value_t rot_angle = result_angle - angle();
	
	sl::quaternion<value_t> q;
	q.from_axis_angle( rot_axis, rot_angle );
	axis_ = sl::linear_map_factory<N_dimension, value_t>::rotation( q ) * axis();
	cos_angle_ = std::cos(result_angle);
        sin_angle_ = std::sin(result_angle);
      }
    }

    /// Merge vector other with this, eventually adjusting axis and angle
    void merge(const vector_t& other) {
      if (is_empty()) {
        to(other);
      } else {
	this_t other_cone( other, sl::scalar_math<value_t>::zero() );
	merge( other_cone );
      }
    }
    
    /// Merge floating_cone other with this, keeping the same axis
    void grow(const this_t& other) {
      if (other.is_empty()) {
        // nothing to do
      } else if (is_empty()) {
        (*this) = other;
      } else {
        value_t result_angle2 = std::acos(sl::median(axis().dot(other.axis()),
                                                     -sl::scalar_math<value_t>::one(), 
                                                     sl::scalar_math<value_t>::one())); // O..Pi
        value_t result_angle = std::min(result_angle2+other.angle(),
                                        sl::scalar_math<value_t>::Pi());
        if (result_angle > angle()) {
          cos_angle_ = std::cos(result_angle);
          sin_angle_ = std::sin(result_angle);
        }
      }
    }
    
    /// Merge vector other with this, keeping the same axis
    void grow(const vector_t& other) {
      if (is_empty()) {
        to(other);
      } else {
#if 0
        const value_t result_angle = std::acos(sl::median(axis().dot(other),
							  -sl::scalar_math<value_t>::one(), 
							  sl::scalar_math<value_t>::one())); // O..Pi
        if (result_angle > angle()) {
          cos_angle_ = std::cos(result_angle);
          sin_angle_ = std::sin(result_angle);
        }
#else
	const value_t dot = sl::median(axis().dot(other),
				       -sl::scalar_math<value_t>::one(), 
				       sl::scalar_math<value_t>::one());
	if (dot<cos_angle_) {
	  cos_angle_ = dot;
	  sin_angle_ = std::sqrt(sl::scalar_math<value_t>::one() - dot*dot);
	}
#endif
      }
    }
    
  public: // Contains

    bool contains(const vector_t& d) const {
      return
        (!is_empty()) &&
        (axis().angle(d) <= angle()); // FIXME speedup, remove trig
    }    
 
   bool contains(const this_t& other) const {
     if ( is_empty() ) {
       return false;
     } else if ( other.is_empty() ) {
       return true;
     }

     value_t delta_angle = std::acos(sl::median( axis().dot( other.axis() ),
						 -sl::scalar_math<value_t>::one(),
						 sl::scalar_math<value_t>::one())); // O..Pi
      if ( angle() >= delta_angle + other.angle() ) {
	return true;
      } else {
	return false;
      }
    }
    
  public: // Direction culling

    /// Is this fully backfacing with respect to view cone?
    bool is_fully_backfacing(const ball_t& bsphere, const point_t& eye) const {
      // From Qsplat
      bool result = false;
      if (is_empty() || !is_acute()) {
	// No check -- assume false
      } else {
	vector_t cam = eye - bsphere.center();
	value_t camdotnorm = cam.dot(axis_);
	if (camdotnorm < -bsphere.radius()) {
	  value_t camdist2 = cam.two_norm_squared();
	  if (sl::sqr(camdotnorm + bsphere.radius()) > camdist2 * sl::sqr(sin_angle_)) {
	    result = true;
	  }
	}
      }
      return result;
    }
    
    /// Is this fully frontfacing with respect to view cone?
    bool is_fully_frontfacing(const ball_t& bsphere, const point_t& eye) const {
      // From Qsplat
      bool result = false;
      if (is_empty() || !is_acute()) {
	// No check -- assume false
      } else {
	vector_t cam = eye - bsphere.center();
	value_t camdotnorm = cam.dot(axis_);
	if (camdotnorm > bsphere.radius()) {
	  value_t camdist2 = cam.two_norm_squared();
	  if (sl::sqr(camdotnorm - bsphere.radius()) > camdist2 * sl::sqr(sin_angle_)) {
	    result = true;
	  }
	}
      }
      return result;
    }
			     
    /// Is this fully backfacing with respect to view cone?
    bool is_fully_backfacing(const this_t& view_cone) const {
      if (is_empty()) {
        return false;
      } else if (view_cone.is_empty()) {
        return false;
      } else if (!view_cone.is_acute()) {
        return false;
      } else if (!is_acute()) {
        return false;
      } else {
        value_t vdotn = view_cone.axis().dot(axis());
        if (sl::is_positive(vdotn)) {
          return vdotn + view_cone.sin_angle() > sin_angle();
        } else {
          // Eye is in frontfacing side
          return false;
        }
      }        
    }

    /// Is this fully frontfacing with respect to view cone?
    bool is_fully_frontfacing(const this_t& view_cone) const {
      if (is_empty()) {
        return false;
      } else if (view_cone.is_empty()) {
        return false;
      } else if (!view_cone.is_acute()) {
        return false;
      } else if (!is_acute()) {
        return false;
      } else {
        value_t vdotn = view_cone.axis().dot(axis());
        if (sl::is_negative(vdotn)) {
          return -vdotn - view_cone.sin_angle() > sin_angle();
        } else {
          // Eye is in backfacing side
          return false;
        }
      }        
    }
    
  public: // Transformations

    /// the floating_cone containing this transformed by map 
    template <class T_tag>
    this_t transformed_by(const sl::linear_map_base<T_tag, dimension, value_t>& m) const {
      SL_USEVAR(m);
      this_t result;
      SL_FAIL("Not implemented");
      return result;
    }

  public: // Comparison

    /// -1 if this < t2, +1 if this > t2, 0 otherwise (sequential element comparison)
    inline int compare(const this_t& t2) const {
      int result = axis_.compare(t2.axis_);
      if (result == 0) {
        if (sin_angle_ < t2.sin_angle_) {
          result = -1;
        } else if (sin_angle_ > t2.sin_angle_) {
          result = 1;
        }
      }
      if (result == 0) {
        if (cos_angle_ < t2.cos_angle_) {
          result = -1;
        } else if (cos_angle_ > t2.cos_angle_) {
          result = 1;
        }
      }
      return result;
    }

    /// is this < t2 (sequential element comparison
    inline bool operator<(const this_t& t2) const {
      return compare(t2) < 0;
    }

    /// is this equal to t2?
    inline bool operator == (const this_t& t2) const {
      return compare(t2) == 0;
    }

    SL_OP_COMPARABLE1(this_t);
    SL_OP_EQUALITY_COMPARABLE1(this_t);
    
  }; // floating_cone

  template <std::size_t DIM, typename OUT_ET>
  class conv_to< floating_cone<DIM, OUT_ET> > {
  public:
    typedef floating_cone<DIM, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const floating_cone<DIM, IN_ET>& in) {
      return result_t(conv_to< fixed_size_vector<column_orientation,DIM,OUT_ET> >(in.axis()),
		      static_cast<OUT_ET> (in.sin_angle()),
		      static_cast<OUT_ET> (in.cos_angle()));
    }
  };

}; // namespace sl

/// the axis_aligned floating_cone b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::floating_cone<N_dimension,T> transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
                                                       const sl::floating_cone<N_dimension,T>& b) {
  return b.transformed_by(map);
}

/// the axis_aligned floating_cone b transformed by the inverse of map
template <class T_tag, size_t N_dimension, class T>
inline sl::floating_cone<N_dimension,T> inverse_transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
                                                               const sl::floating_cone<N_dimension,T>& b) {
  return b.transformed_by(map.inverse());
}

/// the floating_cone b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::floating_cone<N_dimension,T> operator * (const sl::linear_map_base<T_tag, N_dimension, T>& map, 
                                                    const sl::floating_cone<N_dimension,T>& b) {
  return b.transformed_by(map);
}

// I/O

template <size_t N_dimension, class T>
std::ostream& operator <<(std::ostream& s, const sl::floating_cone<N_dimension,T>& b) {
  s << b.axis() << std::endl;;
  s << b.angle() << std::endl;;
  return s;
}
    
template <size_t N_dimension, class T>
std::istream& operator >>(std::istream& s, sl::floating_cone<N_dimension,T>& b) {
  typename sl::floating_cone<N_dimension,T>::vector_t axis;
  typename sl::floating_cone<N_dimension,T>::value_t angle;
  s >> axis;
  s >> angle;
  b.set_axis(axis);
  b.set_angle(angle);
  return s;
}

namespace sl {

  /// An object for building bounding floating_cones from point/floating_cone clouds
  template <size_t N_dimension, class T> 
  class floating_cone_builder {
  public:
    enum { dimension = N_dimension };

    typedef floating_cone_builder<N_dimension, T> this_t;
    typedef floating_cone<dimension,T>            cone_t;
    typedef T                                     value_t;
    
    typedef fixed_size_vector<column_orientation,N_dimension,T> vector_t;
    typedef fixed_size_vector<row_orientation,N_dimension,T>    dual_vector_t;
    typedef fixed_size_point<N_dimension,T>                     point_t;
    
  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  protected: // State

    bool is_quick_and_dirty_;
    
    bool is_building_;

    cone_t last_bounding_cone_;

    std::vector<vector_t> vectors_;
    
  public: // Creation & Destruction

    /// Default init
    inline floating_cone_builder() {
      is_quick_and_dirty_ = true;
      is_building_ = false;
    }

    inline virtual ~floating_cone_builder() {
    }
  
  public: // Construction

   /// is the object currently building a bounding volume?
    inline bool is_building() const {
      return is_building_;
    }

    /// Start point/floating_cone cloud
    virtual inline void begin_model() {
      SL_REQUIRE("Not is building", !is_building());
      is_building_ = true;
      SL_ENSURE("Is building", is_building());
    }
    
    /// Add point to current cloud
    virtual inline void put_vector(const vector_t& p) {
      SL_REQUIRE("Is building", is_building());
      vectors_.push_back(p);
    }

    /// Finish point/floating_cone cloud and build floating_cone
    virtual inline void end_model() {
      SL_REQUIRE("Is building", is_building());
      last_bounding_cone_.to_empty();

      const std::size_t N = vectors_.size();
      if (N>0) {
	// Find axis
	vector_t axis;
	if (!is_quick_and_dirty_) {
	  // Optimal but slow bounding sphere approach
	  sl::ball_builder<dimension,value_t> mbp;
	  mbp.begin_model();
	  for (std::size_t i=0; i<vectors_.size(); ++i) {
	    mbp.put_point(as_point(vectors_[i]));
	  }
	  mbp.end_model();
	  axis = mbp.last_bounding_volume().center().as_vector().ok_normalized();
	} else {
	  // Bbox approximation, Shirman and Abi-Ezzi, 1993
	  vector_t lo = vectors_[0];
	  vector_t hi = lo;
	  for (std::size_t i=1; i<N; ++i) {
	    for (std::size_t d=0; d<dimension; ++d) {
	      lo[d] = std::min(lo[d], vectors_[i][d]);
	      hi[d] = std::max(hi[d], vectors_[i][d]);
	    }
	  }
	  axis = (lo+hi).ok_normalized();
	}
	// Find angle
        last_bounding_cone_ = cone_t(axis);
        for (std::size_t i=0; i<vectors_.size(); ++i) {
          last_bounding_cone_.grow(vectors_[i]);
        }      
        vectors_.clear();
      }
      
      is_building_ = false;
      SL_ENSURE("Not is building", !is_building());
    }
    
    /// the bounding volume constructed by last end_model
    inline const cone_t& last_bounding_cone() const {
      return last_bounding_cone_;
    }

  }; // class floating_cone_builder

}; // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  // A 2D floating_cone with single precision floating point components
  typedef floating_cone<2,float> floating_cone2f;
  // A 3D floating_cone with single precision floating point components
  typedef floating_cone<3,float> floating_cone3f;
  // A 4D floating_cone with single precision floating point components
  typedef floating_cone<4,float> floating_cone4f;

  // A 2D floating_cone with double precision floating point components
  typedef floating_cone<2,double> floating_cone2d;
  // A 3D floating_cone with double precision floating point components
  typedef floating_cone<3,double> floating_cone3d;
  // A 4D floating_cone with double precision floating point components
  typedef floating_cone<4,double> floating_cone4d;
}

#endif
