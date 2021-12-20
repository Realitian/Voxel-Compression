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
#ifndef SL_FIXED_SIZE_LINE_HPP
#define SL_FIXED_SIZE_LINE_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/linear_map.hpp>

namespace sl {

  /**
   *  Lines of fixed dimension: these are defined as
   *  R(t) = O + t*D, where t is real and D is not necessarily a unit
   *  vector
   */
  template <size_t DIMENSION, class T>
  class fixed_size_line {
  public: // Constants and types

    enum { dimension = DIMENSION, hdimension = DIMENSION+1 };

    typedef fixed_size_line<DIMENSION,T>  this_t;
    typedef T       value_t;

    typedef fixed_size_point<dimension, T>                        point_t;
    typedef fixed_size_vector<column_orientation, dimension, T>   vector_t;
    typedef fixed_size_vector<row_orientation, dimension, T>      dual_vector_t;

    typedef fixed_size_point<dimension+1, T>                      hpoint_t;
    typedef fixed_size_vector<column_orientation, dimension+1, T> hvector_t;
    typedef fixed_size_vector<row_orientation, dimension+1, T>    hdual_vector_t;

    typedef hdual_vector_t storage_t;

  protected: // Storage

    point_t  origin_;
    vector_t direction_;


  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << origin_ << direction_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> origin_ >> direction_;
    }

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", dimension > 0);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);
   
  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_line() {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_line(tags::not_initialized tag): origin_(tag), direction_(tag) {
      // Garbage in, use with care
    }

    /// Init from origin and direction
    inline explicit fixed_size_line(const point_t& o,
				   const vector_t& d)
      : origin_(o), direction_(d) {
    }

    /// Init from origin and extremity
    inline explicit fixed_size_line(const point_t& o,
				   const point_t& e)
      : origin_(o), direction_(e-o) {
    }

  public: // Access

    /// The minimum parameter value
    static inline value_t tmin() {
      return scalar_math<value_t>::finite_lower_bound();
    }
      
    /// The maximum parameter value
    static inline value_t tmax() {
      return scalar_math<value_t>::finite_upper_bound();
    }

    /// The origin of the line
    inline point_t& origin() {
      return origin_;
    }

    /// The direction of the line
    inline vector_t& direction() {
      return direction_;
    }

    /// The origin of the line
    inline const point_t& origin() const {
      return origin_;
    }

    /// The direction of the line
    inline const vector_t& direction() const {
      return direction_;
    }

  public: // Distance

    /// The square of the minimum eclidean distance to p
    inline value_t squared_distance_to(const point_t& p) const {
      vector_t d = p - origin();
      value_t l2 = direction().two_norm_squared();
      value_t t = d.dot(direction())/l2;
      d -= t*direction();
      return d.two_norm_squared();
    }

    /// The minimum eclidean distance to p
    inline value_t distance_to(const point_t& p) const {
      return std::sqrt(squared_distance_to(p));
    }
    
  public: // Transformations

    /// This transformed by map
    template <class T_tag>
    inline this_t transformed_by(const sl::linear_map_base<T_tag, dimension, value_t>& lm) const {
      return this_t(transformation(lm,origin()),
		    transformation(lm,direction()));
    }

  public: // Interpolation

    /// Linear interpolation
    template <class T_PARAMETER>
    this_t lerp(const this_t& other, T_PARAMETER t) const {
      return this_t(origin().lerp(other.origin(),t),
		    direction().lerp(other.direction(),t));
    }

  }; // class fixed_size_line

}; // namespace sl

// I/O

template <size_t DIMENSION, class T>
std::ostream& operator <<(std::ostream& s, const sl::fixed_size_line<DIMENSION,T>& a) {
  s << a.origin() << std::endl;
  s << a.direction() << std::endl;
  return s;
}
    
template <size_t DIMENSION, class T>
std::istream& operator >>(std::istream& s, sl::fixed_size_line<DIMENSION,T>& a) {
  s >> a.origin();
  s >> a.direction();
  return s;
}

/// the line x transformed by map lm
template <class T_tag, size_t N_dimension, class T>
inline sl::fixed_size_line<N_dimension,T> transformation(const sl::linear_map_base<T_tag, N_dimension, T>& lm, 
							const sl::fixed_size_line<N_dimension,T>& x) {
  return sl::fixed_size_line<N_dimension,T>(transformation(lm,x.origin()),
					   transformation(lm,x.direction()));
}

/// the line x transformed by the inverse of map lm
template <class T_tag, size_t N_dimension, class T>
inline sl::fixed_size_line<N_dimension,T> inverse_transformation(const sl::linear_map_base<T_tag, N_dimension, T>& lm, 
								const sl::fixed_size_line<N_dimension,T>& x) {
  return sl::fixed_size_line<N_dimension,T>(inverse_transformation(lm,x.origin()),
					   inverse_transformation(lm,x.direction()));
}

/// the line x transformed by map lm
template <class T_tag, size_t N_dimension, class T>
inline sl::fixed_size_line<N_dimension,T> operator * (const sl::linear_map_base<T_tag, N_dimension, T>& lm, 
						     const sl::fixed_size_line<N_dimension,T>& x) {
  return transformation(lm,x);
}

namespace sl {

  /// 2D line with single precision floating line components
  typedef fixed_size_line<2,float> line2f;
  /// 3D line with single precision floating line components
  typedef fixed_size_line<3,float> line3f;
  /// 4D line with single precision floating line components
  typedef fixed_size_line<4,float> line4f;

  /// 2D line with double precision floating line components
  typedef fixed_size_line<2,double> line2d;
  /// 3D line with double precision floating line components
  typedef fixed_size_line<3,double> line3d;
  /// 4D line with double precision floating line components
  typedef fixed_size_line<4,double> line4d;
}

#endif
