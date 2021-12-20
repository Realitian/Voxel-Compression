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
#ifndef SL_FIXED_SIZE_VECTOR_HPP
#define SL_FIXED_SIZE_VECTOR_HPP

#include <sl/conv_to.hpp>
#include <sl/fixed_size_matrix.hpp>

namespace sl {

  /// The orientation of a vector
  enum vector_orientation {
    row_orientation    = 0,
    column_orientation = 1
  };

#define SL_TRANSPOSED_ORIENTATION(o) ((o) == sl::row_orientation ? sl::column_orientation : sl::row_orientation)

  /// The base class for all fixed-size vectors
  template <class SELF_T, class TRANSPOSED_T, enum vector_orientation ORIENTATION, size_t DIMENSION, class T>
  class fixed_size_vector_base: 
    public fixed_size_matrix_base<SELF_T,
    TRANSPOSED_T,
    (ORIENTATION==column_orientation?DIMENSION:1),
    (ORIENTATION==row_orientation?DIMENSION:1),
    T> 
  {
  public: // Constants and types

    enum { dimension = DIMENSION };

    typedef fixed_size_matrix_base<SELF_T, 
      TRANSPOSED_T,
      (ORIENTATION==column_orientation?DIMENSION:1),
      (ORIENTATION==row_orientation?DIMENSION:1),
      T> super_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;
    
  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", DIMENSION > 0);
    
  protected: // Creation, Copy & Destruction

    /// Default initializer (zero)
    inline fixed_size_vector_base() {
    }

    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_vector_base(tags::not_initialized tag): super_t(tag) {
    }

    /// Set this to other
    inline fixed_size_vector_base(const super_t& other): super_t(other) {
    }

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this,  dimension, v);
    }

  public: // Queries on shape
    
    /// Is this a row vector?
    static inline bool is_row_vector() {
      return super_t::row_count() == 1;
    }

    /// Is this a column vector?
    static inline bool is_column_vector() {
      return super_t::column_count() == 1;
    }
      
  public:  // Indexing

    /// Is i a good element index?
    static inline bool good_index(const size_t i) {
      return i< super_t::size();
    }

    /// The i-th component of this
    inline value_t& operator[](size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      return( this->storage_[i] );
    }

    /// The i-th component of this
    inline value_t operator[](size_t i) const {
      SL_REQUIRE("Good index", good_index(i));
      return( this->storage_[i] );
    }

    /// The i-th component of this
    inline value_t& operator()(size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      return( this->storage_[i] );
    }

    /// The i-th component of this
    inline value_t operator()(size_t i) const {
      SL_REQUIRE("Good index", good_index(i));
      return( this->storage_[i] );
    }

  public: // Initializations

    /// Set this to the unit vector aligned with coordinate axis number i
    inline void to_unit(size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      this->to_zero(); (this->storage_)[i] = scalar_math<value_t>::one();
    }

    /// The unit vector aligned with coordinate axis number i
    static inline self_t unit(size_t i) {
      self_t result = tags::not_initialized();
      result.to_unit(i);
      return result;
    }

  public: // Vector operations

    /// The dot product of this and v2
    inline SL_SUMTYPENAME(value_t) dot(const self_t& v2) const {
      return fastest::inner_product<dimension>::apply(this->begin(), v2.begin(), scalar_math<SL_SUMTYPENAME(value_t)>::zero());
    }

    /// The dot product of this and the dual vector v2
    inline SL_SUMTYPENAME(value_t) ddot(const self_t& v2) const {
      return fastest::inner_product<dimension>::apply(this->begin(), v2.begin(), scalar_math<SL_SUMTYPENAME(value_t)>::zero());
    }

    typedef SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) angle_t;

    /// The angle between this and vect
    angle_t angle(const self_t& vect) const {
      angle_t result = scalar_math<angle_t>::zero();
      const angle_t d2 = angle_t(this->two_norm_squared() * vect.two_norm_squared());
      if (d2 > scalar_math<angle_t>::epsilon()) {
        const angle_t one = scalar_math<angle_t>::one();
        const angle_t cos_angle = median(-one, one,
                                         static_cast<angle_t>(this->dot(vect)) / std::sqrt(d2));
        result = std::acos(cos_angle);
      } 
      return result;
    }

    /**
     *  Is v collinear with this within a tolerance eps? For numerical
     *  accuracy, this is implemented as sqr(cos(angle)) in [1-eps .. 1+eps]
     */
    bool is_collinear(const self_t& v2,
		      const value_t& eps = scalar_math<T>::epsilon()) const {
      bool result = true;

      const angle_t d2 = angle_t(this->two_norm_squared() * v2.two_norm_squared());
      if (d2 < scalar_math<value_t>::epsilon()) {
        // Points are coincident, assume collinear
      } else {
        // Angle test
        const angle_t one = scalar_math<angle_t>::one();
        const angle_t sqr_cos_angle = angle_t(sl::sqr(this->dot(v2))) / d2;

        result = (one - angle_t(eps)) <= sqr_cos_angle && sqr_cos_angle <= (one+angle_t(eps));
      }
      
      return result;
    }

    /// Scale this by scalar s
    inline void scale_by(const value_t& s) {
      (*this) *= s;
    }

    /// This scaled by scalar s
    inline self_t scaled_by(const value_t& s) const {
      self_t result(static_cast<const self_t&>(*this));
      result.scale_by(s);
      return result;
    }

  public: // Searching

    /// The index of the largest component
    inline size_t imax() const {
      return (this->storage_).imax();
    }

    /// The index of the component with largest absolute value
    inline size_t iamax() const {
      return (this->storage_).iamax();
    }
    
    /// The index of the smallest component
    inline size_t imin() const {
      return (this->storage_).imin();
    }

    /// The index of the component with smallest absolute value
    inline size_t iamin() const {
      return (this->storage_).iamin();
    }

  public: // Norm

    /// The sum of the absolute values of the elements
    inline SL_SUMTYPENAME(value_t) one_norm() const {
      SL_SUMTYPENAME(value_t) result;
      if (is_column_vector()) {
        result = super_t::infinite_norm();
      } else {
        result = super_t::one_norm();
      }
      return result;
    }

    /// the maximum absolute value of any of the elements
    inline SL_SUMTYPENAME(value_t) infinite_norm() const {
      SL_SUMTYPENAME(value_t) result;
      if (is_column_vector()) {
        result = super_t::one_norm();
      } else {
        result = super_t::infinite_norm();
      }
      return result;
    }
    
    /// Set vo the to the unit vector aligned with this, and ok to false if impossible 
    void normalize_to(self_t& vo, bool *ok) const {
      SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) sum2 = this->two_norm_squared();
      *ok = !is_zero(sum2);
      if (*ok) {
        vo = (*this) * sl::reciprocal_sqrt(sum2);
      }
    }

    /// The unit vector aligned with this
    self_t ok_normalized() const {
      self_t vr = tags::not_initialized();
      bool ok;
      normalize_to(vr, &ok);
      if (!ok) {
        size_t largest = iamax();
        vr = unit(largest) * non_zero_sign((*this)[largest]);
      }
      return vr;
    }

    /// The reflected direction vector. n and this must be normalized
    self_t reflected_around(const self_t& n) const {
      SL_SUMTYPENAME(value_t) twice_this_dot_n = this->dot(n);
      twice_this_dot_n += twice_this_dot_n;
      
      self_t vr = tags::not_initialized();
      for (size_t i=0; i<dimension; ++i) {
        vr[i] = value_t(twice_this_dot_n * n[i] - (*this)[i]);
      }
      return vr;
    }
    
  }; // class fixed_size_vector_base

}; // namespace sl


namespace sl {

  /// A vector of fixed dimension and orientation
  template <enum vector_orientation ORIENTATION, size_t DIMENSION, class T>
  class fixed_size_vector: 
    public fixed_size_vector_base< fixed_size_vector<ORIENTATION, DIMENSION, T>, 
    fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), DIMENSION, T>,
    ORIENTATION, DIMENSION, T > {

    public: // Types

    typedef fixed_size_vector_base< fixed_size_vector<ORIENTATION, DIMENSION, T>, 
      fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), DIMENSION, T>,
      ORIENTATION, DIMENSION, T > super_t;
    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    public: // Creation, Copy & Destruction

    /// Default intializer (zero)
    inline fixed_size_vector() {
      // Storage is already 0-filled
    }

    /// Fast initializer (garbage), handle with care!
    inline fixed_size_vector(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_vector(const super_t& other): super_t(other) {
    }

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

  };

  template <enum vector_orientation OUT_ORIENTATION, std::size_t DIM, typename OUT_ET>
  class conv_to< fixed_size_vector<OUT_ORIENTATION, DIM, OUT_ET> > {
  public:
    typedef fixed_size_vector<OUT_ORIENTATION, DIM, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_array<DIM, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<DIM; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }
    
    // Explicit conversion from vectors of another type
    template <enum vector_orientation IN_ORIENTATION, typename IN_ET> 
    inline static result_t from(const fixed_size_vector<IN_ORIENTATION, DIM, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<DIM; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }    
  }; // class conv_to


}; // namespace sl

// ---------------------------------------------------------------------------
// Matrix/Vector operations overloads
// ---------------------------------------------------------------------------

// TODO: HANDLE TYPE CONVERSIONS!!!

/// the result of vector * matrix
template <class SELF_T, class TRANSPOSED_T, size_t N_ROW, size_t N_COL, class T>
sl::fixed_size_vector<sl::row_orientation, N_COL, T> operator *(const sl::fixed_size_vector<sl::row_orientation, N_ROW, T>& vec,
                                                                const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,N_ROW, N_COL, T>& mat) {
  sl::fixed_size_vector<sl::row_orientation, N_COL, T> result = sl::tags::not_initialized();

  for (size_t i=0; i<N_COL; ++i) {
    result[i] = vec[0] * mat(0,i);
    for (size_t j=1; j<N_ROW; ++j) {
      result[i] += vec[j] * mat(j,i);
    }
  }
  return result;  
}

/// the result of matrix * vector
template <class SELF_T, class TRANSPOSED_T, size_t N_ROW, size_t N_COL, class T>
sl::fixed_size_vector<sl::column_orientation, N_ROW, T> operator *(const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,N_ROW, N_COL, T>& mat,
                                                                   const sl::fixed_size_vector<sl::column_orientation, N_COL, T>& vec) {
  sl::fixed_size_vector<sl::column_orientation, N_ROW, T> result = sl::tags::not_initialized();

  for (size_t i=0; i<N_ROW; ++i) {
    result[i] = mat(i,0) * vec[0];
    for (size_t j=1; j<N_COL; ++j) {
      result[i] += mat(i,j) * vec[j];
    }
  }  

  return result;  
}

// ---------------------------------------------------------------------------
// Matrix/Vector operations manual optimization for low dimensions
// ---------------------------------------------------------------------------

/// the result of vector * matrix
template <class SELF_T, class TRANSPOSED_T, class T>
inline sl::fixed_size_vector<sl::row_orientation, 2, T> operator *(const sl::fixed_size_vector<sl::row_orientation, 2, T>& vec,
                                                                   const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,2,2, T>& mat) {
  return sl::fixed_size_vector<sl::row_orientation, 2, T>(vec[0]*mat(0,0)+vec[1]*mat(1,0),
                                                          vec[0]*mat(0,1)+vec[1]*mat(1,1));
}

/// the result of matrix * vector
template <class SELF_T, class TRANSPOSED_T, class T>
inline sl::fixed_size_vector<sl::column_orientation, 2, T> operator *(const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,2,2, T>& mat,
                                                                      const sl::fixed_size_vector<sl::column_orientation, 2, T>& vec) {
  return sl::fixed_size_vector<sl::column_orientation, 2, T>(mat(0,0)*vec[0]+mat(0,1)*vec[1],
                                                             mat(1,0)*vec[0]+mat(1,1)*vec[1]);
}

/// the result of vector * matrix
template <class SELF_T, class TRANSPOSED_T, class T>
inline sl::fixed_size_vector<sl::row_orientation, 3, T> operator *(const sl::fixed_size_vector<sl::row_orientation, 3, T>& vec,
                                                                   const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,3,3, T>& mat) {
  return sl::fixed_size_vector<sl::row_orientation, 3, T>(vec[0]*mat(0,0)+vec[1]*mat(1,0)+vec[2]*mat(2,0),
                                                          vec[0]*mat(0,1)+vec[1]*mat(1,1)+vec[2]*mat(2,1),
                                                          vec[0]*mat(0,2)+vec[1]*mat(1,2)+vec[2]*mat(2,2));
}

/// the result of matrix * vector
template <class SELF_T, class TRANSPOSED_T, class T>
inline sl::fixed_size_vector<sl::column_orientation, 3, T> operator *(const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,3,3, T>& mat,
                                                                      const sl::fixed_size_vector<sl::column_orientation, 3, T>& vec) {
  return sl::fixed_size_vector<sl::column_orientation, 3, T>(mat(0,0)*vec[0]+mat(0,1)*vec[1]+mat(0,2)*vec[2],
                                                             mat(1,0)*vec[0]+mat(1,1)*vec[1]+mat(1,2)*vec[2],
                                                             mat(2,0)*vec[0]+mat(2,1)*vec[1]+mat(2,2)*vec[2]);
}

/// the result of vector * matrix
template <class SELF_T, class TRANSPOSED_T, class T>
inline sl::fixed_size_vector<sl::row_orientation, 4, T> operator *(const sl::fixed_size_vector<sl::row_orientation, 4, T>& vec,
                                                                   const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,4,4, T>& mat) {
  return sl::fixed_size_vector<sl::row_orientation, 4, T>(vec[0]*mat(0,0)+vec[1]*mat(1,0)+vec[2]*mat(2,0)+vec[3]*mat(3,0),
                                                          vec[0]*mat(0,1)+vec[1]*mat(1,1)+vec[2]*mat(2,1)+vec[3]*mat(3,1),
                                                          vec[0]*mat(0,2)+vec[1]*mat(1,2)+vec[2]*mat(2,2)+vec[3]*mat(3,2),
                                                          vec[0]*mat(0,3)+vec[1]*mat(1,3)+vec[2]*mat(2,3)+vec[3]*mat(3,3));
}

/// the result of matrix * vector
template <class SELF_T, class TRANSPOSED_T, class T>
inline sl::fixed_size_vector<sl::column_orientation, 4, T> operator *(const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,4,4, T>& mat,
                                                                      const sl::fixed_size_vector<sl::column_orientation, 4, T>& vec) {
  return sl::fixed_size_vector<sl::column_orientation, 4, T>(mat(0,0)*vec[0]+mat(0,1)*vec[1]+mat(0,2)*vec[2]+mat(0,3)*vec[3],
                                                             mat(1,0)*vec[0]+mat(1,1)*vec[1]+mat(1,2)*vec[2]+mat(1,3)*vec[3],
                                                             mat(2,0)*vec[0]+mat(2,1)*vec[1]+mat(2,2)*vec[2]+mat(2,3)*vec[3],
                                                             mat(3,0)*vec[0]+mat(3,1)*vec[1]+mat(3,2)*vec[2]+mat(3,3)*vec[3]);
}

// ---------------------------------------------------------------------------
// Useful conversions
// ---------------------------------------------------------------------------

/// Conversion to dual vector representation.
template <size_t DIMENSION, class T>
inline static const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& as_dual(const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& v) {
  /// Tricky implementation, using cast!
  return reinterpret_cast< const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& > (v);
}

/// Conversion to dual vector representation.
template <size_t DIMENSION, class T>
inline static const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& as_dual(const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& v) {
  /// Tricky implementation, using cast!
  return reinterpret_cast< const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& > (v);
}


/// Conversion to homogeneous coordinates, adding a zero in the last coordinate
template <enum sl::vector_orientation ORIENTATION, size_t DIMENSION, class T>
sl::fixed_size_vector<ORIENTATION,DIMENSION+1,T> as_homogeneous(const sl::fixed_size_vector<ORIENTATION,DIMENSION,T>& p) {
  sl::fixed_size_vector<ORIENTATION,DIMENSION+1,T> result = sl::tags::not_initialized();
  for (size_t i = 0; i < DIMENSION; i++) {
    result[i] = p[i];
  }
  result[DIMENSION] = sl::scalar_math<T>::zero();
  return result;
}

/// Conversion from homogeneous coordinates, expects a zero in the last coordinate
template <enum sl::vector_orientation ORIENTATION, size_t DIMENSION, class T>
sl::fixed_size_vector<ORIENTATION,DIMENSION-1,T> from_homogeneous(const sl::fixed_size_vector<ORIENTATION,DIMENSION,T>& p) {
  sl::fixed_size_vector<ORIENTATION,DIMENSION-1,T> result = sl::tags::not_initialized();
  SL_CHECK("Valid homogeneous coordinate", sl::is_zero(p[DIMENSION-1]));
  for (size_t i = 0; i < DIMENSION-1; i++) {
    result[i] = p[i];
  }
  return result;
}

// ---------------------------------------------------------------------------
// Specializations for common vector/matrix dimensions
// ---------------------------------------------------------------------------

// -- fixed_size_vector<ORIENTATION,4, T>

namespace sl {

  /// A 4D vector
  template <enum vector_orientation ORIENTATION, class T>
  class fixed_size_vector <ORIENTATION,4,T>: 
    public fixed_size_vector_base< fixed_size_vector<ORIENTATION, 4, T>, 
    fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), 4, T>,
    ORIENTATION, 4, T > {

    public: // Types

    typedef fixed_size_vector_base< fixed_size_vector<ORIENTATION, 4, T>, 
      fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), 4, T>,
      ORIENTATION, 4, T > super_t;
    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_vector() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_vector(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_vector(const super_t& other): super_t(other) {
    }

    /// Explicit init from components
    inline explicit fixed_size_vector(value_t c0,
                                      value_t c1,
                                      value_t c2,
                                      value_t c3): super_t(tags::not_initialized()) {
      (this->storage_)[0] = c0;
      (this->storage_)[1] = c1;
      (this->storage_)[2] = c2;
      (this->storage_)[3] = c3;
    }


    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

  };

}; // namespace sl

namespace sl {

  /// 3D vectors
  template <enum vector_orientation ORIENTATION, class T>
  class fixed_size_vector <ORIENTATION,3,T>: 
    public fixed_size_vector_base< fixed_size_vector<ORIENTATION, 3, T>, 
    fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), 3, T>,
    ORIENTATION, 3, T > {

    public: // Types

    typedef fixed_size_vector_base< fixed_size_vector<ORIENTATION, 3, T>, 
      fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), 3, T>,
      ORIENTATION, 3, T > super_t;
    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    public: // Creation, Copy & Destruction

    /// Default initializer (zero)
    inline fixed_size_vector() {
      // Storage is already 0-filled
    }

    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_vector(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_vector(const super_t& other): super_t(other) {
    }

    /// Explicit init from components
    inline explicit fixed_size_vector(value_t c0,
                                      value_t c1,
                                      value_t c2): super_t(tags::not_initialized()) {
      (this->storage_)[0] = c0;
      (this->storage_)[1] = c1;
      (this->storage_)[2] = c2;
    }

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

    public: // Special features

    inline self_t cross(const self_t& v2) const { 
      value_t v1x=(*this)[0], v1y=(*this)[1], v1z=(*this)[2];
      value_t v2x=v2[0], v2y=v2[1], v2z=v2[2];
      return self_t(v1y*v2z-v1z*v2y, v1z*v2x-v1x*v2z, v1x*v2y-v1y*v2x);
    }       

    inline transposed_t normal(const self_t& v2) const { 
      value_t v1x=(*this)[0], v1y=(*this)[1], v1z=(*this)[2];
      value_t v2x=v2[0], v2y=v2[1], v2z=v2[2];
      return transposed_t(v1y*v2z-v1z*v2y, v1z*v2x-v1x*v2z, v1x*v2y-v1y*v2x);
    }       

    // Get the signed angle in radian between two vectors
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) signed_angle(const self_t& vect) const {
      typedef SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) angle_t;
      return atan2(angle_t(this->cross(vect).two_norm()), 
                   angle_t(this->dot(vect)));
      
    }
             
  };

}; // namespace sl

namespace sl {

  /// 2D Vectors
  template <enum vector_orientation ORIENTATION, class T>
  class fixed_size_vector <ORIENTATION,2,T>: 
    public fixed_size_vector_base< fixed_size_vector<ORIENTATION, 2, T>, 
    fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), 2, T>,
    ORIENTATION, 2, T > {

    public: // Types

    typedef fixed_size_vector_base< fixed_size_vector<ORIENTATION, 2, T>, 
      fixed_size_vector<SL_TRANSPOSED_ORIENTATION(ORIENTATION), 2, T>,
      ORIENTATION, 2, T > super_t;
    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_vector() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_vector(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_vector(const super_t& other): super_t(other) {
    }

    /// Explicit init from components
    inline explicit fixed_size_vector(value_t c0,
                                      value_t c1): super_t(tags::not_initialized()) {
      (this->storage_)[0] = c0;
      (this->storage_)[1] = c1;
    }

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

  };

  template <enum sl::vector_orientation ORIENTATION, std::size_t dimension, class T>
  sl::fixed_size_vector<ORIENTATION,dimension,T> centroid(std::size_t N, 
							  const sl::fixed_size_vector<ORIENTATION,dimension,T>* v) {
    sl::fixed_size_vector<ORIENTATION,dimension,T> result;
    for (std::size_t i=0; i<N; ++i) {
      result += v[i];
    }
    T scale = sl::reciprocal(T(N));
    for (std::size_t d=0; d<dimension; ++d) {
      result[d] *= scale;
    }
    return result;
  }
  
  template <enum sl::vector_orientation ORIENTATION, std::size_t N, class T>
  inline std::size_t hash_value(const sl::fixed_size_vector<ORIENTATION,N,T>& a) {
    std::size_t seed = 0;
    for (std::size_t i=0; i<N; ++i) {
      hash_combine(seed, a[i]);
    }
    return seed;
  }

} // namespace sl

// ---------------------------------------------------------------------------
// Matrix/Vector operations overloads
// ---------------------------------------------------------------------------

template <enum sl::vector_orientation ORIENTATION, size_t DIMENSION, class T>
inline sl::fixed_size_vector<ORIENTATION,DIMENSION,T> operator*(const T& y, const sl::fixed_size_vector<ORIENTATION,DIMENSION,T>& x) {
  typedef sl::fixed_size_vector<ORIENTATION,DIMENSION,T> self_t;
  self_t result = sl::tags::not_initialized();
   
  const typename self_t::restrict_iterator result_it = result.restrict_begin();
  const typename self_t::const_iterator    src_it    = x.begin();
  for (size_t i=0; i<self_t::element_size; ++i) {
    result_it[i] = src_it[i] * y;
  }
  return result;
}


// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {

  /// 2D column vector with single precision floating point components
  typedef fixed_size_vector<sl::column_orientation,2,float> column_vector2f;
  /// 3D column vector with single precision floating point components
  typedef fixed_size_vector<sl::column_orientation,3,float> column_vector3f;
  /// 4D column vector with single precision floating point components
  typedef fixed_size_vector<sl::column_orientation,4,float> column_vector4f;
  
  /// 2D row vector with single precision floating point components
  typedef fixed_size_vector<sl::row_orientation,2,float> row_vector2f;
  /// 3D row vector with single precision floating point components
  typedef fixed_size_vector<sl::row_orientation,3,float> row_vector3f;
  /// 4D row vector with single precision floating point components
  typedef fixed_size_vector<sl::row_orientation,4,float> row_vector4f;

  /// 2D column vector with single precision floating point components
  typedef column_vector2f vector2f;
  /// 3D column vector with single precision floating point components
  typedef column_vector3f vector3f;
  /// 4D column vector with single precision floating point components
  typedef column_vector4f vector4f;

  /// 2D column vector with double precision floating point components
  typedef fixed_size_vector<sl::column_orientation,2,double> column_vector2d;
  /// 3D column vector with double precision floating point components
  typedef fixed_size_vector<sl::column_orientation,3,double> column_vector3d;
  /// 4D column vector with double precision floating point components
  typedef fixed_size_vector<sl::column_orientation,4,double> column_vector4d;
  
  /// 2D row vector with double precision floating point components
  typedef fixed_size_vector<sl::row_orientation,2,double> row_vector2d;
  /// 3D row vector with double precision floating point components
  typedef fixed_size_vector<sl::row_orientation,3,double> row_vector3d;
  /// 4D row vector with double precision floating point components
  typedef fixed_size_vector<sl::row_orientation,4,double> row_vector4d;

  /// 2D column vector with double precision floating point components
  typedef column_vector2d vector2d;
  /// 3D column vector with double precision floating point components
  typedef column_vector3d vector3d;
  /// 4D column vector with double precision floating point components
  typedef column_vector4d vector4d;

  /// HACK, porting from sl 3.0 - Will be replaced by a color class
  typedef vector3f color3f;
  /// HACK, porting from sl 3.0 - Will be replaced by a color class
  typedef vector4f color4f;

  /// HACK, porting from sl 3.0 - Will be replaced by a color class
  typedef vector3d color3d;
  /// HACK, porting from sl 3.0 - Will be replaced by a color class
  typedef vector4d color4d;

} // namespace sl

#endif 



