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
#ifndef SL_FIXED_SIZE_POINT_HPP
#define SL_FIXED_SIZE_POINT_HPP

#include <sl/conv_to.hpp>
#include <sl/fixed_size_vector.hpp>

namespace sl {


  // The possible orderings of three points
  typedef enum point_ordering { ORDER_CURVED_CCW, 
				ORDER_CURVED_CW, 
				ORDER_COLLINEAR_BEFORE,
				ORDER_COLLINEAR_INBETWEEN,
				ORDER_COLLINEAR_AFTER
  } point_ordering_t;

  /// Base class for points of fixed dimension
  template <class SELF_T, size_t DIMENSION, class T>
  class fixed_size_point_base
  {
  public: // Constants and types

    enum { dimension = DIMENSION };

    typedef SELF_T  self_t;
    typedef T       value_t;

    typedef fixed_size_vector<column_orientation, DIMENSION, T> vector_t;
    typedef fixed_size_vector<row_orientation, DIMENSION, T>    dual_vector_t;

    typedef vector_t storage_t;

  protected: // Storage

    vector_t storage_;

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", DIMENSION > 0);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);
   
  protected: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_point_base() {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_point_base(tags::not_initialized tag): storage_(tag) {
      // Garbage in, use with care
    }

  public: // STL-style iterators

    typedef typename storage_t::iterator          iterator;
    typedef typename storage_t::restrict_iterator restrict_iterator;
    typedef typename storage_t::const_iterator    const_iterator;

    /// The number of components
    static inline size_t size()                { return dimension; }
    
    /// Iterator pointing to first component
    inline iterator begin()                    { return storage_.begin(); }
    /// Restrict iterator pointing to first component
    inline iterator restrict_begin()  { return storage_.restrict_begin(); }
    /// Const iterator pointing to first component 
    inline const_iterator begin() const        { return storage_.begin(); }
    /// Iterator pointing to end of storage
    inline iterator end()                      { return storage_.end(); }
    /// Restrict iterator pointing to end of storage
    inline iterator restrict_end()    { return storage_.restrict_end(); }
    /// Const iterator pointing to end of storage
    inline const_iterator end() const          { return storage_.end(); }

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << storage_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> storage_;
    }

  public:  // Indexing
    
    /// Is i a good component index
    inline bool good_index(const size_t i) const {
      return i< size();
    }

    /// The i-th component
    inline value_t& operator[](size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// The i-th component
    inline value_t operator[](size_t i) const {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// The i-th component
    inline value_t& operator()(size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// The i-th component
    inline value_t operator()(size_t i) const {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// Pointer to storage area
    inline value_t* to_pointer() {
      return( storage_.to_pointer() );
    }
    
    /// Pointer to storage area
    inline const value_t* to_pointer() const {
      return( storage_.to_pointer() );
    }

    /// The storage area
    inline const vector_t& storage() const {
      return storage_.storage();
    } 

    /// this converted to vector
    inline vector_t& as_vector() {
      return storage_;
    }

    /// this converted to vector
    inline const vector_t& as_vector() const {
      return storage_;
    }

  public: // Initialization

    /** Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

    /// fill this with c0
    inline void fill(value_t c0) {
      storage_.fill(c0);
    }

    /// fill this with zero
    inline void to_zero() {
      fill(value_t(0));
    }

    /// set this to the point located at unit distance from the origin on axis i
    inline void to_unit(size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      storage_.to_unit(i);
    }

    /// The origin
    static self_t zero() {
      self_t result = tags::not_initialized();
      result.to_zero();
      return result;
    }

    /// the point located  at unit distance from the origin on axis i
    static self_t unit(size_t i) {
      SL_USEVAR(i);
      self_t result = tags::not_initialized();
      result.to_unit();
      return result;
    }

    /// a point with all coordinates set to s
    static self_t constant(const value_t s) {
      self_t result = tags::not_initialized();
      result.to_constant(s);
      return result;
    }

  public: // Comparison

    /// -1 if this < t2, +1 if this > t2, 0 otherwise (sequential element comparison)
    inline int compare(const self_t& t2) const {
      return storage_.compare(t2.storage_);
    }

    /// is this < t2 (sequential element comparison
    inline bool operator<(const self_t& t2) const {
      return storage_ < t2.storage_;
    }

    /// is this equal to t2?
    inline bool operator == (const self_t& t2) const {
      return storage_ == t2.storage_;
    }

    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

    /// are all coordinates of this with eps of the corresponding coordinate of t2?
    inline bool is_epsilon_equal(const self_t& t2,
			  value_t eps) const {
      return storage_.is_epsilon_equal(t2.storage_, eps);
    }
    
  public: // Algebra
     
    /// Add the vector other to this
    template <class SELF_T2, class TRANSPOSED_T2, class VALUE_T2>
    inline self_t& operator+= (const fixed_size_vector_base<SELF_T2,TRANSPOSED_T2,column_orientation,dimension,VALUE_T2>& other) {
      storage_ += other;
      return static_cast<self_t&>(*this);
    }
	
    /// Subtract vector other from this
    template <class SELF_T2, class TRANSPOSED_T2, class VALUE_T2>
    inline self_t& operator-= (const fixed_size_vector_base<SELF_T2,TRANSPOSED_T2,column_orientation,dimension,VALUE_T2>& other) {
      storage_ -= other;
      return static_cast<self_t&>(*this);
    }
	
    /// The vector resulting from subtracting point other from this
    template <class SELF_T2, class VALUE_T2>
    inline fixed_size_vector<column_orientation,dimension,SL_PROMOTENAME(value_t,VALUE_T2)>
    operator- (const fixed_size_point_base<SELF_T2,dimension,VALUE_T2>& other) const {
       fixed_size_vector<column_orientation,dimension,SL_PROMOTENAME(value_t,VALUE_T2)> result = as_vector();
       return result -= other.as_vector();
    }
   
    // TODO: Handle type conversions
    SL_OP_ADDABLE2(self_t, vector_t);
    SL_OP_SUBTRACTABLE2(self_t, vector_t);

    /// Scale this by scalar s
    inline void scale_by(const value_t& s) {
      storage_ *= s;
    }

    /// This scaled by scalar s
    inline self_t scaled_by(const value_t& s) const {
      self_t result(static_cast<const self_t&>(*this));
      result.scale_by(s);
      return result;
    }

  public: // Interpolation

    /// Linear interpolation
    template <class T_PARAMETER>
    self_t lerp(const self_t& other, T_PARAMETER t) const {
      self_t result = tags::not_initialized();
      T_PARAMETER one_minus_t = sl::one(t) - t;
      for (size_t i=0; i<dimension; ++i) {
	result[i] = one_minus_t * (*this)[i] + t * other[i];
      }
      return result;
    }

  public: // Point operations

    /// point vector inner product (i.e. this converted to vector dot v2)
    inline SL_SUMTYPENAME(value_t) pvdot(const vector_t& v2) const {
      return as_vector().dot(v2);
    }                     

    /// point vector inner product (i.e. this converted to vector dot v2)
    inline SL_SUMTYPENAME(value_t) pvdot(const dual_vector_t& v2) const {
      return pvdot(as_dual(v2));
    }                     

    /// Squared Euclidean distance from this to p2
    inline SL_SUMTYPENAME(value_t) distance_squared_to(const self_t& p2) const {
      return (p2-(*this)).two_norm_squared();
    }
 
    /// Euclidean distance from this to p2
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) distance_to(const self_t& p2) const {
      return (p2-(*this)).two_norm();
    }

    /**
     *  Squared Euclidean distance from this to triangle p0, p2, p2. Optionally returns
     *  also the barycentric coordinates of closest point.
     */
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) distance_squared_to(const self_t& p0,
                                                                         const self_t& p1,
                                                                         const self_t& p2,
                                                                         value_t*      p0p1_parameter = 0,
                                                                         value_t*      p0p2_parameter = 0) const {
      typedef  SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) result_t;
      SL_REQUIRE("Implementation limitation", dimension == 3);

      const result_t zero = sl::scalar_math<result_t>::zero();
      const result_t one  = sl::scalar_math<result_t>::one();
      const result_t two  = sl::scalar_math<result_t>::two();
      
      vector_t p_p0  = p0-*this;
      vector_t p0_p1 = p1-p0;
      vector_t p0_p2 = p2-p0;
      
      result_t a00 = p0_p1.two_norm_squared();
      result_t a01 = p0_p1.dot(p0_p2);
      result_t a11 = p0_p2.two_norm_squared();
      result_t b0  = p_p0.dot(p0_p1);
      result_t b1  = p_p0.dot(p0_p2);
      result_t c   = p_p0.two_norm_squared();
      result_t det = sl::abs(a00*a11-a01*a01);
      result_t alpha = a01*b1-a11*b0;
      result_t beta  = a01*b0-a00*b1;
      result_t d2;

      SL_TRACE_OUT(2) << "P: " << (*this)[0] << " " << (*this)[1] << " " << (*this)[2] << std::endl;
      SL_TRACE_OUT(2) << "P0: " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
      SL_TRACE_OUT(2) << "P1: " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
      SL_TRACE_OUT(2) << "P2: " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
      SL_TRACE_OUT(2) << "a : " << a00 << " " << a01 << " " << a11 << std::endl;
      SL_TRACE_OUT(2) << "b : " << b0 << " " << b1 << std::endl;
      SL_TRACE_OUT(2) << "c : " << c << std::endl;
      
      if ( alpha + beta <= det ) {
        if (sl::is_negative(alpha)) {
          if (sl::is_negative(beta)) { // region 4
            if (sl::is_negative(b0)) {
              beta = zero;
              if (-b0 >= a00) {
                alpha = one;
                d2 = a00+two*b0+c; SL_TRACE_OUT(2) << "C0: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
              } else {
                alpha = -b0/a00;
                d2 = b0*alpha+c; SL_TRACE_OUT(2) << "C1: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
              }
            } else {
              alpha = zero;
              if (sl::is_non_negative(b1)) {
                beta = zero;
                d2 = c; SL_TRACE_OUT(2) << "C2: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
              } else if (-b1 >= a11) {
                beta = one;
                d2 = a11+two*b1+c; SL_TRACE_OUT(2) << "C3: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
              } else {
                beta = -b1/a11;
                d2 = b1*beta+c; SL_TRACE_OUT(2) << "C4: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
              }
            }
          } else {  // region 3
            alpha = zero;
            if (sl::is_non_negative(b1)) {
              beta = zero;
              d2 = c; SL_TRACE_OUT(2) << "C5: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
            } else if (-b1 >= a11) {
              beta = one;
              d2 = a11+two*b1+c; SL_TRACE_OUT(2) << "C6: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
            } else {
              beta = -b1/a11;
              d2 = b1*beta+c; SL_TRACE_OUT(2) << "C7: (" << alpha << " " << beta << ") -> " << d2 << std::endl;   
            }
          }
        } else if (sl::is_negative(beta)) { // region 5
          beta = zero;
          if (sl::is_non_negative(b0)) {
            alpha = zero;
            d2 = c; SL_TRACE_OUT(2) << "C8: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
          } else if ( -b0 >= a00 ) {
            alpha = one;
            d2 = a00+two*b0+c; SL_TRACE_OUT(2) << "C9: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
          } else {
            alpha = -b0/a00;
            d2 = b0*alpha+c; SL_TRACE_OUT(2) << "C10: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
          }
        } else {
          // region 0: minimum at interior point
          if (sl::is_zero(det)) {
            // Degenerate: check segments
            // FIXME: Is this stable enough? Use epsilon above?
            result_t d2_01;
            result_t d2_02;
#if 0
            if (!sl::is_positive(-b0)) {
              d2_01 = c; // Dist to point 0
            } else if (-b0 >= a00) {
              d2_01 = (p_p0 + p0_p1).two_norm_squared(); // Dist to point 0
            } else {
              d2_01 = (p_p0 - b0 * p0_p1).two_norm_squared(); // Dist to seg 01
            }
            if (!sl::is_positive(-b1)) {
              d2_02 = c; // Dist to point 0
            } else if (-b1 >= a11) {
              d2_02 = (p_p0 + p0_p2).two_norm_squared(); // Dist to point 0
            } else {
              d2_02 = (p_p0 - b1 * p0_p2).two_norm_squared(); // Dist to seg 02
            }
            d2 = std::min(d2_01,d2_02);            
            alpha = zero;
            beta = zero;
#else
            if (sl::is_non_negative(b0)) {
              alpha = zero;
              d2_01 = c; 
            } else if ( -b0 >= a00 ) {
              alpha = one;
              d2_01 = a00+two*b0+c; 
            } else {
              alpha = -b0/a00;
              d2_01 = b0*alpha+c; 
            }
            if (sl::is_non_negative(b1)) {
              beta = zero;
              d2_02 = c; 
            } else if (-b1 >= a11) {
              beta = one;
              d2_02 = a11+two*b1+c; 
            } else {
              beta = -b1/a11;
              d2_02 = b1*beta+c;
            }
            if (d2_02<d2_01) {
              alpha = zero; d2 = d2_02;
            } else {
              beta = zero; d2 = d2_01;
            }
#endif
            SL_TRACE_OUT(2) << "C11: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
          } else {
            result_t inv_det = sl::reciprocal(det);
            alpha *= inv_det;
            beta *= inv_det;
            d2 =
              alpha*(a00*alpha+a01*beta+two*b0)+
              beta*(a01*alpha+a11*beta+two*b1)+
              c; SL_TRACE_OUT(2) << "C12: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
          }
        }
      } else {
        result_t tmp0, tmp1, numerator, denominator;
        if (sl::is_negative(alpha)) {  // region 2
          tmp0 = a01 + b0;
          tmp1 = a11 + b1;
          if (tmp1 > tmp0) {
            numerator = tmp1 - tmp0;
            denominator = a00-two*a01+a11;
            if ( numerator >= denominator ) {
              alpha = one;
              beta = zero;
              d2 = a00+two*b0+c; SL_TRACE_OUT(2) << "C13: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else {
              alpha = numerator/denominator;
              beta = one-alpha;
              d2 =
                alpha*(a00*alpha+a01*beta+two*b0)+
                beta*(a01*alpha+a11*beta+two*b1)+
                c;  SL_TRACE_OUT(2) << "C14: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            }
          } else {
            alpha = zero;
            if (sl::is_non_positive(tmp1)) {
              beta = one;
              d2 = a11+two*b1+c; SL_TRACE_OUT(2) << "C15: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else if (sl::is_non_negative(b1)) {
              beta = zero;
              d2 = c; SL_TRACE_OUT(2) << "C16: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else {
              beta = -b1/a11;
              d2 = b1*beta+c; SL_TRACE_OUT(2) << "C17: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            }
          }
        } else if (sl::is_negative(beta)) {  // region 6
          tmp0 = a01 + b1;
          tmp1 = a00 + b0;
          if ( tmp1 > tmp0 ) {
            numerator = tmp1 - tmp0;
            denominator = a00-two*a01+a11;
            if ( numerator >= denominator ) {
              beta = one;
              alpha = zero;
              d2 = a11+two*b1+c; SL_TRACE_OUT(2) << "C18: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else {
              beta = numerator/denominator;
              alpha = one - beta;
              d2 =
                alpha*(a00*alpha+a01*beta+two*b0)+
                beta*(a01*alpha+a11*beta+two*b1)+
                c; SL_TRACE_OUT(2) << "C19: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            }
          } else {
            beta = zero;
            if (sl::is_non_positive(tmp1)) {
              alpha = one;
              d2 = a00+two*b0+c; SL_TRACE_OUT(2) << "C20: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else if (sl::is_non_negative(b0)) {
              alpha = zero;
              d2 = c; SL_TRACE_OUT(2) << "C21: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else {
              alpha = -b0/a00;
              d2 = b0*alpha+c; SL_TRACE_OUT(2) << "C22: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            }
          }
        } else {  // region 1
          numerator = a11 + b1 - a01 - b0;
          if (sl::is_non_positive(numerator)) {
            alpha = zero;
            beta = one;
            d2 = a11+two*b1+c; SL_TRACE_OUT(2) << "C23: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
          } else {
            denominator = a00-two*a01+a11;
            if ( numerator >= denominator ) {
              alpha = one;
              beta = zero;
              d2 = a00+two*b0+c; SL_TRACE_OUT(2) << "C24: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            } else {
              alpha = numerator/denominator;
              beta = one-alpha;
              d2 =
                alpha*(a00*alpha+a01*beta+two*b0)+
                beta*(a01*alpha+a11*beta+two*b1)+
                c; SL_TRACE_OUT(2) << "C25: (" << alpha << " " << beta << ") -> " << d2 << std::endl;
            }
          }
        }
      }
      
      if (p0p1_parameter) *p0p1_parameter = alpha;
      if (p0p2_parameter) *p0p2_parameter = beta;
      return sl::abs(d2);
    }

    /**
     *  Euclidean distance from this to triangle p0, p2, p2. Optionally returns
     *  also the barycentric coordinates of closest point.
     */
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) distance_to(const self_t& p0,
                                                                 const self_t& p1,
                                                                 const self_t& p2,
                                                                 value_t*      p0p1_parameter = 0,
                                                                 value_t*      p0p2_parameter = 0) const {
      return std::sqrt(distance_squared_to(p0,p1,p2, p0p1_parameter, p0p2_parameter));
    }

    /**
     *  Projection of this point on plane defined by (p, n)
     */
    inline void project_on_plane(const self_t& p, const vector_t& n) {
      vector_t p2q  = storage_-p.as_vector();
      vector_t dpq = p2q - (p2q.dot(n)) * n;
      storage_= (p + dpq).storage_;
    }

    /// This projected on plane defined by (p, n)
    inline self_t projected_on_plane(const self_t& p, const vector_t& n) const {
      self_t result(static_cast<const self_t&>(*this));
      result.project_on_plane(p,n);
      return result;
    }

    /**
     *  Projection of this point on plane defined by (p, n)
     */
    inline void project_on_plane(const self_t& p, const dual_vector_t& n) {
      return project_on_plane(p, as_dual(n));
    }
    
    /// This projected on plane defined by (p, n)
    inline self_t projected_on_plane(const self_t& p, const dual_vector_t& n) const {
      self_t result(static_cast<const self_t&>(*this));
      result.project_on_plane(p,n);
      return result;
    }
    
  }; // class fixed_size_point_base
		     
}; // namespace sl

// I/O

template <class SELF_T, size_t DIMENSION, class T>
std::ostream& operator <<(std::ostream& s, const sl::fixed_size_point_base<SELF_T,DIMENSION,T>& a) {
  for (size_t i = 0; i<a.size(); i++) {
      s << a(i) << std::endl;;
  }
  return s;
}
    
template <class SELF_T, size_t DIMENSION, class T>
std::istream& operator >>(std::istream& s, sl::fixed_size_point_base<SELF_T,DIMENSION,T>& a) {
  for (size_t i = 0; i<a.size(); i++) {
      s >> a(i) ;
  }
  return s;
}

namespace sl {

  /** 
   *  Are P0, P1, P2 collinear within a tolerance eps?
   */
  template <class SELF_T, size_t DIMENSION, class T>
  inline  bool are_epsilon_collinear(const fixed_size_point_base<SELF_T,DIMENSION,T>& p0,
                                     const fixed_size_point_base<SELF_T,DIMENSION,T>& p1,
                                     const fixed_size_point_base<SELF_T,DIMENSION,T>& p2,
                                     const T& eps = scalar_math<T>::epsilon()) {
    SL_USEVAR(eps);
    return (p1-p0).is_collinear(p2-p0);
  }
  
  /// the square of the area of the n-dimensional triangle p1,p2,p3 
  template <class SELF_T, size_t DIMENSION, class T>
  inline T triangle_area_squared(const fixed_size_point_base<SELF_T,DIMENSION,T>& p1, 
                                 const fixed_size_point_base<SELF_T,DIMENSION,T>& p2, 
                                 const fixed_size_point_base<SELF_T,DIMENSION,T>& p3) {
    // Stable floating point Heron formula, see Kahan
    T a = (p2-p1).two_norm();
    T b = (p3-p2).two_norm();
    T c = (p1-p3).two_norm();
    if (b>a) std::swap(a,b); // a>=b
    if (c>a) std::swap(a,c); // a>= b ; a>=c 
    if (c>b) std::swap(c,b); // a>= b >= c
    return std::max(T(0.0), T(0.0625)*(a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)));
  }

  /// the area of the n-dimensional triangle p1,p2,p3 
  template <class SELF_T, size_t DIMENSION, class T>
  inline T triangle_area(const fixed_size_point_base<SELF_T,DIMENSION,T>& p1, 
                         const fixed_size_point_base<SELF_T,DIMENSION,T>& p2, 
                         const fixed_size_point_base<SELF_T,DIMENSION,T>& p3) {
    return std::sqrt(triangle_area_squared(p1,p2,p3));
  }
  
  /**
   *  the compactness of the triangle p1,p2,p3: one for an equilateral,
   *  triangle, zero for a triangle whose vertices are colinear
   */
  template <class SELF_T, size_t DIMENSION, class T>
  inline T triangle_compactness(const fixed_size_point_base<SELF_T,DIMENSION,T>& p1, 
                                const fixed_size_point_base<SELF_T,DIMENSION,T>& p2, 
                                const fixed_size_point_base<SELF_T,DIMENSION,T>& p3) {
    T FOUR_ROOT3 = T(6.928203230275509);

    T L1 = (p2-p1).two_norm_squared();
    T L2 = (p3-p2).two_norm_squared();
    T L3 = (p1-p3).two_norm_squared();
    T L = L1+L2+L3;
    if (is_positive(L)) {
      return FOUR_ROOT3 * triangle_area(p1, p2, p3) / L;
    } else {
      return scalar_math<T>::one();
    }
  }             
  
  /// cross product of p2-p2 with p3-p1
  template <class SELF_T, class T>
  inline typename fixed_size_point_base<SELF_T,3,T>::vector_t cross(const fixed_size_point_base<SELF_T,3,T>& p1, 
                                                                    const fixed_size_point_base<SELF_T,3,T>& p2, 
                                                                    const fixed_size_point_base<SELF_T,3,T>& p3) {
    return (p2-p1).cross(p3-p1);
  }          
  
  /// the normal to the triangle p1,p2,p3 (not normalized)
  template <class SELF_T, class T>
  inline typename fixed_size_point_base<SELF_T,3,T>::dual_vector_t normal(const fixed_size_point_base<SELF_T,3,T>& p1, 
                                                                          const fixed_size_point_base<SELF_T,3,T>& p2, 
                                                                          const fixed_size_point_base<SELF_T,3,T>& p3) {
    return (p2-p1).normal(p3-p1);
  }     
  
  /**
   *  The order of p wrt to the segment q0 -> q1, using a tolerance eps
   *  to determine collinearity. The possible result values are:
   *  <ul>
   *  <li>ORDER_CURVED_CCW: The segment q0q1 is oriented CCW when seed from p
   *  <li>ORDER_CURVED_CW : The segment q0q1 is oriented CW when seed from p
   *  <li>ORDER_COLLINEAR_BEFORE: The point p is on the line q0q1 positioned before q0
   *  <li>ORDER_COLLINEAR_AFTER : The point p is on the line q0q1 positioned after q1
   *  <li>ORDER_COLLINEAR_INBETWEEN : The point p is on the segment q0q1 
   *  </ul>
   *
   *  TODO: Generalize to N-dimensions
   */
   template <class SELF_T, class T>
   point_ordering_t order(const fixed_size_point_base<SELF_T,2,T>& p,
			  const fixed_size_point_base<SELF_T,2,T>& q0,
			  const fixed_size_point_base<SELF_T,2,T>& q1,
			  const T& eps = scalar_math<T>::epsilon()) {
     
     typedef T                                                    value_t;
     typedef typename fixed_size_point_base<SELF_T,2,T>::vector_t vector_t;

     point_ordering_t result = ORDER_COLLINEAR_INBETWEEN;

     const vector_t d = q1 - q0;
     const vector_t a = p - q0;
     const value_t d_dot_d = d.dot(d);
     const value_t a_dot_a = a.dot(a);
     const value_t det_da = d[0]*a[1] - d[1]*a[0];

     const value_t r = det_da*det_da - eps * d_dot_d * a_dot_a;
     
     if (is_positive(r)) {
       // Curved
       if (is_positive(det_da)) {
	 result = ORDER_CURVED_CCW;
       } else {
	 SL_CHECK("Negative determinant", is_negative(det_da));
	 result = ORDER_CURVED_CW;
       }
     } else {
       // Collinear
       const value_t d_dot_a = d.dot(a);
       if (is_negative(d_dot_a)) {
	 result = ORDER_COLLINEAR_BEFORE;
       } else if (d_dot_a > d_dot_d) {
	 result = ORDER_COLLINEAR_AFTER;
       }
     }
     
     return result;
   }
 
}

namespace sl {

  /// Points of fixed dimension
  template <size_t DIMENSION, class T>
  class fixed_size_point: 
    public fixed_size_point_base< fixed_size_point<DIMENSION, T>, DIMENSION, T > {

  public: // Types

    typedef fixed_size_point_base< fixed_size_point<DIMENSION, T>, DIMENSION, T > super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_point() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_point(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /** Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }
  };

  template <std::size_t DIM, typename OUT_ET>
  class conv_to< fixed_size_point<DIM, OUT_ET> > {
  public:
    typedef fixed_size_point<DIM, OUT_ET> result_t;

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
    
    // Explicit conversion from points of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_point<DIM, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<DIM; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }    
  }; // class conv_to
 
} // namespace sl

namespace sl {

  /// 4D points
  template <class T>
  class fixed_size_point<4,T>: 
    public fixed_size_point_base< fixed_size_point<4, T>, 4, T > {

  public: // Types

    typedef fixed_size_point_base< fixed_size_point<4, T>, 4, T > super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;

  public: // Creation, Copy & Destruction

    /// Default initializer (zero)
    inline fixed_size_point() {
      // Storage is already 0-filled
    }

    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_point(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from components
    inline explicit fixed_size_point(const value_t c0,
				     const value_t c1,
				     const value_t c2,
				     const value_t c3): super_t(tags::not_initialized()) {
      this->storage_[0] = c0;
      this->storage_[1] = c1;
      this->storage_[2] = c2;
      this->storage_[3] = c3;
    }
    

    /** Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }
  };

}; // namespace sl

namespace sl {

  /// 3D points
  template <class T>
  class fixed_size_point<3,T>: 
    public fixed_size_point_base< fixed_size_point<3, T>, 3, T > {

  public: // Types

    typedef fixed_size_point_base< fixed_size_point<3, T>, 3, T > super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;

  public: // Creation, Copy & Destruction

    /// Default initializer (zero)
    inline fixed_size_point() {
      // Storage is already 0-filled
    }

    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_point(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from components
    inline explicit fixed_size_point(value_t c0,
				     value_t c1,
				     value_t c2): super_t(tags::not_initialized()) {
      this->storage_[0] = c0;
      this->storage_[1] = c1;
      this->storage_[2] = c2;
    }
    
    /** Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

  };

}; // namespace sl

namespace sl {

  /// 2D points
  template <class T>
  class fixed_size_point<2,T>: 
    public fixed_size_point_base< fixed_size_point<2, T>, 2, T > {

  public: // Types

    typedef fixed_size_point_base< fixed_size_point<2, T>, 2, T > super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;

  public: // Creation, Copy & Destruction

    /// Default initializer (zero)
    inline fixed_size_point() {
      // Storage is already 0-filled
    }

    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_point(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from components
    inline explicit fixed_size_point(value_t c0,
				value_t c1): super_t(tags::not_initialized()) {
      this->storage_[0] = c0;
      this->storage_[1] = c1;
    }
        
    /** Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }

  public: // Specific 2D operations

  };

  /**
   *  Intersection of ray with a triangle.
   */
  template <class T>
  void triangle_intersection_in(bool *hit,
				T         *hit_t,
				T         *hit_u,
				T         *hit_v,
				const sl::fixed_size_point<3,T>& org,
				const sl::fixed_size_vector<sl::column_orientation,3,T>& dir,
				T          t_min,
				T         t_max,
				const sl::fixed_size_point<3,T>& p0,
				const sl::fixed_size_point<3,T>& p1,
				const sl::fixed_size_point<3,T>& p2) {
    typedef T value_t;
    // Moller Trumbore 97 algorithm
    // compute u ,v barycentric coordinates of the intersection.
    // ray		r(t)  = o + d*t
    // triangle point	p(u,v)= (1-u-v)*v0 + u*v1 + v*v2
    
    const value_t Epsilon = std::numeric_limits<T>::epsilon();
    
    sl::fixed_size_vector<sl::column_orientation,3,T> e1 = p1 - p0;
    sl::fixed_size_vector<sl::column_orientation,3,T> e2 = p2 - p0;
    sl::fixed_size_vector<sl::column_orientation,3,T> p = dir.cross( e2 );
    value_t        det = e1.dot( p );
    
    // check if determinant is equal to 0
    if ( det > -Epsilon && det < Epsilon ) {
      *hit = false;
      return;
    }
    
    value_t inv_det = 1.0f/det;
    
    // find and check first barycentric coordinate
    sl::fixed_size_vector<sl::column_orientation,3,T> s = org - p0;
    value_t u = inv_det * s.dot( p ) ;
    if ( u < 0.0f || u > 1.0f ) {
      *hit = false;
      return;
    }
    
    // find and check second barycentric coordinate
    sl::fixed_size_vector<sl::column_orientation,3,T> q = s.cross( e1 );
    value_t v = inv_det * dir.dot( q );
    if ( v < 0.0f || u + v > 1.0f ) {
      *hit = false;
	return;
    }
    
    // intersect: find ray parameter
    value_t t = inv_det * e2.dot( q );
    if (t < t_min || t > t_max) {
      *hit = false;
      return;
    }
    
    *hit = true;
    *hit_t = t;
    *hit_u = u;
    *hit_v = v;
  }

  template <std::size_t dimension, class T>
  sl::fixed_size_point<dimension,T> centroid(std::size_t N, 
					     const sl::fixed_size_point<dimension,T>* p) {
    sl::fixed_size_point<dimension,T> result;
    for (std::size_t i=0; i<N; ++i) {
      result += p[i].as_vector();
    }
    T scale = sl::reciprocal(T(N));
    for (std::size_t d=0; d<dimension; ++d) {
      result[d] *= scale;
    }
    return result;
  }
  
  template <size_t N, class T>
  inline std::size_t hash_value(const fixed_size_point<N,T>& a) {
    std::size_t seed = 0;
    for (std::size_t i=0; i<N; ++i) {
      hash_combine(seed, a[i]);
    }
    return seed;
  }

}; // namespace sl

/// conversion to homogeneous coordinates (last component set to 1)
template <size_t DIMENSION, class T>
sl::fixed_size_point<DIMENSION+1,T> as_homogeneous(const sl::fixed_size_point<DIMENSION,T>& p) {
  sl::fixed_size_point<DIMENSION+1,T> result = sl::tags::not_initialized();
  for (size_t i = 0; i < DIMENSION; i++) {
    result[i] = p[i];
  }
  result[DIMENSION] = sl::one(T());;
  return result;
}

/// conversion from homogeneous coordinates (division by last component)
template <size_t DIMENSION, class T>
sl::fixed_size_point<DIMENSION-1,T> from_homogeneous(const sl::fixed_size_point<DIMENSION,T>& p) {
  sl::fixed_size_point<DIMENSION-1,T> result = sl::tags::not_initialized();
  T w = p[DIMENSION-1];
  SL_CHECK("Valid homogeneous coordinate", !sl::is_zero(w));
  w = sl::reciprocal(w);
  for (size_t i = 0; i < DIMENSION-1; i++) {
    result[i] = p[i] * w;
  }
  return result;
}

/// conversion of vector v to the point origin + v (data is shared!)
template <enum sl::vector_orientation ORIENTATION, size_t DIMENSION, class T> 
inline static const sl::fixed_size_point<DIMENSION,T>& as_point(const sl::fixed_size_vector<ORIENTATION,DIMENSION,T>& v) {
  // tricky implementation
  return reinterpret_cast< const sl::fixed_size_point<DIMENSION,T>& >(v);
}

/// conversion of point p to the vector p - origin (data is shared!)
template <size_t DIMENSION, class T> 
inline static const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& as_column_vector(const sl::fixed_size_point<DIMENSION,T>& p) {
  // tricky implementation
  return reinterpret_cast< const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& > (p);
}

/// conversion of point p to the dual of the vector p - origin (data is shared!)
template <size_t DIMENSION, class T> 
inline static const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& as_row_vector(const sl::fixed_size_point<DIMENSION,T>& p) {
  return reinterpret_cast< const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& > (p);
}

/// conversion of point p to the vector p - origin (alias of as_colum_vector)  (data is shared!)
template <size_t DIMENSION, class T> 
inline static const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& as_vector(const sl::fixed_size_point<DIMENSION,T>& p) {
  return reinterpret_cast< const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& > (p);
}

/// conversion of point p to the vector p - origin (alias of as_row_vector)  (data is shared!)
template <size_t DIMENSION, class T> 
inline static const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& as_dual_vector(sl::fixed_size_point<DIMENSION,T>& p) {
  return reinterpret_cast< const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& > (p);
}

SL_OP_ADDABLE2_OVERLOADS(SL_OP_TEMPLATE2(template <size_t DIMENSION, class T>),
			 SL_OP_TEMPLATE2(sl::fixed_size_point<DIMENSION,T>),
			 SL_OP_TEMPLATE3(sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>));
SL_OP_SUBTRACTABLE2_OVERLOADS(SL_OP_TEMPLATE2(template <size_t DIMENSION, class T>),
			      SL_OP_TEMPLATE2(sl::fixed_size_point<DIMENSION,T>),
			      SL_OP_TEMPLATE3(sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>));


namespace sl {

  /// 2D point with single precision floating point components
  typedef fixed_size_point<2,float> point2f;
  /// 3D point with single precision floating point components
  typedef fixed_size_point<3,float> point3f;
  /// 4D point with single precision floating point components
  typedef fixed_size_point<4,float> point4f;

  /// 2D point with double precision floating point components
  typedef fixed_size_point<2,double> point2d;
  /// 3D point with double precision floating point components
  typedef fixed_size_point<3,double> point3d;
  /// 4D point with double precision floating point components
  typedef fixed_size_point<4,double> point4d;
}

#endif
