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
#ifndef SL_FIXED_SIZE_QUADRIC_MATRIX_HPP
#define SL_FIXED_SIZE_QUADRIC_MATRIX_HPP

#include <sl/fixed_size_plane.hpp>
#include <sl/fixed_size_packed_matrix.hpp>

namespace sl {
  
  /// Base class for (hyper)quadric matrices of fixed dimension
  template <class SELF_T, size_t DIMENSION, class T>
  class fixed_size_quadric_matrix_base
  {
  public: // Constants and types

    enum { dimension = DIMENSION, hdimension = DIMENSION+1 };

    typedef SELF_T  self_t;
    typedef T       value_t;

    typedef fixed_size_point<dimension, T>                        point_t;
    typedef fixed_size_vector<column_orientation, dimension, T>   vector_t;
    typedef fixed_size_vector<row_orientation, dimension, T>      dual_vector_t;
    typedef fixed_size_plane<dimension, T>                        plane_t;

    typedef fixed_size_point<dimension+1, T>                      hpoint_t;
    typedef fixed_size_vector<column_orientation, dimension+1, T> hvector_t;
    typedef fixed_size_vector<row_orientation, dimension+1, T>    hdual_vector_t;

    typedef fixed_size_packed_symmetric_matrix<dimension,T>       matrix_t;
    typedef fixed_size_square_matrix<dimension,T>                 hmatrix_t;

  protected: // Storage

    matrix_t A_;
    vector_t b_;
    value_t  c_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << A_ << b_ << c_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> A_ >> b_ >> c_;
    }

  protected: // Helpers

    static void symmetric_subfrom(matrix_t& A, const vector_t& a, const vector_t& b) {
      for(int i=0; i<A.dimension; ++i) {
        for(int j=i; j<A.dimension; ++j) {
          A(i,j) -= a[i]*b[j];
        }
      }
    }

    static void symmetric_from(matrix_t& A, const vector_t& a, const vector_t& b) {
      for(int i=0; i<A.dimension; ++i) {
        for(int j=i; j<A.dimension; ++j) {
          A(i,j) = a[i]*b[j];
        }
      }
    }

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", dimension > 0);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);
   
  protected: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_quadric_matrix_base() : c_(scalar_math<value_t>::zero()) {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_quadric_matrix_base(tags::not_initialized tag): A_(tag), b_(tag) {
      // Garbage in, use with care
    }

    inline fixed_size_quadric_matrix_base(const point_t& p1,const point_t& p2,const point_t& p3)
        : A_(sl::tags::not_initialized()), b_(sl::tags::not_initialized()) {
      vector_t e1=(p2-p1).ok_normalized();
      vector_t e2=((p3-p1)-(e1.dot(p3-p1))*e1).ok_normalized();
      
      A_.to_identity();
      symmetric_subfrom(A_, e1,e1);
      symmetric_subfrom(A_, e2,e2);

      value_t p1e1 = p1.pvdot(e1);
      value_t p1e2 = p1.pvdot(e2);
      
      b_ = e1*p1e1 + e2*p1e2 - p1.as_vector();
      c_ = p1.as_vector().two_norm_squared() - p1e1*p1e1 - p1e2*p1e2;
    }
    
    inline fixed_size_quadric_matrix_base(const vector_t& p1,const vector_t& p2,const vector_t& p3) 
        : A_(sl::tags::not_initialized()), b_(sl::tags::not_initialized()) {
      vector_t e1=(p2-p1).ok_normalized();
      vector_t e2=((p3-p1)-(e1.dot(p3-p1))*e1).ok_normalized();
      
      A_.to_identity();
      symmetric_subfrom(A_, e1,e1);
      symmetric_subfrom(A_, e2,e2);

      value_t p1e1 = p1.dot(e1);
      value_t p1e2 = p1.dot(e2);
      
      b_ = e1*p1e1 + e2*p1e2 - p1;
      c_ = p1.two_norm_squared() - p1e1*p1e1 - p1e2*p1e2;
    }

    inline fixed_size_quadric_matrix_base(const plane_t& h) 
        : A_(sl::tags::not_initialized()), b_(sl::tags::not_initialized()) {
      plane_t hn = h; hn.normalize();
      const vector_t& n = as_dual(hn.normal());
      value_t         d = hn.offset();
      
      symmetric_from(A_,n,n);
      b_ = d*n;
      c_ = d*d;
    }

    inline fixed_size_quadric_matrix_base(const matrix_t& A,
                                          const vector_t& b,
                                          const value_t& c)
        : A_(A), b_(b), c_(c) {
    }
    
    inline fixed_size_quadric_matrix_base(const self_t& other)
        : A_(other.A_), b_(other.b_), c_(other.c_) {
    }
    
  public:
    
    const matrix_t& tensor() const {
      return A_;
    }

    const vector_t& vector() const {
      return b_;
    }

    const value_t& offset() const {
      return c_;
    }

    matrix_t& tensor()  {
      return A_;
    }

    vector_t& vector()  {
      return b_;
    }

    value_t& offset()  {
      return c_;
    }
    
    hmatrix_t homogeneous() const {
      hmatrix_t H(tags::not_initialized());
      for(std::size_t i=0; i<A_.dimension; i++) {
        for(std::size_t j=0; j<A_.dimension; j++) {
          H(i,j) = A_(i,j);
        }
      }
      for(std::size_t i=0; i<b_.dimension; i++) {
        H(i, b_.dimension) = H(b_.dimension, i) = b_[i];
      }
      H(b_.dimension,b_.dimension) = c_;
      return H;
    }

    void fill(const value_t& val) {
      A_.fill(val);
      b_.fill(val);
      c_=val;
    }
    
    self_t& operator=(const self_t& Q) {
      A_=Q.A_; b_=Q.b_; c_=Q.c_; 
      return static_cast<self_t&>(*this);
    }
    
    self_t& operator+=(const self_t& Q) {
      A_+=Q.A_; b_+=Q.b_; c_+=Q.c_; 
      return static_cast<self_t&>(*this);
    }
    
    self_t& operator-=(const self_t& Q) {
      A_-=Q.A_; b_-=Q.b_; c_-=Q.c_; 
      return static_cast<self_t&>(*this);
    }

    /// change sign to all components
    inline self_t operator- () const {
      A_ = -A_; b_ = -b_; c_ = -c_; 
      return static_cast<self_t&>(*this);
    }
    
    self_t& operator*=(const value_t& s) {
      A_*=s; b_*=s; c_*=s;
      return static_cast<self_t&>(*this);
    }

    self_t& operator/=(const value_t& s) {
      A_/=s; b_/=s; c_/=s;
      return static_cast<self_t&>(*this);
    }
    
    SL_OP_LINEAR_SPACE(self_t,value_t); // Extend to generic types!

  public: // Evaluation
    
    value_t evaluate(const vector_t& v) const {
      const value_t v_A_v = v.dot(A_*v);
      const value_t b_v   = b_.dot(v);
      return v_A_v + b_v + b_v + c_;
    }

    value_t evaluate(const point_t& v) const {
      return evaluate(v.as_vector());
    }
    
    value_t operator()(const vector_t& v) const {
      return evaluate(v);
    }

  public: //Optimization
    
    void optimize_in(vector_t& v, bool *ok) const {
      matrix_t Ainv = sl::tags::not_initialized();
      A_.invert_to(Ainv,ok);
      if (ok) {
        v = -(Ainv*b_);
      }
    }

    void optimize_in(point_t& v, bool *ok) const {
      optimize_in(v.as_vector(),ok);
    }

    point_t optimized() const {
      point_t result = tags::not_initialized();
      bool ok;
      optimize_in(result, &ok);
      SL_CHECK("Works", ok);
      return result;
    }
    
    /// Optimize along line
    void optimize_in(vector_t& v, bool *ok, const vector_t& v1, const vector_t& v2) const {
      vector_t d = v1 - v2;
      const matrix_t& A = tensor();
      
      vector_t Av2 = A*v2;
      vector_t Ad  = A*d;
      
      value_t denom = sl::scalar_math<value_t>::two()*d.dot(Ad);
      *ok = !sl::is_zero(denom);
      if (ok) {
        value_t a =  sl::median( (-sl::scalar_math<value_t>::two()*vector().dot(d) - d.dot(Av2) - v2.dot(Ad) ) / denom,
                                 sl::scalar_math<value_t>::zero(),
                                 sl::scalar_math<value_t>::one());
        v = a*d + v2;
      }
    }

    /// Optimize along line
    void optimize_in(point_t& v, bool *ok, const point_t& v1, const point_t& v2) const {
      optimize_in(v.as_vector(),ok, v1.as_vector(), v2.as_vector());
    }

    point_t optimized(const point_t& v1, const point_t& v2) const {
      point_t result = tags::not_initialized();
      bool ok;
      optimize_in(result, &ok);
      if (!ok) result = v1.lerp(v2,value_t(1.0f/2.0f));
      return result;
    }

    /// Optimize inside triangle
    void optimize_in(vector_t& v, bool *ok, const vector_t& v1, const vector_t& v2, const vector_t& v3) const {
      vector_t d13 = v1 - v3;
      vector_t d23 = v2 - v3;
      const matrix_t& A = tensor();
      const vector_t& B = vector();
      vector_t Ad13 = A*d13;
      vector_t Ad23 = A*d23;
      vector_t Av3  = A*v3;

      value_t d13_d23 = d13.dot(Ad23) + d23.dot(Ad13);
      value_t v3_d13  = d13.dot(Av3)  + v3.dot(Ad13);
      value_t v3_d23  = d23.dot(Av3)  + v3.dot(Ad23);
      
      value_t d23Ad23 = d23.dot(Ad23);
      value_t d13Ad13 = d13.dot(Ad13);

      value_t denom = d13Ad13*d23Ad23 - scalar_math<value_t>::two()*d13_d23;
      ok = !sl::is_zero(denom);
      if (ok) {
        value_t a = median(( d23Ad23*(2*B.dot(d13) + v3_d13) -
                             d13_d23*(2*B.dot(d23) + v3_d23) ) / -denom,
                           scalar_math<value_t>::zero(),
                           scalar_math<value_t>::one());
        
        value_t b =  median(( d13Ad13*(2*B.dot(d23) + v3_d23) -
                              d13_d23*(2*B.dot(d13) + v3_d13) ) / -denom,
                            scalar_math<value_t>::zero(),
                            scalar_math<value_t>::one());
        v = a*d13 + b*d23 + v3;
      }    
    };

    /// Optimize inside triangle
    void optimize_in(point_t& v, bool *ok, const point_t& v1, const point_t& v2, const point_t& v3) const {
      optimize_in(v.as_vector(),ok, v1.as_vector(), v2.as_vector(), v3.as_vector());
    }

    point_t optimized(const point_t& v1, const point_t& v2, const point_t& v3) const {
      point_t result = tags::not_initialized();
      bool ok;
      optimize_in(result, &ok);
      if (!ok) result = v1.lerp(v2,value_t(1.0f/2.0f)).lerp(v3,value_t(1.0f/2.0f));
      return result;
    }

    /// fill this with zero
    inline void to_zero() {
      fill(value_t(0));
    }

    /// The "zero" origin
    static self_t zero() {
      self_t result = tags::not_initialized();
      result.to_zero();
      return result;
    }

  public: // Comparison

    /// -1 if this < t2, +1 if this > t2, 0 otherwise (sequential element comparison)
    inline int compare(const self_t& t2) const {
      int result = A_.compare(t2.A_);
      if (result == 0) {
        result = b_.compare(t2.b_);
        if (result == 0) {
          if (c_ < t2.c_) {
            result = -1;
          } else if (c_ > t2.c_) {
            result = 1;
          } else {
            result = 0;
          }
        }
      }
      return result;
    }

    /// is this < t2 (sequential element comparison
    inline bool operator<(const self_t& t2) const {
      return compare(t2) < 0;
    }

    /// is this equal to t2?
    inline bool operator == (const self_t& t2) const {
      return compare(t2) == 0;
    }

    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

    /// are all values of this within eps of the corresponding coordinate of t2?
    bool is_epsilon_equal(const self_t& t2,
                          value_t eps) const {
      return A_.is_epsilon_equal(t2.A_) &&
        b_.is_epsilon_equal(t2.b_) &&
        sl::abs(c_-t2.c_) <= eps;
    }

  }; // class fixed_size_quadric_matrix_base

}; // namespace sl

// I/O

template <class SELF_T, size_t DIMENSION, class T>
std::ostream& operator <<(std::ostream& s, const sl::fixed_size_quadric_matrix_base<SELF_T,DIMENSION,T>& a) {
  // FIXME
  return s;
}
    
template <class SELF_T, size_t DIMENSION, class T>
std::istream& operator >>(std::istream& s, sl::fixed_size_quadric_matrix_base<SELF_T,DIMENSION,T>& a) {
  //FIXME
  return s;
}


namespace sl {

  /// Quadric matrices of fixed dimension
  template <size_t DIMENSION, class T>
  class fixed_size_quadric_matrix: 
    public fixed_size_quadric_matrix_base< fixed_size_quadric_matrix<DIMENSION, T>, DIMENSION, T > {

  public: // Types

    typedef fixed_size_quadric_matrix_base< fixed_size_quadric_matrix<DIMENSION, T>, DIMENSION, T > super_t;

    enum { dimension = super_t::dimension, hdimension = super_t::hdimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::point_t       point_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;
    typedef typename super_t::plane_t       plane_t;
    typedef typename super_t::hpoint_t       hpoint_t;
    typedef typename super_t::hvector_t      hvector_t;
    typedef typename super_t::hdual_vector_t hdual_vector_t;

    typedef typename super_t::matrix_t      matrix_t;
    typedef typename super_t::hmatrix_t     hmatrix_t;

  public: // Creation, Copy & Destruction
    
    /// Default init (zero)
    inline fixed_size_quadric_matrix() {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_quadric_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    inline fixed_size_quadric_matrix(const point_t& p1,const point_t& p2,const point_t& p3) 
        : super_t(p1,p2,p3) {
    }
    
    inline fixed_size_quadric_matrix(const vector_t& p1,const vector_t& p2,const vector_t& p3)
        : super_t(p1,p2,p3) {
    }

    inline fixed_size_quadric_matrix(const plane_t& h)
        : super_t(h) {
    }

    inline fixed_size_quadric_matrix(const matrix_t& A,
                                     const vector_t& b,
                                     const value_t& c)
        : super_t(A,b,c) {
    }

    inline fixed_size_quadric_matrix(const super_t& other)
        : super_t(other) {
    }
    
  };

}; // namespace sl

namespace sl {

  /// Quadric matrices of dimension 3
  template <class T>
  class fixed_size_quadric_matrix<3,T>: 
    public fixed_size_quadric_matrix_base< fixed_size_quadric_matrix<3, T>, 3, T > {

  public: // Types

    typedef fixed_size_quadric_matrix_base< fixed_size_quadric_matrix<3, T>, 3, T > super_t;

    enum { dimension = super_t::dimension, hdimension = super_t::hdimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::point_t       point_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;
    typedef typename super_t::plane_t       plane_t;
    typedef typename super_t::hpoint_t       hpoint_t;
    typedef typename super_t::hvector_t      hvector_t;
    typedef typename super_t::hdual_vector_t hdual_vector_t;

    typedef typename super_t::matrix_t      matrix_t;
    typedef typename super_t::hmatrix_t     hmatrix_t;

  public: // Creation, Copy & Destruction
    
    /// Default init (zero)
    inline fixed_size_quadric_matrix() {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_quadric_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    inline fixed_size_quadric_matrix(const point_t& p1, const point_t& p2, const point_t& p3)
        : super_t(p1,p2,p3) {
    }
    
    inline fixed_size_quadric_matrix(const vector_t& p1, const vector_t& p2, const vector_t& p3)
        : super_t(p1,p2,p3) {
    }

    inline fixed_size_quadric_matrix(const plane_t& h)
        : super_t(h) {
    }

    inline fixed_size_quadric_matrix(const matrix_t& A,
                                     const vector_t& b,
                                     const value_t& c)
        : super_t(A,b,c) {
    }

    inline fixed_size_quadric_matrix(const super_t& other)
        : super_t(other) {
    }

  };

}; // namespace sl

namespace sl {

  /// 2D quadric_matrix with single precision floating quadric_matrix components
  typedef fixed_size_quadric_matrix<2,float> quadric_matrix2f;
  /// 3D quadric_matrix with single precision floating quadric_matrix components
  typedef fixed_size_quadric_matrix<3,float> quadric_matrix3f;
  /// 4D quadric_matrix with single precision floating quadric_matrix components
  typedef fixed_size_quadric_matrix<4,float> quadric_matrix4f;

  /// 2D quadric_matrix with double precision floating quadric_matrix components
  typedef fixed_size_quadric_matrix<2,double> quadric_matrix2d;
  /// 3D quadric_matrix with double precision floating quadric_matrix components
  typedef fixed_size_quadric_matrix<3,double> quadric_matrix3d;
  /// 4D quadric_matrix with double precision floating quadric_matrix components
  typedef fixed_size_quadric_matrix<4,double> quadric_matrix4d;
}

#endif
