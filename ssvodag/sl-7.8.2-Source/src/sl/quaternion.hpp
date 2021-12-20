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
#ifndef SL_QUATERNION_HPP
#define SL_QUATERNION_HPP

#include <sl/conv_to.hpp>
#include <sl/fixed_size_vector.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/math.hpp>

namespace sl {

  // --------------------------------------------------------------------
  // -- quaternion<T>
  // --------------------------------------------------------------------

  /**
   * Quaternions, i.e. a four dimensional complex number that
   * can represent three-dimensional rotations.
   */
  template <class T>
  class quaternion : 
    public fixed_size_vector_base< quaternion<T>,
                                   quaternion<T>,
                                   column_orientation,
                                   4, SL_FLOATTYPENAME(T) > 
  {
  public: // Constants and types

    typedef fixed_size_vector_base< quaternion<T>,
                               quaternion<T>,
                               column_orientation,
                               4, SL_FLOATTYPENAME(T) > super_t;

    enum { dimension = super_t::dimension };

    typedef fixed_size_vector<column_orientation, 3, SL_FLOATTYPENAME(T)> axis_t;
    typedef fixed_size_square_matrix<4, SL_FLOATTYPENAME(T)>              matrix_t;

    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

  public:

    SL_COMPILE_TIME_CHECK("Good dimension", dimension == 4);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);

  public: // Element access

    inline value_t x() const { return (*this)[0]; }
    inline value_t y() const { return (*this)[1]; }
    inline value_t z() const { return (*this)[2]; }
    inline value_t w() const { return (*this)[3]; }

  public: // Creation, Copy & Destruction

    /// Default initialization (identity)
    inline quaternion(): super_t(tags::not_initialized()) {
      (*this)[0] = sl::zero(value_t());
      (*this)[1] = sl::zero(value_t());
      (*this)[2] = sl::zero(value_t());
      (*this)[3] =  sl::one(value_t());       
    }

    /// Null initialization, use with care!
    inline quaternion(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Copy constructor
    inline quaternion(const super_t& other): super_t(other) {
    }

    /// Explicit initialization
    inline quaternion(const value_t& x,
		      const value_t& y,
		      const value_t& z,
		      const value_t& w): super_t(tags::not_initialized()) {
      (*this)[0] = x;
      (*this)[1] = y;
      (*this)[2] = z;
      (*this)[3] = w;
    }


  public: // Checks

    /// Is this norm equal to one?
    inline bool is_normalized() const {
      return sl::abs(this->two_norm_squared() - sl::one(value_t())) < static_cast<value_t>(0.001);
    }

  public: // Conversions

    /// Set this to the quaternion representation of a rotation of angle around axis. Axis must be normalized.
    void from_axis_angle(const axis_t& axis, 
			 const value_t& angle) {
      SL_REQUIRE("Unit axis", sl::abs(axis.two_norm_squared() - sl::one(value_t())) < static_cast<value_t>(0.001));

      const value_t half_angle = angle/two(value_t());
      const value_t sinus   = std::sin(half_angle);
      const value_t cosinus = std::cos(half_angle);
      (*this)[0] = axis[0] * sinus;
      (*this)[1] = axis[1] * sinus;
      (*this)[2] = axis[2] * sinus;
      (*this)[3] = cosinus;

      SL_ENSURE("Unit quaternion", is_normalized());
    }

    /// Set this to the quaternion representation of the rotation matrix rot
    void from_rotation_matrix(const matrix_t& rot) {
      // Ken Shoemake, "Quaternion Calculus and Fast Animation", course notes, 
      // SIGGRAPH, 1987
      value_t s = rot(0,0)+rot(1,1)+rot(2,2);
      if (is_positive(s)) {
	s = std::sqrt(s+scalar_math<value_t>::one());
	(*this)[3] = s * static_cast<value_t>(0.5);
	s = static_cast<value_t>(0.5)/s;
	(*this)[0] = (rot(2,1)-rot(1,2))*s;
	(*this)[1] = (rot(0,2)-rot(2,0))*s;
	(*this)[2] = (rot(1,0)-rot(0,1))*s;
      } else {
	int i=0, j=1, k=2;
	if (rot(1,1) > rot(0,0)) {
	  i=1; j=2; k=0;
	}
	if (rot(2,2) > rot(i,i)) {
	  i=2; j=0; k=1;
	}
	s = std::sqrt(rot(i,i)-rot(j,j)-rot(k,k)+scalar_math<value_t>::one());
	(*this)[i]= s*static_cast<value_t>(0.5);
	s=static_cast<value_t>(0.5)/s;
	(*this)[3] = (rot(k,j)-rot(j,k))*s;
	(*this)[j] = (rot(j,i)+rot(i,j))*s;
	(*this)[k] = (rot(k,i)+rot(i,k))*s;
      }

      SL_ENSURE("Unit quaternion", is_normalized());
    }      

    /// Set m to the rotation matrix representation of this
    void rotation_matrix_in(matrix_t& m) const {
      SL_REQUIRE("Unit quaternion", is_normalized());

      // Ken Shoemake, "Quaternion Calculus and Fast Animation", course notes, 
      // SIGGRAPH, 1987

      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();

      const value_t x2 = x() * two;
      const value_t y2 = y() * two;
      const value_t z2 = z() * two;

      const value_t xx2= x() * x2;  const value_t yy2= y() * y2;  const value_t zz2= z() * z2;
      const value_t wx2= w() * x2;  const value_t wy2= w() * y2;  const value_t wz2= w() * z2;
      const value_t xy2= x() * y2;  const value_t yz2= y() * z2;  const value_t zx2= z() * x2;

      m = 
	one - yy2 - zz2, xy2 - wz2,       zx2 + wy2,       zero,
	xy2 + wz2,       one - xx2 - zz2, yz2 - wx2,       zero,
	zx2 - wy2,       yz2 + wx2,       one - xx2 - yy2, zero,
	zero,            zero,            zero,            one;
    }

    /// Set axis and angle to the axis-angle representation of this
    void axis_angle_in(axis_t& axis, 
		       value_t& angle) const {
      SL_REQUIRE("Unit quaternion", is_normalized());

      const value_t l2 = sqr(x()) + sqr(y()) + sqr(z());
      if (is_positive(l2)) {
	angle = (two(value_t()) * 
		 std::acos(median(w(), -one(value_t()), one(value_t()))));
	const value_t one_over_l = reciprocal_sqrt(l2);

	axis  = x() * one_over_l, y() * one_over_l, z() * one_over_l;
      } else {
	angle = zero(angle);
	axis.to_unit(1);
      }
    }

    /// Set axis to the i-th axes of the rotated frame
    void frame_axis_in(size_t i, axis_t& axis) const {
      SL_REQUIRE("Good index", i<3);
      SL_REQUIRE("Unit quaternion", is_normalized());

      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();

      switch (i) {
      case 0:
	{      
	  const value_t x2 = x() * two;
	  const value_t y2 = y() * two;
	  const value_t z2 = z() * two;
	  const value_t yy2= y() * y2;  const value_t zz2= z() * z2;
	  const value_t wy2= w() * y2;  const value_t wz2= w() * z2;
	  const value_t xy2= x() * y2;  const value_t zx2= z() * x2;
	  
	  axis = one - yy2 - zz2, xy2 + wz2, zx2 - wy2;
	}
	break;
      case 1:
	{
	  const value_t x2 = x() * two;
	  const value_t y2 = y() * two;
	  const value_t z2 = z() * two;
	  const value_t xx2= x() * x2;  const value_t zz2= z() * z2;
	  const value_t wx2= w() * x2;  const value_t wz2= w() * z2;
	  const value_t xy2= x() * y2;  const value_t yz2= y() * z2;

	  axis = xy2 - wz2, one - xx2 - zz2, yz2 + wx2;
	}
	break;
      case 2:
	{
	  const value_t x2 = x() * two;
	  const value_t y2 = y() * two;
	  const value_t z2 = z() * two;
	  const value_t xx2= x() * x2;  const value_t yy2= y() * y2;
	  const value_t wx2= w() * x2;  const value_t wy2= w() * y2;
	  const value_t yz2= y() * z2;  const value_t zx2= z() * x2;

	  axis = zx2 + wy2, yz2 - wx2, one - xx2 - yy2;
	}
	break;
      default:
	SL_CHECK("Wrong axis", i<3);
      }
    }
    
  public: // Converting constructors

    /// Init to the quaternion representation of a rotation of angle around axis
    inline quaternion(const axis_t& axis, 
		      const value_t& angle): super_t(tags::not_initialized()) {
      SL_REQUIRE("Unit axis", sl::abs(axis.two_norm_squared() - sl::one(value_t())) < static_cast<value_t>(0.001));
      from_axis_angle(axis,angle);
    }                  

    /// Init to the quaternion representation of the rotation matrix rot
    inline explicit quaternion(const matrix_t& rot): super_t(tags::not_initialized()) {
      from_rotation_matrix(rot);
    }                  
    
  public: 

    /// To identity
    inline void to_identity() {
      (*this)[0] = scalar_math<value_t>::zero();
      (*this)[1] = scalar_math<value_t>::zero();
      (*this)[2] = scalar_math<value_t>::zero();
      (*this)[3] = scalar_math<value_t>::one();
    }

    /// The identity quaternion
    static inline self_t identity() {
      return self_t(scalar_math<value_t>::zero(),
		    scalar_math<value_t>::zero(),
		    scalar_math<value_t>::zero(),
		    scalar_math<value_t>::one());
    }

  public: 

    /// The axis of rotation
    inline axis_t axis() const {
      return axis_t(this->storage_[0],this->storage_[1],this->storage_[2]).ok_normalized();
    }

    /// The angle of rotation 
    inline value_t angle() const {
      return 
	(two(value_t()) * 
	 std::acos(median((*this)[3], -one(value_t()), one(value_t())))); // Avoid numerical err.
    }

  public: // Division algebra operations

    /// The multiplication of this by other, i.e. the composition of rotations.
    self_t operator*(const self_t& other) const {
      return self_t((*this)[0]*other[3] + (*this)[3]*other[0] - (*this)[2]*other[1] + (*this)[1]*other[2],
		    (*this)[1]*other[3] + (*this)[2]*other[0] + (*this)[3]*other[1] - (*this)[0]*other[2],
		    (*this)[2]*other[3] - (*this)[1]*other[0] + (*this)[0]*other[1] + (*this)[3]*other[2],
		    (*this)[3]*other[3] - (*this)[0]*other[0] - (*this)[1]*other[1] - (*this)[2]*other[2]);
    }

    /// Post-multiply this by other
    inline self_t& operator *= (const self_t& other) {
      *this = *this * other;
      return static_cast<self_t&>(*this);
    }
 
    /// Invert to mo, setting ok to true if this is invertible
    void invert_to(self_t& mo, bool* ok) const {
      mo[0] = -(*this)[0];
      mo[1] = -(*this)[1];
      mo[2] = -(*this)[2];
      mo[3] =  (*this)[3];
      *ok = true;
    }

    /// The inverse of this.
    inline self_t inverse() const {
      return self_t(-(*this)[0],-(*this)[1],-(*this)[2],(*this)[3]);
    }

    /// The inverse of this.
    inline self_t operator~() const {
      return inverse();
    }                                                                                                                    
    
    /// Multiply this by inverse of other.
    inline self_t& operator /= (const self_t& other) {
      *this *= other.inverse();
      return static_cast<self_t&>(*this);
    }


  public: // G++ Workarounds

    /// Scale this by scalar
    inline self_t& operator *= (const value_t& scalar) { 
      return super_t::operator *= (scalar); 
    }

    /// Divide this by scalat
    inline self_t& operator /= (const value_t& scalar) { 
      return super_t::operator /= (scalar); 
    }

  public: // Other operations

    /// The conjugate of this, i.e. a rotation about the same axis by 2*Pi-angle() 
    inline self_t conjugate() const {
      return self_t(axis(), angle()-scalar_math<value_t>::two()*scalar_math<value_t>::Pi());
    }

    /// This normalized to unit length.
    self_t ok_normalized() const {
      self_t vr;
      bool ok;
      this->normalize_to(vr, &ok);
      if (!ok) {
	vr = identity(); // junk
      }
      return vr;
    }

    /// The closest to q between this and its conjugate
    self_t closest(const self_t& q) const {
      self_t result = (*this);
      self_t dq = ~q * result;
      self_t conj = conjugate();
      self_t dq2 = ~q * conj;
      if (dq.angle() > dq2.angle()) {
	result = conj;
      }
      return result;
    }

  public: // Linear interpolation

    /// Linear spherical interpolation
    self_t lerp(const self_t& qq2, const value_t& t) const {
      self_t q2 = qq2.closest((*this)); // HACK ### CHECK LATEST VERSION ###
      value_t q1_norm = (*this).two_norm();
      value_t q2_norm = q2.two_norm();
      value_t proj = median((*this).dot(q2)/(q1_norm*q2_norm),-sl::one(value_t()),sl::one(value_t()));

      value_t theta = std::acos(proj);

      value_t alpha, beta;
      if (sl::abs(theta) <= std::numeric_limits<value_t>::epsilon()) {
	alpha = sl::one(t)-t;
	beta  = t;
      } else {
	value_t sintheta = std::sin(theta);
	alpha = std::sin((sl::one(t)-t)*theta)/sintheta;
	beta  = std::sin(            t *theta)/sintheta;
      }
      self_t q(alpha*(*this)[0]+beta*q2[0],
	       alpha*(*this)[1]+beta*q2[1],
	       alpha*(*this)[2]+beta*q2[2],
	       alpha*(*this)[3]+beta*q2[3]);
      
      q = q.ok_normalized();
      
      return q;
    }

  }; // class quaternion

  template <typename OUT_ET>
  class conv_to< quaternion<OUT_ET> > {
  public:
    typedef quaternion<OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_array<4, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<4; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }
    
    // Explicit conversion from vectors of another type
    template <enum vector_orientation IN_ORIENTATION, typename IN_ET> 
    inline static result_t from(const fixed_size_vector<IN_ORIENTATION, 4, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<4; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }
    
    // Explicit conversion from points of another type
    template <typename IN_ET> 
    inline static result_t from(const quaternion<IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<4; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }    
  }; // class conv_to

} // namespace sl

// Arithmetic operators overloads

SL_OP_DIVISION_ALGEBRA_OVERLOADS(SL_OP_TEMPLATE1(template <class T>),
				 SL_OP_TEMPLATE1(sl::quaternion<T>),
				 T);

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  ///  A quaternion with single precision floating point components.
  typedef quaternion<float> quaternionf;
  
  /// A quaternion with double precision floating point components.
  typedef quaternion<double> quaterniond; 
}

#endif
