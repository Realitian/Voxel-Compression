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
#ifndef SL_AFFINE_MAP_HPP
#define SL_AFFINE_MAP_HPP

#include <sl/conv_to.hpp>
#include <sl/linear_map.hpp>
#include <sl/quaternion.hpp>

// --------------------------------------------------------------------
// sl::affine_map<DIMENSION,T>
// --------------------------------------------------------------------

namespace sl {

  /**
   *  An N-dimensional affine maps
   */
  template <size_t N_dimension, class T>
  class affine_map: 
    public linear_map_base<affine_tag, N_dimension, T>
  {
  public: // Constants and types
    typedef linear_map_base<affine_tag, N_dimension, T> super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::concrete_t    concrete_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;
    typedef typename super_t::point_t       point_t;
    typedef typename super_t::plane_t       plane_t;
    typedef typename super_t::matrix_t      matrix_t;

    enum { is_identity   = super_t::is_identity };
    enum { is_rigid_body = super_t::is_rigid_body };
    enum { is_affine     = super_t::is_affine };

  public: // Creation, Copy & Destruction

    /// Default init (identity)
    inline affine_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Fast init (garbage), handle with care!
    inline affine_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from matrix m, assuming m is affine
    inline affine_map(const matrix_t& m): super_t(m) {
      SL_REQUIRE("Compatible matrix", this->is_compatible(m));
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    template <class T_other_tag >
    inline affine_map(const linear_map_base<
		      T_other_tag,
		      dimension, 
		      value_t>& other): super_t(other.storage()) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_affine);
      SL_INVARIANT(this->invariant());
    }

    /// Assign from compatible map
    template <class T_other_tag >
    inline concrete_t& operator = (const linear_map_base<
			           T_other_tag,
				   dimension, 
				   value_t>& other) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_affine);
      this->storage_ = other.storage();
      SL_INVARIANT(this->invariant());
      return *this;
    }

  };
  
  template <std::size_t DIM, typename OUT_ET>
  class conv_to< affine_map<DIM, OUT_ET> > {
  public:
    typedef affine_map<DIM, OUT_ET> result_t;
    
    // Explicit conversion from maps of another type
    template <typename IN_TAG, typename IN_ET>
    inline static result_t from(const linear_map_base<IN_TAG,
						      DIM,
						      IN_ET>& in) {
      return result_t(conv_to<typename result_t::matrix_t>::from(in.as_matrix()));
    }
    
  }; // class conv_to
  
} // namespace sl

// --------------------------------------------------------------------
// sl::affine_map<3,T>
// --------------------------------------------------------------------

namespace sl {

  /// 3D affine maps
  template <class T>
  class affine_map<3,T>: 
    public linear_map_base<affine_tag, 3, T> {
  public: // Constants and types

    typedef linear_map_base<affine_tag, 3, T> super_t;
    enum { dimension = super_t::dimension };

    typedef typename super_t::concrete_t    concrete_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;
    typedef typename super_t::point_t       point_t;
    typedef typename super_t::plane_t       plane_t;
    typedef typename super_t::matrix_t      matrix_t;
    typedef quaternion<T>                   quaternion_t;

    enum { is_identity   = super_t::is_identity };
    enum { is_rigid_body = super_t::is_rigid_body };
    enum { is_affine     = super_t::is_affine };

  public: // Creation, Copy & Destruction

    /// Default init, identity
    inline affine_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Fast init (garbage), handle with care!
    inline affine_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from matrix m, assuming m is affine 
    inline affine_map(const matrix_t& m): super_t(m) {
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    template <class T_other_tag >
    inline affine_map(const linear_map_base<
		      T_other_tag,
		      dimension, 
		      value_t>& other): super_t(other.storage()) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_affine);
      SL_INVARIANT(this->invariant());
    }

    /// Assign from compatible map
    template <class T_other_tag >
    inline concrete_t& operator = (const linear_map_base<
			           T_other_tag,
				   dimension, 
				   value_t>& other) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_affine);
      this->storage_ = other.storage();
      SL_INVARIANT(this->invariant());
      return *this;
    }

    /// Init from rotation specified by quaternion q
    affine_map(const quaternion_t& q): super_t(tags::not_initialized()) {
      quaternion_t qt = q.ok_normalized();

      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();

      value_t xx2= qt[0] * qt[0] * two;  value_t yy2= qt[1] * qt[1] * two;  value_t zz2= qt[2] * qt[2] * two; 
      value_t wx2= qt[3] * qt[0] * two;  value_t wy2= qt[3] * qt[1] * two;  value_t wz2= qt[3] * qt[2] * two; 
      value_t xy2= qt[0] * qt[1] * two;  value_t yz2= qt[1] * qt[2] * two;  value_t zx2= qt[2] * qt[0] * two; 
  
      value_t m11= one - yy2 - zz2;   value_t m21= xy2 + wz2;        value_t m31= zx2 - wy2;
      value_t m12= xy2 - wz2;         value_t m22= one - xx2 - zz2;  value_t m32= yz2 + wx2;
      value_t m13= zx2 + wy2;         value_t m23= yz2 - wx2;        value_t m33= one - xx2 - yy2;
#if 0
      value_t n1= std::sqrt(sqr(m11)+sqr(m12)+sqr(m13));
      value_t n2= std::sqrt(sqr(m21)+sqr(m22)+sqr(m23));
      value_t n3= std::sqrt(sqr(m31)+sqr(m32)+sqr(m33));
    
      this->storage_ = matrix_t(m11/n1, m12/n1, m13/n1, zero,
			  m21/n2, m22/n2, m23/n2, zero,
			  m31/n3, m32/n3, m33/n3, zero,
			  zero,     zero,   zero,  one);
#else
      this->storage_ = matrix_t( m11, m12, m13, zero,
			   m21, m22, m23, zero,
			   m31, m32, m33, zero,
			  zero,zero,zero, one);
#endif
      SL_INVARIANT(this->invariant());
    } 

    /// Init from components 
    affine_map(const vector_t& scaling_value,
	       const vector_t& shearing_value,
	       const quaternion_t& rotation_value,
	       const vector_t &translation_value): super_t(tags::not_initialized()) {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();

      value_t sx, sy, sz;
      value_t sxy, sxz, syz;

      sx  = scaling_value[0];   sy  = scaling_value[1];  sz  = scaling_value[2];
      sxy = shearing_value[0]; sxz = shearing_value[1]; syz = shearing_value[2];

      this->storage_ = 
	matrix_t(one, zero, zero,  translation_value[0],
		 zero, one, zero,  translation_value[1],
		 zero, zero, one,  translation_value[2],
		 zero, zero, zero, one) *
	affine_map(rotation_value).as_matrix() *
	matrix_t(sx,       zero, zero, zero,
		 sx*sxy,   sy,   zero, zero,
		 sx*sxz,   sy*syz, sz, zero,
		 zero,     zero, zero, one);
      SL_INVARIANT(this->invariant());
    }

  public: // Additional features

    /// rotation part converted to quaternion
    quaternion_t as_quaternion() const {
      quaternion_t result;

      value_t s = (*this)(0,0)+(*this)(1,1)+(*this)(2,2);
      if (s  >= std::numeric_limits<value_t>::epsilon()) {
	s = std::sqrt(s+sl::scalar_math<value_t>::one());
	result[3] = s * static_cast<value_t>(0.5);
	s = static_cast<value_t>(0.5)/s;
	result[0] = ((*this)(2,1)-(*this)(1,2))*s;
	result[1] = ((*this)(0,2)-(*this)(2,0))*s;
	result[2] = ((*this)(1,0)-(*this)(0,1))*s;
      } else {
	int i=0, j=1, k=2;
	if ((*this)(1,1) > (*this)(0,0)) {
	  i=1; j=2; k=0;
	}
	if ((*this)(2,2) > (*this)(i,i)) {
	  i=2; j=0; k=1;
	}
	s = std::sqrt((*this)(i,i)-(*this)(j,j)-(*this)(k,k)+sl::scalar_math<value_t>::one());
	result[i]= s*static_cast<value_t>(0.5);
	s=static_cast<value_t>(0.5)/s;
	result[3] = ((*this)(k,j)-(*this)(j,k))*s;
	result[j] = ((*this)(j,i)+(*this)(i,j))*s;
	result[k] = ((*this)(k,i)+(*this)(i,k))*s;
      }
      return result.ok_normalized();
    }

    /// rotation part converted to Euler angles
    vector_t as_euler_angles() const {
      value_t cx, cy, cz;
      value_t sx, sy, sz;

      sy = -(*this)(2,0);
      cy = sl::scalar_math<value_t>::one() - (sy * sy);
      if (cy > std::numeric_limits<value_t>::epsilon()) {
	cy = std::sqrt(cy);
	cx =  (*this)(2,2)/cy;
	sx =  (*this)(2,1)/cy;
	cz =  (*this)(0,0)/cy;
	sz =  (*this)(1,0)/cy;
      } else {
	cy = sl::scalar_math<value_t>::zero();
	cx =  (*this)(1,1);
	sx = -(*this)(1,2);
	cz = sl::scalar_math<value_t>::one();
	sz = sl::scalar_math<value_t>::zero();
      }
      return vector_t(atan2(sx, cx), atan2(sy, cy), atan2(sz, cz)); 
    } 

    /// Factorize this to components
    void factorize_to(vector_t& scaling_value,
		      vector_t& shearing_value,
		      quaternion_t& rotation_value,
		      vector_t &translation_value) const {
      const value_t one    = sl::scalar_math<value_t>::one();
      const value_t zero   = sl::scalar_math<value_t>::zero();
      if (is_identity) {
	translation_value    = vector_t::zero();
	rotation_value       = quaternion_t::identity();
	scaling_value        = vector_t(one,one,one);
	shearing_value       = vector_t::zero();
      } else {
	value_t xx, xy, xz;
	value_t yx, yy, yz;
	value_t zx, zy, zz;
	value_t sx, sy, sz;
	value_t sxy, sxz, syz;

	translation_value = vector_t((*this)(0,3), (*this)(1,3), (*this)(2,3));

	xx = (*this)(0,0); xy = (*this)(1,0); xz = (*this)(2,0);
	yx = (*this)(0,1); yy = (*this)(1,1); yz = (*this)(2,1);
	zx = (*this)(0,2); zy = (*this)(1,2); zz = (*this)(2,2);
	
	sx = std::sqrt(sqr(xx)+sqr(xy)+sqr(xz));
	xx = xx/sx; xy = xy/sx; xz = xz/sx;
	
	sxy = xx*yx+xy*yy+xz*yz;
	yx = yx-sxy*xx;
	yy = yy-sxy*xy;
	yz = yz-sxy*xz;
    
	sy = std::sqrt(sqr(yx)+sqr(yy)+sqr(yz));
	yx = yx/sy; yy = yy/sy; yz = yz/sy;
	sxy = sxy/sy;
	
	sxz= xx*zx+xy*zy+xz*zz;
	syz= yx*zx+yy*zy+yz*zz;
	zx = zx-sxz*xx-syz*yx;
	zy = zy-sxz*xy-syz*yy;
	zz = zz-sxz*xz-syz*yz;
    
	sz= std::sqrt(sqr(zx)+sqr(zy)+sqr(zz));
	zx= zx/sz; 
	zy= zy/sz;
	zz= zz/sz;
	sxz= sxz/sz;
	syz= syz/sz;
	
	scaling_value  = vector_t(sx, sy, sz);
	shearing_value = vector_t(sxy, sxz, syz);
	matrix_t m = matrix_t(xx, yx, zx, zero, 
			      xy, yy, zy, zero,
			      xz, yz, zz, zero,
			      zero, zero, zero, one);
	if (is_negative(m.determinant())) {
	  m *= -scalar_math<value_t>::one(); m(3,3) = scalar_math<value_t>::one();
	  scaling_value = scaling_value;
	}
	
	rotation_value = concrete_t(m).as_quaternion();
      }
    } 

    /// Linear interpolation (from factorization)
    concrete_t lerp(const concrete_t& f2, value_t t) const {
      vector_t p1, p2;
      quaternion_t q1, q2;
      vector_t sh1, sh2, sc1, sc2;
      
      (*this).factorize_to(sc1,sh1,q1,p1);
      f2.factorize_to(sc2,sh2,q2,p2);
      return concrete_t(sc1.lerp(sc2,t),
                        sh1.lerp(sh2,t),
                        q1.lerp(q2,t),
                        p1.lerp(p2,t));
    }
      
  };

  
} // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  /// A 2D affine map with single precision floating point components.
  typedef affine_map<2,float> affine_map2f;
  /// A 3D affine map with single precision floating point components.
  typedef affine_map<3,float> affine_map3f;
  /// A 4D affine map with single precision floating point components.
  typedef affine_map<4,float> affine_map4f;

  /// A 2D affine map with double precision floating point components.
  typedef affine_map<2,double> affine_map2d;
  /// A 3D affine map with double precision floating point components.
  typedef affine_map<3,double> affine_map3d;
  /// A 4D affine map with double precision floating point components.
  typedef affine_map<4,double> affine_map4d;
}

#endif




