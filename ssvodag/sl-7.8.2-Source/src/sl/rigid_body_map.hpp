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
#ifndef SL_RIGID_BODY_MAP_HPP
#define SL_RIGID_BODY_MAP_HPP

#include <sl/conv_to.hpp>
#include <sl/linear_map.hpp>
#include <sl/quaternion.hpp>
#include <cassert>

// --------------------------------------------------------------------
// sl::rigid_body_map<DIMENSION,T>
// --------------------------------------------------------------------

namespace sl {

  /**
   *  N-dimensional rigid_body maps
   */
  template <size_t N_dimension, class T>
  class rigid_body_map: 
    public linear_map_base<sl::rigid_body_tag, N_dimension, T>
  {
  public: // Constants and types
    typedef linear_map_base<sl::rigid_body_tag, N_dimension, T> super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::concrete_t     concrete_t;
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
    inline rigid_body_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Fast init (garbage), handle with care!
    inline rigid_body_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from matrix m, assuming m is rigid_body
    inline rigid_body_map(const matrix_t& m): super_t(m) {
      SL_REQUIRE("Compatible matrix", this->is_compatible(m));
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    template <class T_other_tag >
    inline rigid_body_map(const linear_map_base<
			  T_other_tag,
			  dimension, 
			  value_t>& other): super_t(other.storage()) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_rigid_body);
      SL_INVARIANT(this->invariant());
    }

    /// Assign from compatible map
    template <class T_other_tag >
    inline concrete_t& operator = (const linear_map_base<
			           T_other_tag,
				   dimension, 
				   value_t>& other) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_rigid_body);
      this->storage_ = other.storage();
      SL_INVARIANT(this->invariant());
      return *this;
    }

  };

    
  template <std::size_t DIM, typename OUT_ET>
  class conv_to< rigid_body_map<DIM, OUT_ET> > {
  public:
    typedef rigid_body_map<DIM, OUT_ET> result_t;
    
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
// sl::rigid_body_map<3,T>
// --------------------------------------------------------------------

namespace sl {

  /// 3D rigid_body maps
  template <class T>
  class rigid_body_map<3,T>: 
    public linear_map_base<rigid_body_tag, 3, T> {
  public: // Constants and types

    typedef linear_map_base<rigid_body_tag, 3, T> super_t;

    enum { dimension = super_t::dimension };

    typedef quaternion<T>                   quaternion_t;

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

    /// Default init, identity
    inline rigid_body_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Fast init (garbage), handle with care!
    inline rigid_body_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from matrix m, assuming m is rigid_body 
    inline rigid_body_map(const matrix_t& m): super_t(m) {
      SL_REQUIRE("Compatible matrix", this->is_compatible(m));
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    template <class T_other_tag >
    inline rigid_body_map(const linear_map_base<
			  T_other_tag,
			  dimension, 
			  value_t>& other): super_t(other.storage()) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_rigid_body);
      SL_INVARIANT(this->invariant());
    }

    /// Assign from compatible map
    template <class T_other_tag >
    inline concrete_t& operator = (const linear_map_base<
			           T_other_tag,
				   dimension, 
				   value_t>& other) {
      SL_COMPILE_TIME_CHECK("Compatible", T_other_tag::is_rigid_body);
      this->storage_ = other.storage();
      SL_INVARIANT(this->invariant());
      return *this;
    }

    /// Init from rotation specified by quaternion q
    inline rigid_body_map(const quaternion_t& q): super_t(tags::not_initialized()) {
      q.rotation_matrix_in(this->storage_);
    }

    /// Init from components 
    rigid_body_map(const quaternion_t& rotation_value,
		   const vector_t &translation_value): super_t(tags::not_initialized()) {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();

      this->storage_ = 
        matrix_t(one, zero, zero,  translation_value[0],
                 zero, one, zero,  translation_value[1],
                 zero, zero, one,  translation_value[2],
                 zero, zero, zero, one) *
        rigid_body_map(rotation_value).as_matrix();
      SL_INVARIANT(this->invariant());
    }

    /// Init from translation
    rigid_body_map(const vector_t &translation_value): super_t(tags::not_initialized()) {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();

      this->storage_ = 
        matrix_t(one, zero, zero,  translation_value[0],
                 zero, one, zero,  translation_value[1],
                 zero, zero, one,  translation_value[2],
                 zero, zero, zero, one);
      SL_INVARIANT(this->invariant());
    }

  public: // Additional features

    /// rotation part converted to quaternion
    quaternion_t as_quaternion() const {
      quaternion_t result = tags::not_initialized();
      result.from_rotation_matrix(this->storage_);
      return result;
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
    void factorize_to(quaternion_t& rotation_value,
		      vector_t &translation_value) const {
      if (is_identity) {
        translation_value    = vector_t::zero();
        rotation_value       = quaternion_t::identity();
      } else {
        translation_value = vector_t((*this)(0,3), (*this)(1,3), (*this)(2,3));
        rotation_value = as_quaternion();
      }
    } 

    /// Linear interpolation (from factorization)
    concrete_t lerp(const concrete_t& f2, value_t t) const {
      vector_t p1, p2;
      quaternion_t q1, q2;
      
      (*this).factorize_to(q1,p1);
      f2.factorize_to(q2,p2);
      return concrete_t(q1.lerp(q2,t),
                        p1.lerp(p2,t));
    }
      
  };

}; // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  /// A 2D rigid_body map with single precision floating point components.
  typedef rigid_body_map<2,float> rigid_body_map2f;
  /// A 3D rigid_body map with single precision floating point components.
  typedef rigid_body_map<3,float> rigid_body_map3f;
  /// A 4D rigid_body map with single precision floating point components.
  typedef rigid_body_map<4,float> rigid_body_map4f;

  /// A 2D rigid_body map with double precision floating point components.
  typedef rigid_body_map<2,double> rigid_body_map2d;
  /// A 3D rigid_body map with double precision floating point components.
  typedef rigid_body_map<3,double> rigid_body_map3d;
  /// A 4D rigid_body map with double precision floating point components.
  typedef rigid_body_map<4,double> rigid_body_map4d;
}

#endif



