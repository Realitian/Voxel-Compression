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
#ifndef SL_PROJECTIVE_MAP_HPP
#define SL_PROJECTIVE_MAP_HPP

#include <sl/conv_to.hpp>
#include <sl/linear_map.hpp>

// --------------------------------------------------------------------
// -- sl::projective_map<DIMENSION,T>
// --------------------------------------------------------------------

namespace sl {

  /**
   *  An N-dimensional projective map
   */
  template <size_t N_dimension, class T>
  class projective_map: 
    public linear_map_base< general_tag, N_dimension, T> {
  public: // Constants and types
    typedef linear_map_base< general_tag, N_dimension, T> super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::concrete_t concrete_t;
    typedef typename super_t::value_t  value_t;
    typedef typename super_t::vector_t vector_t;
    typedef typename super_t::point_t  point_t;
    typedef typename super_t::plane_t  plane_t;
    typedef typename super_t::matrix_t matrix_t;

    enum { is_identity   = super_t::is_identity };
    enum { is_rigid_body = super_t::is_rigid_body };
    enum { is_affine     = super_t::is_affine };

  public: // Creation, Copy & Destruction

    /// Default creation (identity)
    inline projective_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Creator without initialization, handle with care!
    inline projective_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Construct this from matrix m.
    inline projective_map(const matrix_t& m): super_t(m) {
      SL_REQUIRE("Compatible matrix", is_compatible(m));
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    template <class T_other_tag >
    inline projective_map(const linear_map_base<
			  T_other_tag,
			  dimension, 
			  value_t>& other): super_t(other.storage()) {
      SL_INVARIANT(this->invariant());
    }

    /// Assign from compatible map
    template <class T_other_tag >
    inline concrete_t& operator = (const linear_map_base<
			           T_other_tag,
				   dimension, 
				   value_t>& other) {
      this->storage_ = other.storage();
      SL_INVARIANT(this->invariant());
      return *this;
    }

  };
    
  template <std::size_t DIM, typename OUT_ET>
  class conv_to< projective_map<DIM, OUT_ET> > {
  public:
    typedef projective_map<DIM, OUT_ET> result_t;
    
    // Explicit conversion from maps of another type
    template <typename IN_TAG, typename IN_ET>
    inline static result_t from(const linear_map_base<IN_TAG,
						      DIM,
						      IN_ET>& in) {
      return result_t(conv_to<typename result_t::matrix_t>::from(in.as_matrix()));
    }
    
  }; // class conv_to

}

// --------------------------------------------------------------------
// -- sl::projective_map<3,T>
// --------------------------------------------------------------------

namespace sl {

  /**
   * A 3D projective map (general case specialization)
   */
  template <class T>
  class projective_map<3, T>: 
    public linear_map_base< general_tag, 3, T> {
  public: // Constants and types
    typedef linear_map_base< general_tag, 3, T> super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::concrete_t concrete_t;
    typedef typename super_t::value_t  value_t;
    typedef typename super_t::vector_t vector_t;
    typedef typename super_t::point_t  point_t;
    typedef typename super_t::plane_t  plane_t;
    typedef typename super_t::matrix_t matrix_t;

    enum { is_identity   = super_t::is_identity };
    enum { is_rigid_body = super_t::is_rigid_body };
    enum { is_affine     = super_t::is_affine };

  public: // Creation, Copy & Destruction

    /// Default creator (identity)
    inline projective_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Creator without initialization, handle with care!
    inline projective_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Construct this from matrix m.
    inline projective_map(const matrix_t& m): super_t(m) {
      SL_REQUIRE("Compatible matrix", this->is_compatible(m));
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    template <class T_other_tag >
    inline projective_map(const linear_map_base<
			  T_other_tag,
			  dimension, 
			  value_t>& other): super_t(other.storage()) {
      SL_INVARIANT(this->invariant());
    }

    /// Assign from compatible map
    template <class T_other_tag >
    inline concrete_t& operator = (const linear_map_base<
			           T_other_tag,
				   dimension, 
				   value_t>& other) {
      (this->storage_) = other.storage();
      SL_INVARIANT(this->invariant());
      return *this;
    }

  public: // Specific 4D matrix features

    /// Helper: Compute a 3 by 3 determinant.
    inline value_t det(const value_t& m00,const value_t& m01,const value_t& m02,
                       const value_t& m10,const value_t& m11,const value_t& m12,
                       const value_t& m20,const value_t& m21,const value_t& m22) const {
      return
        m00*m11*m22 + m01*m12*m20 + m02*m10*m21 -
        m20*m11*m02 - m10*m01*m22 - m00*m21*m12;
    }

    /// The i-th OpenGL style clip plane (coord = +-1)
    plane_t clip_plane(std::size_t i) const {
      SL_REQUIRE("Good index", i<6);
      std::size_t ii = i/2;
      const matrix_t& P = this->as_matrix();
      plane_t result = tags::not_initialized();
      if (i%2 == 0) {
        // Left / Bottom / Near
        result[0] = P(3,0) + P(ii,0);
        result[1] = P(3,1) + P(ii,1);
        result[2] = P(3,2) + P(ii,2);
        result[3] = P(3,3) + P(ii,3);
      } else {
        // Right / Up / Far
        result[0] = P(3,0) - P(ii,0);
        result[1] = P(3,1) - P(ii,1);
        result[2] = P(3,2) - P(ii,2);
        result[3] = P(3,3) - P(ii,3);
      }
      result.normalize();
      return result;
    }

    /// The eye position
    point_t eye() const {
      // The 3 lines of the matrix are the normals to the planes x=0, y=0, z=0
      // in the camera CS. As we normalize them, we do not need the 4th coordinate.
      const matrix_t& P = this->as_matrix();
      
      vector_t line_0(P(0,0),P(0,1),P(0,2));
      vector_t line_1(P(1,0),P(1,1),P(1,2));
      vector_t line_2(P(2,0),P(2,1),P(2,2));

      line_0 = line_0.ok_normalized();
      line_1 = line_1.ok_normalized();
      line_2 = line_2.ok_normalized();
  
      // The camera position is at (0,0,0) in the camera CS so it is the
      // intersection of the 3 planes. It can be seen as the kernel
      // of the 3x4 projection matrix. We calculate it through 4 dimensional
      // vectorial product. We go directly into 3D that is to say we directly
      // divide the first 3 coordinates by the 4th one.

      // We derive the 4 dimensional vectorial product formula from the
      // computation of a 4x4 determinant that is developped according to
      // its 4th column. This implies some 3x3 determinants.
      value_t w = -det(P(0,0),P(0,1),P(0,2),
                     P(1,0),P(1,1),P(1,2),
                     P(2,0),P(2,1),P(2,2));
      value_t x =  det(P(0,1),P(0,2),P(0,3),
                     P(1,1),P(1,2),P(1,3),
                     P(2,1),P(2,2),P(2,3));
      value_t y = -det(P(0,0),P(0,2),P(0,3),
                     P(1,0),P(1,2),P(1,3),
                     P(2,0),P(2,2),P(2,3));
      value_t z =  det(P(0,0),P(0,1),P(0,3),
                     P(1,0),P(1,1),P(1,3),
                     P(2,0),P(2,1),P(2,3));

      return point_t(x/w, y/w, z/w);
    }

    /// The field of view in the x direction
    value_t fovx() const {
      point_t pnl = inverse_transformation(*this, point_t(-1.0f, 0.0f, -1.0f));
      point_t pnr = inverse_transformation(*this, point_t( 1.0f, 0.0f, -1.0f));
      point_t pfl = inverse_transformation(*this, point_t(-1.0f, 0.0f,  1.0f));
      point_t pfr = inverse_transformation(*this, point_t( 1.0f, 0.0f,  1.0f));
      vector_t l = pfl - pnl;
      vector_t r = pfr - pnr;
      return l.angle(r);
    }

    /// The field of view in the y direction
    value_t fovy() const {
      point_t pnl = inverse_transformation(*this, point_t(0.0f,-1.0f, -1.0f));
      point_t pnh = inverse_transformation(*this, point_t(0.0f, 1.0f, -1.0f));
      point_t pfl = inverse_transformation(*this, point_t(0.0f,-1.0f,  1.0f));
      point_t pfh = inverse_transformation(*this, point_t(0.0f, 1.0f,  1.0f));
      vector_t l = pfl - pnl;
      vector_t h = pfh - pnh;
      return l.angle(h);
    }

    sl::fixed_size_array<6,value_t> perspective_frustum_lrbtnf() const {
      typedef sl::fixed_size_vector<sl::column_orientation,4,value_t> vector4_t;

      sl::fixed_size_array<6,value_t> result;

      const sl::fixed_size_square_matrix<4,value_t> Pinv = (~(*this)).as_matrix();
      vector4_t lbnw = Pinv * vector4_t(-1,-1,-1, 1);
      vector4_t rtnw = Pinv * vector4_t( 1, 1,-1, 1);
      vector4_t rtfw = Pinv * vector4_t( 1, 1, 1, 1);
      result[0] = lbnw[0]/lbnw[3]; // l
      result[1] = rtnw[0]/rtnw[3]; // r
      result[2] = lbnw[1]/lbnw[3]; // b
      result[3] = rtnw[1]/rtnw[3]; // t
      result[4] = -lbnw[2]/lbnw[3]; // n
      result[5] = -rtfw[2]/rtfw[3]; // f
      
      return result;
    }

    sl::fixed_size_array<6,value_t> ortho_frustum_lrbtnf() const {
      typedef sl::fixed_size_vector<sl::column_orientation,4,value_t> vector4_t;

      sl::fixed_size_array<6,value_t> result;

      const sl::fixed_size_square_matrix<4,value_t> Pinv = (~(*this)).as_matrix();
      vector4_t lbnw = Pinv * vector4_t(-1,-1,-1, 1);
      vector4_t rtnw = Pinv * vector4_t( 1, 1,-1, 1);
      vector4_t rtfw = Pinv * vector4_t( 1, 1, 1, 1);
      result[0] = lbnw[0]/lbnw[3]; // l
      result[1] = rtnw[0]/rtnw[3]; // r
      result[2] = lbnw[1]/lbnw[3]; // b
      result[3] = rtnw[1]/rtnw[3]; // t
      result[4] = lbnw[2]/lbnw[3]; // n
      result[5] = rtfw[2]/rtfw[3]; // f
      
      return result;
    }

  public: // Standard 3D perspective (obsolete!)
    
    /// True iff this represents a standard 3D perspective (possibly off center)
    inline bool is_std_3d_perspective() const {
      return 
	((*this)(1,0) == sl::zero(value_t())) && ((*this)(2,0) == sl::zero(value_t())) && ((*this)(3,0) ==  sl::zero(value_t())) &&
	((*this)(0,1) == sl::zero(value_t())) && ((*this)(2,1) == sl::zero(value_t())) && ((*this)(3,1) ==  sl::zero(value_t())) &&
	/* ((*this)(0,2) == sl::zero(value_t())) && ((*this)(1,2) == sl::zero(value_t())) && */ ((*this)(3,2) == -sl::one(value_t())) &&
	((*this)(0,3) == sl::zero(value_t())) && ((*this)(1,3) == sl::zero(value_t())) && ((*this)(3,3) ==  sl::zero(value_t()));      
    }

    /// True iff this represents a standard 3D perspective (possibly off center)
    inline bool is_std_centered_3d_perspective() const {
      return
	is_std_3d_perspective() &&
	((*this)(0,2) == sl::zero(value_t())) &&
	((*this)(1,2) == sl::zero(value_t()));
    }
 
    /// The field of view of this (requires is_std_3d_perspective)
    inline value_t fov_from_std_3d_perspective() const {
      SL_REQUIRE("Is standard perspective", is_std_3d_perspective());
      return sl::two(value_t()) * std::atan(sl::one(value_t())/(*this)(1,1));
    }
    
    /// The aspect of this (requires is_std_3d_perspective)
    inline value_t aspect_from_std_3d_perspective() const {
      SL_REQUIRE("Is standard perspective", is_std_3d_perspective());
      return (*this)(1,1)/(*this)(0,0);
    }
      
    /// The near plane of this (requires is_std_3d_perspective)
    inline value_t near_from_std_3d_perspective() const {
      SL_REQUIRE("Is standard perspective", is_std_3d_perspective());
      return (*this)(2,3)/(-sl::one(value_t())+(*this)(2,2));
    }
 
    /// The far plane of this (requires is_std_3d_perspective)
    inline value_t far_from_std_3d_perspective() const {
      SL_REQUIRE("Is standard perspective", is_std_3d_perspective());
      return (*this)(2,3)/(sl::one(value_t())+(*this)(2,2));
    }                                                                                                                    
      
  };
  
  // ---------------------------------------------------------------------------
  // Distortion
  // ---------------------------------------------------------------------------

  /**
   *  Perform projection using a second order model of radial distortion.
   */
  template <class T>
  fixed_size_vector<column_orientation,4,T> barrel_distorted_projection(const projective_map<3,T>& P,
									const T& k1,
									const T& k2,
									const fixed_size_vector<column_orientation,4,T>& p) {
    SL_REQUIRE("Is standard perspective", P.is_std_3d_perspective());
    
    typedef T value_t;
    typedef fixed_size_vector<column_orientation,4,T> hvector_t;
    
    const value_t fx = P(0,0);
    const value_t fy = P(1,1);
    const value_t cx = P(0,2);
    const value_t cy = P(1,2);
    const value_t zs = P(2,2);
    const value_t zt = P(2,3);

    value_t eta = value_t(1.0);
    if (p[2]!=value_t(0.0)) {
      // Point is at finite distance, compute distortion
      const float R2 = (p[0]*p[0]+p[1]*p[1])/(p[2]*p[2]);
      eta += k1 * R2 + k2 * R2 * R2;
    }

    hvector_t result = tags::not_initialized();
    result[0] = eta * fx * p[0] + cx * p[2];
    result[1] = eta * fy * p[1] + cy * p[2];
    result[2] =       zs * p[2] + zt * p[3];
    result[3] =          - p[2];
    
    return result;	
  }

} // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  /// A 2D projective map with single precision floating point components.
  typedef projective_map<2,float> projective_map2f;
  /// A 3D projective map with single precision floating point components.
  typedef projective_map<3,float> projective_map3f;
  /// A 4D projective map with single precision floating point components.
  typedef projective_map<4,float> projective_map4f;

  /// A 2D projective map with double precision floating point components.
  typedef projective_map<2,double> projective_map2d;
  /// A 3D projective map with double precision floating point components.
  typedef projective_map<3,double> projective_map3d;
  /// A 4D projective map with double precision floating point components.
  typedef projective_map<4,double> projective_map4d;
}

#endif





