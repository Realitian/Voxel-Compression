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
#ifndef SL_LINEAR_MAP_FACTORY_HPP
#define SL_LINEAR_MAP_FACTORY_HPP

#include <sl/projective_map.hpp>
#include <sl/affine_map.hpp>
#include <sl/rigid_body_map.hpp>
#include <sl/utility.hpp>

namespace sl {

  /// Base class for objects that build N-dimensional linear maps
  template <size_t DIMENSION, class T>
  class linear_map_factory_base {
    
  public: // Constants and types 

    enum { dimension = DIMENSION };

    typedef projective_map<dimension,T>          projective_map_t;
    typedef affine_map<dimension, T>             affine_map_t;
    typedef rigid_body_map<dimension,T>          rigid_body_map_t;

    typedef typename affine_map_t::value_t       value_t;
    typedef typename affine_map_t::matrix_t      matrix_t;
    typedef typename affine_map_t::vector_t      vector_t;
    typedef typename affine_map_t::dual_vector_t dual_vector_t;
    typedef typename affine_map_t::point_t       point_t;
    typedef typename affine_map_t::plane_t       plane_t;

  public: // Creators

    /// Identity map
    static rigid_body_map_t identity() {
      rigid_body_map_t result;
      return result;
    }

    /// Translation
    static rigid_body_map_t translation(const vector_t& vv) {
      matrix_t m = tags::not_initialized();
      m.to_identity();
      for (size_t i = 0; i<dimension; i++) {
        m(i,dimension) = vv[i];
      }
      return rigid_body_map_t(m);
    } 
  
    /// Scaling
    static affine_map_t scaling(const vector_t& vv) {
      matrix_t m;
      for (size_t i = 0; i<dimension; i++) {
        m(i,i) = vv[i];
      }
      m(dimension,dimension) = scalar_math<value_t>::one();
      return affine_map_t(m);
    } 

    /// Uniform scaling
    static affine_map_t scaling(value_t s) {
      matrix_t m;
      for (size_t i = 0; i<dimension; i++) {
        m(i,i) = s;
      }
      m(dimension,dimension) = scalar_math<value_t>::one();
      return affine_map_t(m);
    } 
  }; // class linear_map_factory_base<DIMENSION,T>
      
}; // namespace sl


// ---------------------------------------------------------------------
// sl::linear_map_factory<DIMENSION,T>
// ---------------------------------------------------------------------

namespace sl {

  /// Objects that build N-dimensional affine maps
  template <size_t DIMENSION, class T>
  class linear_map_factory: 
    public linear_map_factory_base<DIMENSION,T> {
  public:
    typedef linear_map_factory_base<DIMENSION,T> super_t;

    typedef typename super_t::value_t           value_t;
    typedef typename super_t::projective_map_t  projective_map_t;
    typedef typename super_t::rigid_body_map_t  rigid_body_map_t;
    typedef typename super_t::affine_map_t      affine_map_t;
    typedef typename super_t::vector_t          vector_t;
    typedef typename super_t::dual_vector_t     dual_vector_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::plane_t           plane_t;
    typedef typename super_t::matrix_t          matrix_t;

  }; // class linear_map_factory<DIMENSION,T>

}; // namespace sl

// ---------------------------------------------------------------------
// sl::linear_map_factory<4,T>
// ---------------------------------------------------------------------

namespace sl {

  /// Objects that build 4D affine maps
  template <class T>
  class linear_map_factory<4,T>: 
    public linear_map_factory_base<4,T> {
  public:
    typedef linear_map_factory_base<4,T> super_t;

    typedef typename super_t::value_t           value_t;
    typedef typename super_t::projective_map_t  projective_map_t;
    typedef typename super_t::rigid_body_map_t  rigid_body_map_t;
    typedef typename super_t::affine_map_t      affine_map_t;
    typedef typename super_t::vector_t          vector_t;
    typedef typename super_t::dual_vector_t     dual_vector_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::plane_t           plane_t;
    typedef typename super_t::matrix_t          matrix_t;

  public:

    /// Identity
    static rigid_body_map_t identity() {
      return super_t::identity();
    }

    /// Translation
    static rigid_body_map_t translation(const vector_t& v) {
      return super_t::translation(v);
    }

    /// Scaling
    static affine_map_t scaling(const vector_t& vv) {
      return super_t::scaling(vv);
    }

    /// Uniform scaling
    static affine_map_t scaling(value_t s) {
      return super_t::scaling(s);
    }

    /**
     * A rigid body map that has the unit vector "u" as first axis.
     * John Hughes and Tomas Moeller, "Efficiently Building a
     * Matrix to Rotate One Vector to Another", ACM Journal of
     * Graphics Tools, 4(4): 33--35.
     */
    static rigid_body_map_t basis_from(const vector_t& u) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();

#ifndef NDEBUG
      const value_t eps = scalar_math<value_t>::epsilon();
      SL_REQUIRE("u is normalized", sl::abs(u.two_norm() - one) < 100.0f*eps);
#endif
      
      matrix_t m=tags::not_initialized();
      m = 
        u[0], -u[1], -u[2], -u[3], zero,
        u[1],  u[0],  u[3], -u[2], zero,
        u[2], -u[3],  u[0],  u[1], zero,
        u[3],  u[2], -u[1],  u[0], zero,
        zero,  zero,  zero,  zero, one;

      return rigid_body_map_t(m);
    } 
    
  }; // class linear_map_factory<4,T>

} // namespace sl

// ---------------------------------------------------------------------
// sl::linear_map_factory<3,T>
// ---------------------------------------------------------------------

namespace sl {

  /// Objects that build 3D affine maps
  template <class T>
  class linear_map_factory<3,T>: 
    public linear_map_factory_base<3,T> {
  public:
    typedef linear_map_factory_base<3,T> super_t;

    typedef typename super_t::value_t           value_t;
    typedef typename super_t::projective_map_t  projective_map_t;
    typedef typename super_t::rigid_body_map_t  rigid_body_map_t;
    typedef typename super_t::affine_map_t      affine_map_t;
    typedef typename super_t::vector_t          vector_t;
    typedef typename super_t::dual_vector_t     dual_vector_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::plane_t           plane_t;
    typedef typename super_t::matrix_t          matrix_t;
    typedef typename affine_map_t::quaternion_t quaternion_t;
  public:

    /// identity
    static rigid_body_map_t identity() {
      return super_t::identity();
    }

    /// Translation
    static rigid_body_map_t translation(const vector_t& v) {
      return super_t::translation(v);
    }

    /// Translation
    static rigid_body_map_t translation(value_t x,
                                        value_t y,
                                        value_t z) {
      return translation(vector_t(x,y,z));
    }

    /// Scaling
    static affine_map_t scaling(const vector_t& vv) {
      return super_t::scaling(vv);
    }

    /// Uniform scaling
    static affine_map_t scaling(value_t s) {
      return super_t::scaling(s);
    }

    /// Scaling
    static affine_map_t scaling(value_t x,
                                value_t y,
                                value_t z) {
      return scaling(vector_t(x,y,z));
    }

    /// Shear
    static affine_map_t shear(value_t xy, 
                              value_t xz, 
                              value_t yz) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();
      return affine_map_t(matrix_t(one, zero, zero, zero,
                                   xy, one, zero, zero,
                                   xz, yz, one, zero,
                                   zero, zero, zero, one));
    } 

    /**
     * A rigid body map that has the unit vector "u" as first axis.
     * John Hughes and Tomas Moeller, "Efficiently Building a
     * Matrix to Rotate One Vector to Another", ACM Journal of
     * Graphics Tools, 4(4): 33--35.
     */
    static rigid_body_map_t basis_from(const vector_t& u) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();

#ifndef NDEBUG
      const value_t eps = scalar_math<value_t>::epsilon();
      SL_REQUIRE("u is normalized", sl::abs(u.two_norm() - one) < 100.0f*eps);
#endif
      //const value_t x_abs = sl::abs(u[0]);
      //const value_t y_abs = sl::abs(u[1]);
      //const value_t z_abs = sl::abs(u[2]);

      // Zero out the smallest entry of u and swap the other two, negating the first one
      const size_t u_iamin0 = u.iamin();
      const size_t u_iamin1 = min((u_iamin0+1)%3, (u_iamin0+2)%3); 
      const size_t u_iamin2 = max((u_iamin0+1)%3, (u_iamin0+2)%3); 

      vector_t v_hat = tags::not_initialized();
      v_hat[u_iamin0] = zero;
      v_hat[u_iamin1] = -u[u_iamin2];
      v_hat[u_iamin2] =  u[u_iamin1];

      // Normalize and build third basis vector with cross product
      const vector_t v = v_hat / v_hat.two_norm(); // guaranteed not null
      const vector_t w = u.cross(v);
     
      matrix_t m=tags::not_initialized();

      m = 
        u[0], v[0], w[0], zero,
        u[1], v[1], w[1], zero,
        u[2], v[2], w[2], zero,
        zero, zero, zero, one;

      return rigid_body_map_t(m);
    }

    
    /**
     * A rigid body map that has the unit vector "u" as axis i,
     * the vector "v" as axis j, and o as origin. u and v must
     * be orthonormal.
     */
    static rigid_body_map_t basis_from(std::size_t i, const vector_t& u,
				       std::size_t j, const vector_t& v,
				       const vector_t& o = vector_t::zero()) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();

#ifndef NDEBUG
      const value_t eps = scalar_math<value_t>::epsilon();
      SL_REQUIRE("u is normalized", sl::abs(u.two_norm() - one) < 100.0f*eps);
      SL_REQUIRE("v is normalized", sl::abs(v.two_norm() - one) < 100.0f*eps);
      SL_REQUIRE("u,v are orthogonal",  sl::abs(u.dot(v)) < 100.0f*eps);
      SL_REQUIRE("i, j valid", i<3 && j<3 and i!=j);
#endif
      
      const vector_t w = u.cross(v);

      matrix_t m=tags::not_initialized();

      switch (i) {
      case 0:
        switch (j) {
        case 1:
          m =
            u[0], v[0], w[0], o[0],
            u[1], v[1], w[1], o[1],
            u[2], v[2], w[2], o[2],
            zero, zero, zero, one;
          break;
        case 2:
          m =
            u[0], -w[0], v[0], o[0],
            u[1], -w[1], v[1], o[1],
            u[2], -w[2], v[2], o[2],
            zero, zero, zero, one;
          break;
        }
        break;
      case 1:
        switch (j) {
        case 0:
         m =
            v[0], u[0], -w[0], o[0],
            v[1], u[1], -w[1], o[1],
            v[2], u[2], -w[2], o[2],
            zero, zero, zero, one;
          break;
        case 2:
          m =
            w[0], u[0], v[0], o[0],
            w[1], u[1], v[1], o[1],
            w[2], u[2], v[2], o[2],
            zero, zero, zero, one;
          break;
        }
        break;
      case 2:
        switch (j) {
        case 0:
         m =
            v[0], w[0], u[0], o[0],
            v[1], w[1], u[1], o[1],
            v[2], w[2], u[2], o[2],
            zero, zero, zero, one;
          break;
        case 1:
          m =
            -w[0], v[0], u[0], o[0],
            -w[1], v[1], u[1], o[1],
            -w[2], v[2], u[2], o[2],
            zero, zero, zero, one;
          break;
        }
        break;
      }
      return rigid_body_map_t(m);
    }
      
      
    /**
     * A rigid body map that rotates the normalized "from" vector 
     * into the normalized "to" vector. Based on 
     * Tomas Moeller and John Hughes, "Efficiently Building a
     * Matrix to Rotate One Vector to Another", ACM Journal of
     * Graphics Tools, 4(4): 1--4.
     */
    static rigid_body_map_t rotation(const vector_t& from,
                                     const vector_t& to) {
      const value_t one = scalar_math<value_t>::one();
      const value_t two = scalar_math<value_t>::two();
      const value_t zero = scalar_math<value_t>::zero();
      const value_t eps = scalar_math<value_t>::epsilon();

#ifndef NDEBUG
      SL_REQUIRE("from is normalized", sl::abs(from.two_norm() - one) < 100.0f*eps);
      SL_REQUIRE("to is normalized", sl::abs(to.two_norm() - one) < 100.0f*eps);
#endif
      
      matrix_t m=tags::not_initialized();

      vector_t v = from.cross(to);
      const value_t e = from.dot(to);

      if (sl::abs(e) > one - eps) {
        // "from" and "to"-vector almost parallel 
        
        // coordinate axis most nearly orthogonal to "from"
        const vector_t x = vector_t::unit(from.iamin());
        
        vector_t u = x - from;
        v = x - to;
        
        const value_t c1 = two / u.dot(u);
        const value_t c2 = two / v.dot(v);
        const value_t c3 = c1 * c2 * u.dot(v);
        
        m.to_identity();
        for (size_t i = 0; i < 3; i++) {
          for (size_t j = 0; j < 3; j++) {
            m(j,i) +=  
              - c1 * u[i] * u[j]
              - c2 * v[i] * v[j]
              + c3 * v[i] * u[j];
          }
        }
      } else {
        // the most common case, unless "from"="to", or "from"= -"to"
        const value_t h = (one - e)/v.dot(v);
        m(0,0) = e + h * v[0] * v[0];  
        m(0,1) = h * v[0] * v[1] - v[2]; 
        m(0,2) = h * v[0] * v[2] + v[1];
        m(0,3) = zero;

        m(1,0) = h * v[0] * v[1] + v[2]; 
        m(1,1) = e + h * v[1] * v[1];    
        m(1,2) = h * v[1] * v[2] - v[0];
        m(1,3) = zero;
	
        m(2,0) = h * v[0] * v[2] - v[1]; 
        m(2,1) = h * v[1] * v[2] + v[0]; 
        m(2,2) = e + h * v[2] * v[2];
        m(2,3) = zero;

        m(3,0) = zero;
        m(3,1) = zero;
        m(3,2) = zero;
        m(3,3) = one;
#if 0
        /* ...otherwise use this hand optimized version (9 mults less) */
        const value_t hvx=h*v[0];
        const value_t hvz=h*v[2];
        const value_t hvxy=hvx*v[1];
        const value_t hvxz=hvx*v[2];
        const value_t hvyz=hvz*v[1];
        m(0,0)=e+hvx*v[0]; m(1,0)=hvxy-v[2];     m(2,0)=hvxz+v[1];  m(3,0)=zero;
        m(0,1)=hvxy+v[2];  m(1,1)=e+h*v[1]*v[1]; m(2,1)=hvyz-v[0];  m(3,1)=zero;
        m(0,2)=hvxz-v[1];  m(1,2)=hvyz+v[0];     m(2,2)=e+hvz*v[2]; m(3,2)=zero;
        m(0,3)=zero;       m(1,3)=zero;          m(2,3)=zero;       m(3,3)=one;
        m.transpose();
#endif
      }
      
      return rigid_body_map_t(m);
    }
    
    /// Rotation specifyd by axisidx and angle
    static rigid_body_map_t rotation(size_t axisidx, 
                                     value_t angle) {
      SL_REQUIRE("Good axisidx", axisidx <= 2);
      matrix_t m=tags::not_initialized();
      m.to_identity();
      value_t c=std::cos(angle),s=std::sin(angle);
      switch (axisidx) {
      case 0:
        m(0, 0)=1; m(1, 1)=c; m(2, 1)=s; m(1, 2)= -s; m(2, 2)=c;
        break;
      case 1:
        m(1, 1)=1; m(2, 2)=c; m(0, 2)=s; m(2, 0)= -s; m(0, 0)=c;
        break;
      case 2:
        m(2, 2)=1; m(0, 0)=c; m(1, 0)=s; m(0, 1)= -s; m(1, 1)=c;
        break;
      }
      return rigid_body_map_t(m);
    } 

    /// Euler rotation: rot_z(roll)*rot_y(yaw)*rot_x(pitch)
    static rigid_body_map_t rotation(value_t pitch, value_t yaw, value_t roll) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();
      value_t cx = std::cos(pitch), sx = std::sin(pitch);
      value_t cy = std::cos(yaw),   sy = std::sin(yaw);
      value_t cz = std::cos(roll),  sz = std::sin(roll);
      
      return rigid_body_map_t(matrix_t(cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz, zero,
                                       cy*sz, sx*sy*sz+cx*cz, cx*sy*sz-sx*cz, zero,
                                       -sy,   sx*cy,          cx*cy,          zero,
                                       zero,  zero,           zero,           one));
    } 

    /// rotation from quaternion
    static rigid_body_map_t rotation(const quaternion_t& q) {
      return rigid_body_map_t(q);
    }

    /// rotation from axis, angle. Axis must be normalized.
    static rigid_body_map_t rotation(const vector_t& axis, const value_t& angle) {
      quaternion_t q; q.from_axis_angle(axis, angle);
      return rigid_body_map_t(q);
    }

    /**
     *  Look at map, for a camera positioned at vp and looking 
     *  at rp with a given up vector
     *  Undefined when up is parallel to vp-rp or vp is equal to rp.
     */
    static rigid_body_map_t lookat(const point_t& vp, 
                                   const point_t& rp, 
                                   const vector_t& up) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();

      vector_t cz  = (vp - rp).ok_normalized();
      vector_t cx  = (up.cross(cz)).ok_normalized();
      vector_t cy  = cz.cross(cx);

      rigid_body_map_t result = rigid_body_map_t(matrix_t(cx[0], cx[1], cx[2], zero,
							  cy[0], cy[1], cy[2], zero,
							  cz[0], cz[1], cz[2], zero,
							  zero,  zero,  zero,  one)); // Inv rot!
      result *= translation(-vp[0], -vp[1], -vp[2]); // Inv transl!

      return result;
    }

    /// look at map, for a camera positioned at vp and looking at rp
    static rigid_body_map_t lookat(const point_t& vp, 
                                   const point_t& rp, 
                                   value_t twist) {
      
      rigid_body_map_t result;
      
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();

      value_t pvx,pvy,pvz;
      value_t sin_roty,cos_roty,sin_rotx,cos_rotx,rxz,rxyz;
      
      pvx = sqr(rp[0]-vp[0]);
      pvy = sqr(rp[1]-vp[1]);
      pvz = sqr(rp[2]-vp[2]);
      rxz = std::sqrt(pvx+pvz);
      rxyz = std::sqrt(pvx+pvy+pvz);
    
      result = rotation(2,-twist);
      
      if (rxyz >= std::numeric_limits<value_t>::epsilon()) {
	
        sin_rotx = (vp[1]-rp[1])/ rxyz;
        cos_rotx = rxz/ rxyz;
	
        result *= rigid_body_map_t(matrix_t(one, zero,       zero,        zero,
                                            zero, cos_rotx, -sin_rotx, zero,
                                            zero, sin_rotx, cos_rotx,  zero,
                                            zero, zero,       zero,        one));
      }
      
      if (rxz >= std::numeric_limits<value_t>::epsilon()) {
        sin_roty = (rp[0]-vp[0])/ rxz;
        cos_roty = (vp[2]-rp[2])/ rxz;
	
        result *= rigid_body_map_t(matrix_t(cos_roty,  zero, sin_roty, zero,
                                            zero,        one, zero,       zero,
                                            -sin_roty, zero, cos_roty, zero,
                                            zero,        zero, zero,       one));
      }

      result *= translation(-vp[0], -vp[1], -vp[2]);
      
      return result;
    } 
    /**
     * An orthographic projection transformation 
     * @param l   expects the coordinate for the left vertical clipping plane.
     * @param r   expects the coordinate for the right vertical clipping plane.
     * @param b   expects the coordinate for the bottom horizontal clipping plane.
     * @param t   expects the coordinate for the top horizontal clipping plane.
     * @param n   expects the distance to the nearer depth clipping plane.
     * @param f   expects the distance to the farther depth clipping plane. 
     */
    static projective_map_t ortho(value_t l, value_t r,
                                  value_t b, value_t t,
                                  value_t n, value_t f) {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();

      return 
        projective_map_t
        (matrix_t
         (two/(r-l), zero,       zero,      -(r+l)/(r-l),
          zero,      two/(t-b),  zero,      -(t+b)/(t-b),
          zero,      zero,      -two/(f-n), -(f+n)/(f-n),
          zero,      zero,       zero,              one));
    }                                    
       
    /**
     * A perspective projection transformation
     * @param fov    expects the field-of-view angle in the y direction. 
     * @param aspect expects the aspect ratio which determines the field of 
     *               view in the x direction.  The aspect ratio is the 
     *               ratio of x (width) to y (height).
     * @param v_near expects the distance from the viewer to the closest clipping
     *               plane (always positive).
     * @param v_far  expects the distance from the viewer to the farthest clipping
     *               plane (always positive).
     */                                   
    static projective_map_t perspective(value_t fov, 
                                        value_t aspect, 
                                        value_t v_near, 
                                        value_t v_far) {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();
      const value_t cot = one/std::tan(fov/two);
      return 
        projective_map_t
        (matrix_t
         (cot/aspect,     zero,    zero,                                            zero,  
          zero,            cot,    zero,                                            zero,
          zero,            zero,    -(v_far + v_near)/(v_far - v_near),           -(two * v_far * v_near)/(v_far - v_near),
          zero,            zero,    -one,                                            zero));
    }

    /// A perspective projection transformation with an axpect ratio of one
    static projective_map_t perspective(value_t fov, 
                                        value_t v_near, 
                                        value_t v_far) {
      const value_t one  = scalar_math<value_t>::one();
      return perspective(fov,one,v_near,v_far);
    }

    /**
     * A perspective projection transformation
     * @param l   expects the coordinate for the left vertical clipping plane.
     * @param r   expects the coordinate for the right vertical clipping plane.
     * @param b   expects the coordinate for the bottom horizontal clipping plane.
     * @param t   expects the coordinate for the top horizontal clipping plane.
     * @param n   expects the distance to the nearer depth clipping plane.
     * @param f   expects the distance to the farther depth clipping plane. 
     */
    static projective_map_t frustum(value_t l, 
                                    value_t r, 
                                    value_t b, 
                                    value_t t, 
                                    value_t n, 
                                    value_t f) {
      const value_t zero = scalar_math<value_t>::zero();
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();
      return 
        projective_map_t
        (matrix_t
         (two*n/(r-l),    zero,             (r+l)/(r-l), zero,
          zero,                two*n/(t-b), (t+b)/(t-b), zero,
          zero,                   zero,            -(f+n)/(f-n), -two*n*f/(f-n),
          zero,                   zero,             -one,         zero));
    }

    /**
     * A stereoscopic perspective projection transformation
     * @param fovy   expects the field-of-view angle in the y direction. 
     * @param ratio  expects the aspect ratio which determines the field of 
     *               view in the x direction.  The aspect ratio is the 
     *               ratio of x (width) to y (height).
     * @param v_near expects the distance from the viewer to the closest clipping
     *               plane (always positive).
     * @param v_far  expects the distance from the viewer to the farthest clipping
     *               plane (always positive).
     * @param dist   expects the distance to the zero-parallax plane
     * @param eye    expects the inter-pupillary distance (positive for right, negative for left) 
     */                                   
    static projective_map_t stereo_perspective(value_t fovy, 
                                               value_t ratio, 
                                               value_t v_near, 
                                               value_t v_far,
                                               value_t dist, 
                                               value_t eye) {
      const value_t zero = scalar_math<value_t>::zero();
      SL_REQUIRE("Good aspect ratio", ratio > zero);
      const value_t two  = scalar_math<value_t>::two();
      
      const value_t gltan  = std::tan(fovy/two);
      const value_t gltan2 = std::tan(fovy*ratio/two);
      // the projection planes.
      const value_t e = eye / dist * v_near;

      return 
        frustum(-gltan2*v_near-e, 
                gltan2*v_near-e, 
                -gltan*v_near, 
                gltan*v_near, 
                v_near, 
                v_far) * 
        translation(-eye,zero,zero);
    } 
     
    /**
     * A default stereoscopic perspective projection transformation. 
     * The zero-parallax plane to near (eveything will appear behind the screen),
     * and the inter-pupillary distance is automatically computed to get a "good" 
     * stereo effect.
     * @param fovy   expects the field-of-view angle in the y direction. 
     * @param ratio  expects the aspect ratio which determines the field of 
     *               view in the x direction.  The aspect ratio is the 
     *               ratio of x (width) to y (height).
     * @param v_near expects the distance from the viewer to the closest clipping
     *               plane (always positive).
     * @param v_far  expects the distance from the viewer to the farthest clipping
     *               plane (always positive).
     * @param right  true iff this is the right eye.
     */                                   
    static projective_map_t stereo_perspective(value_t v_fovy, 
                                               value_t ratio, 
                                               value_t v_near, 
                                               value_t v_far,
                                               bool right) {
#ifndef NDEBUG
      const value_t zero = scalar_math<value_t>::zero();
      SL_REQUIRE("Good aspect ratio", ratio > zero);
#endif
      const value_t one  = scalar_math<value_t>::one();
      const value_t two  = scalar_math<value_t>::two();

      value_t h_fovy = v_fovy * ratio;

      value_t zmin = -v_far;
      value_t zmax = -v_near;
 
      // We do not allow negative parallax. Everything will appear
      // behind the screen
      
      value_t dx = two * v_near * std::tan(h_fovy/two);
      value_t d  = v_near;
 
      value_t parallax_pos = std::tan(static_cast<value_t>(0.6) * scalar_math<value_t>::Pi() / static_cast<value_t>(180.0) / two)/std::tan(h_fovy/two)*dx;
 
      value_t tc = parallax_pos * (one + d / (zmax - zmin));

      return stereo_perspective(v_fovy, ratio, v_near, v_far, -zmax, tc/(right?two:-two));
    }             

     
    /**
     * Transform the standard centered perspective map projection into a "good"
     * stereoscopic map. The position of the clipping planes and the
     * viewing angle is extracted from the projection, 
     * the zero-parallax plane is set to near (eveything will appear 
     * behind the screen), and the inter-pupillary distance is 
     * automatically computed to get a "good" stereo effect.
     * @param projection expects a standard 3D perspective
     * @param ratio  expects the aspect ratio which determines the field of 
     *               view in the x direction.  The aspect ratio is the 
     *               ratio of x (width) to y (height).
     * @param right  true iff this is the right eye.
     */                                   
    static projective_map_t stereo_perspective(const projective_map_t& projection,
                                               value_t ratio, 
                                               bool right) {
      SL_REQUIRE("Is std perspective", projection.is_std_centered_3d_perspective());

      value_t v_near = projection.near_from_std_3d_perspective();
      value_t v_far  = projection.far_from_std_3d_perspective();
      value_t v_fovy = projection.fov_from_std_3d_perspective();

      return stereo_perspective(v_fovy, ratio, v_near, v_far,right);
    }
    
  }; // class linear_map_factory<3,T>

}; // namespace sl

// ---------------------------------------------------------------------
// sl::linear_map_factory<2,T>
// ---------------------------------------------------------------------

namespace sl {

  /// Objects that build 2D affine maps
  template <class T>
  class linear_map_factory<2,T>: 
    public linear_map_factory_base<2,T> {
  public:
    typedef linear_map_factory_base<2,T> super_t;

    typedef typename super_t::value_t           value_t;
    typedef typename super_t::projective_map_t  projective_map_t;
    typedef typename super_t::rigid_body_map_t  rigid_body_map_t;
    typedef typename super_t::affine_map_t      affine_map_t;
    typedef typename super_t::vector_t          vector_t;
    typedef typename super_t::dual_vector_t     dual_vector_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::plane_t           plane_t;
    typedef typename super_t::matrix_t          matrix_t;

  public:

    /// Identity
    static rigid_body_map_t identity() {
      return super_t::identity();
    }

    /// Translation
    static rigid_body_map_t translation(const vector_t& v) {
      return super_t::translation(v);
    }

    /// translation
    static rigid_body_map_t translation(value_t x,
                                        value_t y) {
      return translation(vector_t(x,y));
    }

    /// Scaling
    static affine_map_t scaling(const vector_t& vv) {
      return super_t::scaling(vv);
    }

    /// Uniform scaling
    static affine_map_t scaling(value_t s) {
      return super_t::scaling(s);
    }

    /// Scaling
    static affine_map_t scaling(value_t x,
                                value_t y) {
      return scaling(vector_t(x,y));
    }

    /// Shear
    static affine_map_t shear(value_t xy) { 
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();
      return affine_map_t(matrix_t( one,  zero, zero,
                                    xy,   one, zero,
                                    zero,  zero,  one));
    } 

    /// Rotation
    static rigid_body_map_t rotation(value_t angle) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();
      const value_t c=std::cos(angle),s=std::sin(angle);
      return rigid_body_map_t(matrix_t(   c,   -s, zero,
                                          s,    c, zero,
                                          zero, zero, one));
    }

    /**
     * A rigid body map that has the unit vector "u" as first axis.
     * John Hughes and Tomas Moeller, "Efficiently Building a
     * Matrix to Rotate One Vector to Another", ACM Journal of
     * Graphics Tools, 4(4): 33--35.
     */
    static rigid_body_map_t basis_from(const vector_t& u) {
      const value_t one = scalar_math<value_t>::one();
      const value_t zero = scalar_math<value_t>::zero();

#ifndef NDEBUG
      const value_t eps = scalar_math<value_t>::epsilon();
      SL_REQUIRE("u is normalized", sl::abs(u.two_norm() - one) < 100.0f*eps);
#endif
      
      matrix_t m=tags::not_initialized();
      m = 
         u[0], -u[1], zero,
         u[1],  u[0], zero,
         zero,  zero,  one;

      return rigid_body_map_t(m);
    } 
    
  }; // class linear_map_factory<2,T>

}; // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  /// A 2D linear map factory with single precision floating point components.
  typedef linear_map_factory<2,float> linear_map_factory2f;
  /// A 3D linear map factory with single precision floating point components.
  typedef linear_map_factory<3,float> linear_map_factory3f;
  /// A 4D linear map factory with single precision floating point components.
  typedef linear_map_factory<4,float> linear_map_factory4f;

  /// A 2D linear map factory with double precision floating point components.
  typedef linear_map_factory<2,double> linear_map_factory2d;
  /// A 3D linear map factory with double precision floating point components.
  typedef linear_map_factory<3,double> linear_map_factory3d;
  /// A 4D linear map factory with double precision floating point components.
  typedef linear_map_factory<4,double> linear_map_factory4d;
}

#endif
