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
#ifndef SL_FIXED_SIZE_TETRAHEDRON_HPP
#define SL_FIXED_SIZE_TETRAHEDRON_HPP

#include <sl/fixed_size_plane.hpp>
#include <cassert>

namespace sl {

  /// Base class for tetrahedra of fixed dimension
  template <class SELF_T, size_t DIMENSION, class T>
  class fixed_size_tetrahedron_base
  {
  public: // Constants and types

    enum { dimension = DIMENSION };

    typedef SELF_T  this_t;
    typedef T       value_t;

    typedef fixed_size_point<dimension, T>                        point_t;
    typedef fixed_size_plane<dimension, T>                        plane_t;
    typedef fixed_size_vector<column_orientation, dimension, T>   vector_t;
    typedef fixed_size_vector<row_orientation, dimension, T>      dual_vector_t;

    typedef fixed_size_array<4, point_t> storage_t;
    
  protected: // Storage

    storage_t corner_; // Corners

  public:

    inline fixed_size_tetrahedron_base() {
    }
    
    inline fixed_size_tetrahedron_base(const point_t& p0,
                                       const point_t& p1,
                                       const point_t& p2,
                                       const point_t& p3) {
      corner_[0] = p0;
      corner_[1] = p1;
      corner_[2] = p2;
      corner_[3] = p3;
    }

    inline const point_t& corner(std::size_t i) const {
      assert(i<4);
      return corner_[i];
    }

    inline point_t& corner(std::size_t i) {
      assert(i<4);
      return corner_[i];
    }

    inline int longest_edge() const {
      std::size_t ilongest = 0;
      value_t dlongest = corner(0).distance_squared_to(corner(1));
      
      value_t dcurrent = corner(0).distance_squared_to(corner(2));
      if (dcurrent>dlongest) { ilongest=1; dlongest = dcurrent; }
      dcurrent = corner(0).distance_squared_to(corner(3));
      if (dcurrent>dlongest) { ilongest=2; dlongest = dcurrent; }
      dcurrent = corner(1).distance_squared_to(corner(2));
      if (dcurrent>dlongest) { ilongest=3; dlongest = dcurrent; }
      dcurrent = corner(1).distance_squared_to(corner(3));
      if (dcurrent>dlongest) { ilongest=4; dlongest = dcurrent; }
      dcurrent = corner(2).distance_squared_to(corner(3));
      if (dcurrent>dlongest) { ilongest=5; dlongest = dcurrent; }

      return ilongest;
    }
    
    /// The center of the tetrahedron (unique identifier at a given level)
    inline point_t center() const {
      return
        corner(0).lerp(corner(1).lerp(corner(2).lerp(corner(3), 0.5f), 0.5f), 0.5f);
    }

    inline int orientation() const {
      return orientation(corner(0), corner(1), corner(2), corner(3));
    }
    
    inline bool contains(const point_t& p) const {
      // FIXME: works only for 3D
      const int isign_abcd = 1;
      assert(isign_abcd == orientation(corner(0), corner(1), corner(2), corner(3)));

      vector_t ap =         p - corner(0);
      vector_t ac = corner(2) - corner(0);
      vector_t ad = corner(3) - corner(0);
      int result = (-isign_abcd != det3x3_sign(ap, ac, ad));

      if (result) {
        vector_t ab = corner(1) - corner(0);
        result = (-isign_abcd != det3x3_sign(ab, ap, ad));

        if (result) {
          result = (-isign_abcd != det3x3_sign(ab, ac, ap));
          
          if (result) {
            vector_t pb = corner(1) - p;
            vector_t pc = corner(2) - p;
            vector_t pd = corner(3) - p;

            result = (-isign_abcd != det3x3_sign(pb, pc, pd));
          }
        }
      }
      
      return result;
    }
    
    /// The tetrahedron's bounding planes 
    void planes_in(plane_t planes[4]) const {
      point_t p0 = corner(0);
      point_t p1 = corner(1);
      point_t p2 = corner(2);
      point_t p3 = corner(3);
      planes[0] = stable_plane(p2,p1,p0); 
      planes[1] = stable_plane(p3,p1,p2); 
      planes[2] = stable_plane(p1,p3,p0); 
      planes[3] = stable_plane(p3,p2,p0); 
    }

  public: // Longest edge refinement - assumes edge 0 is the longest

    // True iff edge 0 is the longest
    inline bool is_canonical() const {
      return longest_edge() == 0;
    }

    /// The midpoint of longest edge (unique diamond identifier)
    inline point_t canonical_split_vertex() const {
      assert(is_canonical());
      return corner(0).lerp(corner(1), 0.5f);
    }

    /// The tetrahedron's split plane. Positive side is child 1
    plane_t canonical_split_plane() const {
      assert(is_canonical());
      point_t p0 = corner(0);
      point_t p2 = corner(2);
      point_t p3 = corner(3);
      point_t p01 = canonical_split_vertex();
      plane_t result = stable_plane(p01,p2,p3);
      if (result.value(p0)>0.0f) {
        result = stable_plane(p01,p3,p2); 
      }
      return result;
    }
    
    /// The children cells obtained by longest edge refinement
    inline this_t canonical_child(std::size_t i) const {
      assert(is_canonical());
      assert(i<2);
      point_t pnew = canonical_split_vertex();

      if (i==0) {
        if (corner(0).distance_squared_to(corner(3)) > corner(0).distance_squared_to(corner(2))) {
          return this_t(corner(0), corner(3), pnew, corner(2));
        } else {
          return this_t(corner(0), corner(2), corner(3), pnew);
        }
      } else {
        if (corner(1).distance_squared_to(corner(2)) > corner(1).distance_squared_to(corner(3))) {
          return this_t(corner(1), corner(2), pnew, corner(3));
        } else {
          return this_t(corner(1), corner(3), corner(2), pnew);
        }
      }
    }
    
    /// The children cells obtained by longest edge refinement
    inline std::size_t canonical_child_index(const point_t& p) const {
      assert(is_canonical());
      point_t pnew = canonical_split_vertex();
      int isign_abcd=orientation(pnew, p, corner(2), corner(3));
      return (isign_abcd <= 0) ? 0 : 1;
    }

    /// The children cells obtained by longest edge refinement
    inline this_t canonical_child(const point_t& p) const {
      return canonical_child(canonical_child_index(p));
    }

  public:
    
    inline bool operator==(const this_t& other) const {
      return
        corner(0) == other.corner(0) &&
        corner(1) == other.corner(1) &&
        corner(2) == other.corner(2) &&
        corner(3) == other.corner(3);
    }

    inline bool operator!=(const this_t& other) const {
      return !((*this) == other);
    }
    
  protected: // Helpers

    static inline plane_t stable_plane(const point_t& p0,
                                       const point_t& p1,
                                       const point_t& p2) {
      if (p0<p1) {
        if (p0<p2) {
          //012
          return plane_t(sl::normal(p0,p1,p2),p0);
        } else {
          // 201
          return plane_t(sl::normal(p2,p0,p1), p2);
        }
      } else {
        if (p1<p2) {
          // 120
          return plane_t(sl::normal(p0,p1,p2),p0);
        } else {
          // 201
          return plane_t(sl::normal(p2,p0,p1), p2);
        }
      }
    }
    
    static inline int det3x3_sign(value_t ax, value_t ay, value_t az,
                                  value_t bx, value_t by, value_t bz,
                                  value_t cx, value_t cy, value_t cz) {
      value_t det = 0;
      det += ax * (by * cz  -  bz * cy);
      det += ay * (bz * cx  -  bx * cz);
      det += az * (bx * cy  -  by * cx);
      int result = ((det>0) ? 1 : ((det<0) ? -1 : 0));
      
      return result;
    }
    
    static inline int det3x3_sign(const vector_t& ab,
                                   const vector_t& ac,
                                   const vector_t& ad) {
      return det3x3_sign(ab[0], ab[1], ab[2],
                         ac[0], ac[1], ac[2],
                         ad[0], ad[1], ad[2]);
    }

    static inline int orientation(const point_t& a,
                                  const point_t& b,
                                  const point_t& c,
                                  const point_t& d) {
      // FIXME: works only for 3D
      return det3x3_sign(b-a, c-a, d-a);
    }
  };

} // namespace sl

namespace sl {

  /// Tetrahedra of fixed dimension
  template <size_t DIMENSION, class T>
  class fixed_size_tetrahedron: 
    public fixed_size_tetrahedron_base< fixed_size_tetrahedron<DIMENSION, T>, DIMENSION, T > {

  public: // Types

    typedef fixed_size_tetrahedron_base< fixed_size_tetrahedron<DIMENSION, T>, DIMENSION, T > super_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::this_t        this_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::point_t       point_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;
    typedef typename super_t::plane_t       plane_t;

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_tetrahedron() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_tetrahedron(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    inline fixed_size_tetrahedron(const point_t& p0,
                                  const point_t& p1,
                                  const point_t& p2,
                                  const point_t& p3)
        : super_t(p0, p1, p2, p3) {
    }

  };

}; // namespace sl

template <class SELF_T, size_t DIMENSION, class T>
std::ostream& operator <<(std::ostream& s, const sl::fixed_size_tetrahedron_base<SELF_T,DIMENSION,T>& a) {
  for (size_t i = 0; i<4; i++) {
      s << a.corner(i) << std::endl;;
  }
  return s;
}
    
template <class SELF_T, size_t DIMENSION, class T>
std::istream& operator >>(std::istream& s, sl::fixed_size_tetrahedron_base<SELF_T,DIMENSION,T>& a) {
  for (size_t i = 0; i<4; i++) {
      s >> a.corner(i) ;
  }
  return s;
}


namespace sl {

  /// 2D tetrahedron with single precision floating tetrahedron components
  typedef fixed_size_tetrahedron<2,float> tetrahedron2f;
  /// 3D tetrahedron with single precision floating tetrahedron components
  typedef fixed_size_tetrahedron<3,float> tetrahedron3f;
  /// 4D tetrahedron with single precision floating tetrahedron components
  typedef fixed_size_tetrahedron<4,float> tetrahedron4f;

  /// 2D tetrahedron with double precision floating tetrahedron components
  typedef fixed_size_tetrahedron<2,double> tetrahedron2d;
  /// 3D tetrahedron with double precision floating tetrahedron components
  typedef fixed_size_tetrahedron<3,double> tetrahedron3d;
  /// 4D tetrahedron with double precision floating tetrahedron components
  typedef fixed_size_tetrahedron<4,double> tetrahedron4d;
}


#endif
