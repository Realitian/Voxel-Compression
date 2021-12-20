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
#ifndef SL_IDENTITY_MAP_HPP
#define SL_IDENTITY_MAP_HPP

#include <sl/conv_to.hpp>
#include <sl/linear_map.hpp>
#include <sl/quaternion.hpp>

// --------------------------------------------------------------------
// sl::identity_map<DIMENSION,T>
// --------------------------------------------------------------------

namespace sl {

  /**
   *  N-dimensional identity maps
   */
  template <size_t N_dimension, class T>
  class identity_map: 
    public linear_map_base<sl::identity_tag, N_dimension, T>
  {
  public: // Constants and types
    typedef linear_map_base<sl::identity_tag, N_dimension, T> super_t;

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
    inline identity_map() {
      SL_INVARIANT(this->invariant());
    }

    /// Fast init (garbage), handle with care!
    inline identity_map(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Init from matrix m, assuming m is identity
    inline identity_map(const matrix_t& m): super_t(m) {
      SL_REQUIRE("Compatible matrix", is_compatible(m));
      SL_INVARIANT(this->invariant());
    }

    /// Convert from compatible linear map
    inline identity_map(const linear_map_base<
			sl::identity_tag,
			dimension, 
			value_t>& other): super_t(other) {
      SL_INVARIANT(this->invariant());
    }

  };

  template <std::size_t DIM, typename OUT_ET>
  class conv_to< identity_map<DIM, OUT_ET> > {
  public:
    typedef identity_map<DIM, OUT_ET> result_t;
    
    // Explicit conversion from maps of another type
    template <typename IN_ET>
    inline static result_t from(const identity_map<DIM,IN_ET>& /*in*/) {
      return result_t();
    }
    
  }; // class conv_to

} // namespace sl


// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  /// A 2D identity map with single precision floating point components.
  typedef identity_map<2,float> identity_map2f;
  /// A 3D identity map with single precision floating point components.
  typedef identity_map<3,float> identity_map3f;
  /// A 4D identity map with single precision floating point components.
  typedef identity_map<4,float> identity_map4f;

  /// A 2D identity map with double precision floating point components.
  typedef identity_map<2,double> identity_map2d;
  /// A 3D identity map with double precision floating point components.
  typedef identity_map<3,double> identity_map3d;
  /// A 4D identity map with double precision floating point components.
  typedef identity_map<4,double> identity_map4d;
}

#endif




