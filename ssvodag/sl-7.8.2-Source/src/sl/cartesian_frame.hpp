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
#ifndef SL_CARTESIAN_FRAME_HPP
#define SL_CARTESIAN_FRAME_HPP

#include <sl/affine_map.hpp>

// --------------------------------------------------------------------
// -- cartesian_frame_base<SELF_T,DIMENSION,T>
// --------------------------------------------------------------------

namespace sl {

  /// Base class for cartesian frames
  template <class SELF_T, size_t DIMENSION, class T>
  class cartesian_frame_base
  {
  public: // Constants and types

    enum { dimension = DIMENSION };

    typedef SELF_T  self_t;
    typedef T       value_t;

    typedef affine_map<dimension,T>                        affine_map_t;
    typedef typename affine_map_t::vector_t                vector_t;
    typedef typename affine_map_t::dual_vector_t           dual_vector_t;
    typedef typename affine_map_t::point_t                 point_t;
    typedef typename affine_map_t::matrix_t                matrix_t;
    
  protected: // Storage

    affine_map_t to_parent_frame_;
    affine_map_t from_parent_frame_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << to_parent_frame_ << from_parent_frame_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> to_parent_frame_ >> from_parent_frame_;
    }

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", DIMENSION > 0);
    
  public: // Invariant

    bool invariant() const {
      return true; // CHECK: from_parent_frame_ * to_parent_frame_ == identity
    }

  protected: // Creation, Copy & Destruction

    /// Default init (identity maps)
    inline cartesian_frame_base() {
      SL_INVARIANT(invariant());
    }

    /// Fast init (garbage), handle with care!
    inline cartesian_frame_base(tags::not_initialized tag): from_parent_frame_(tag), to_parent_frame_(tag) {
      // Garbage in, use with care
    }

    /// Explicint init from  parameters
    inline cartesian_frame_base(const affine_map_t& to_parent_frame,
				const affine_map_t& from_parent_frame): 
      to_parent_frame_(to_parent_frame), from_parent_frame_(from_parent_frame) {
      SL_INVARIANT(invariant());
    }

    /// Explicint init from parameter
    inline cartesian_frame_base(const affine_map_t& to_parent_frame):
      to_parent_frame_(to_parent_frame), from_parent_frame_(~to_parent_frame) {
      SL_INVARIANT(invariant());
    }

  public: // Element Access

    /// the map from local to global coordinates
    inline const affine_map_t& to_parent_frame_map() const {
      return to_parent_frame_;
    }

    /// the map from global to local coordinates
    inline const affine_map_t& from_parent_frame_map() const {
      return from_parent_frame_;
    }

  public: // Setting

    /// Set local to global map to other
    void set_to_parent_frame_map(const affine_map_t& other) {
      to_parent_frame_ = other;
      from_parent_frame_ = ~other;
      SL_INVARIANT(invariant());
    }

    /// Set global to local map to other
    void set_from_parent_frame_map(const affine_map_t& other) {
      from_parent_frame_ = other;
      to_parent_frame_ = ~other;
      SL_INVARIANT(invariant());
    }

    /// Set both maps 
    void set_maps(const affine_map_t& to_parent_frame,
		  const affine_map_t& from_parent_frame) {
      to_parent_frame_ = to_parent_frame;
      from_parent_frame_ = from_parent_frame;
      SL_INVARIANT(invariant());
    }
      
  public: // Multiplicative group operations

    /// Compose this frame with a frame in local coordinates
    self_t& operator *= (const self_t& local_frame) const {
      to_parent_frame_   = to_parent_frame_ * local_frame.to_parent_frame_map();
      from_parent_frame_ = local_frame.from_parent_frame_map() * from_parent_frame_;
      SL_INVARIANT(invariant());
      return *this;
    }

    /// The reference frame constructed from the inverse of this maps
    self_t inverse() const {
      self_t result(from_parent_frame_map(),to_parent_frame_map());
      return result;
    }

    /// The reference frame constructed from the inverse of this maps
    inline self_t operator~() const {
      return inverse();
    }                                                                                                                    

    /// Compose this frame with the inverse of a frame in local coordinates
    self_t& operator /= (const self_t& local_frame) const {
      *this *= inverse();
      return *this;
    }

    SL_OP_MULTIPLICATIVE_GROUP(self_t);

  }; 

}; // namespace sl


// --------------------------------------------------------------------
// -- cartesian_frame<DIMENSION,T>
// --------------------------------------------------------------------

namespace sl {

  /// An N-dimensional Cartesian frame defined by an affine map and its inverse
  template <size_t DIMENSION, class T>
  class cartesian_frame: 
    public cartesian_frame_base< cartesian_frame<DIMENSION,T>,
                                 DIMENSION,
                                 T> {
  public: // Constants & types
    
    typedef cartesian_frame_base< cartesian_frame<DIMENSION,T>,
                                 DIMENSION,
                                 T> super_t;

    typedef typename super_t::self_t              self_t;
    typedef typename super_t::value_t             value_t;
    typedef typename super_t::affine_map_t        affine_map_t;
    typedef typename super_t::vector_t            vector_t;
    typedef typename super_t::dual_vector_t       dual_vector_t;
    typedef typename super_t::point_t             point_t;
    typedef typename super_t::matrix_t            matrix_t;

  public: // Creation, copy & destruction

    /// Default init (identity)
    inline cartesian_frame() {
      SL_INVARIANT(this->invariant());
    }

    /// Fast init (garbage), handle with care!
    inline cartesian_frame(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }
    
    /// Init from params
    inline cartesian_frame( const affine_map_t& to_parent_frame,
			    const affine_map_t& from_parent_frame): 
      super_t(to_parent_frame,from_parent_frame) {
    }

    /// Init from paramss
    inline cartesian_frame( const affine_map_t& to_parent_frame):
      super_t(to_parent_frame) {
    }

  };

  
  template <std::size_t DIM, typename OUT_ET>
  class conv_to< cartesian_frame<DIM, OUT_ET> > {
  public:
    typedef cartesian_frame<DIM, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const cartesian_frame<DIM, IN_ET>& in) {
      return result_t(conv_to< typename cartesian_frame<DIM, IN_ET>::affine_map_t > (in.to_parent_frame_map()),
		      conv_to< typename cartesian_frame<DIM, IN_ET>::affine_map_t > (in.from_parent_frame_map()));
    }
  };

}; // namespace sl

/// the point p in the parent coordinates of f
template <class SELF_T, size_t DIMENSION, class T>
inline static sl::fixed_size_point<DIMENSION,T> operator * (const sl::cartesian_frame_base<SELF_T, DIMENSION, T>& f, 
							    const sl::fixed_size_point<DIMENSION,T>& p) {
  return 
    f.to_parent_frame_map() * p;
}

/// the vector v in the parent coordinates of f
template <class SELF_T, size_t DIMENSION, class T>
inline sl::fixed_size_vector<sl::column_orientation,DIMENSION,T> operator * (const sl::cartesian_frame_base<SELF_T, DIMENSION, T>& f, 
									     const sl::fixed_size_vector<sl::column_orientation,DIMENSION,T>& v) {
  return 
    f.to_parent_frame_map() * v;
}

/// the dual vector v in the parent coordinates of f
template <class SELF_T, size_t DIMENSION, class T>
inline sl::fixed_size_vector<sl::row_orientation, DIMENSION,T> operator * (const sl::cartesian_frame_base<SELF_T, DIMENSION, T>& f,
									   const sl::fixed_size_vector<sl::row_orientation,DIMENSION,T>& dual_v) {
  return 
    dual_v * f.from_parent_frame_map().as_matrix();
}

/// the (hyper)plane h in the parent coordinates of f
template <class SELF_T, size_t DIMENSION, class T>
inline sl::fixed_size_plane<DIMENSION,T> operator * (const sl::cartesian_frame_base<SELF_T, DIMENSION, T>& f,
						     const sl::fixed_size_plane<DIMENSION,T>& hp) {
  return 
    sl::fixed_size_plane<DIMENSION,T>( hp.storage() * f.from_parent_frame_map().as_matrix() );
}

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  // A 2D cartesian frame with single precision floating point components
  typedef cartesian_frame<2,float> frame2f;
  // A 3D cartesian frame with single precision floating point components
  typedef cartesian_frame<3,float> frame3f;
  // A 4D cartesian frame with single precision floating point components
  typedef cartesian_frame<4,float> frame4f;

  // A 2D cartesian frame with double precision floating point components
  typedef cartesian_frame<2,double> frame2d;
  // A 3D cartesian frame with double precision floating point components
  typedef cartesian_frame<3,double> frame3d;
  // A 4D cartesian frame with double precision floating point components
  typedef cartesian_frame<4,double> frame4d;
};

#endif
