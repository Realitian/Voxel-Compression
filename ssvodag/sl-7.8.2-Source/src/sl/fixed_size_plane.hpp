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
#ifndef SL_FIXED_SIZE_PLANE_HPP
#define SL_FIXED_SIZE_PLANE_HPP

#include <sl/fixed_size_point.hpp>

namespace sl {

  /// Base class for (hyper)planes of fixed dimension
  template <class SELF_T, size_t DIMENSION, class T>
  class fixed_size_plane_base
  {
  public: // Constants and types

    enum { dimension = DIMENSION, hdimension = DIMENSION+1 };

    typedef SELF_T  self_t;
    typedef T       value_t;

    typedef fixed_size_point<dimension, T>                        point_t;
    typedef fixed_size_vector<column_orientation, dimension, T>   vector_t;
    typedef fixed_size_vector<row_orientation, dimension, T>      dual_vector_t;

    typedef fixed_size_point<hdimension, T>                       hpoint_t;
    typedef fixed_size_vector<column_orientation, hdimension, T>  hvector_t;
    typedef fixed_size_vector<row_orientation, hdimension, T>     hdual_vector_t;

    typedef hdual_vector_t storage_t;

  protected: // Storage

    hdual_vector_t storage_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << storage_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> storage_;
    }

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", dimension > 0);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);
   
  protected: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_plane_base() {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_plane_base(tags::not_initialized tag): storage_(tag) {
      // Garbage in, use with care
    }

    /// Init from homogeneous vector
    inline explicit fixed_size_plane_base(const hdual_vector_t& hv): storage_(hv) {
    }

  public: // STL-style iterators

    typedef typename storage_t::iterator          iterator;
    typedef typename storage_t::restrict_iterator restrict_iterator;
    typedef typename storage_t::const_iterator    const_iterator;

    /// The number of components
    static inline size_t size()                { return hdimension; }
    
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

    /// The homogeneous vector representing this
    inline const hdual_vector_t& as_vector() const {
      return storage_;
    } 

    /// The storage area
    inline const hdual_vector_t& storage() const {
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

    /// The "zero" origin
    static self_t zero() {
      self_t result = tags::not_initialized();
      result.to_zero();
      return result;
    }

    /// The plane passing throuth the origin and whose normal is the i-th coordinate axis
    static self_t perpendicular_to(size_t i) {
      SL_REQUIRE("Good index", i<dimension);
      self_t result = zero();
      result[i] = sl::one(value_t());
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
    bool is_epsilon_equal(const self_t& t2,
			  value_t eps) const {
      SL_USEVAR(eps);
      return storage_.is_epsilon_equal(t2.storage_);
    }

  public: // Algebra

    // Nothing here for now

  public: // Interpolation

    /// Linear interpolation
    template <class T_PARAMETER>
    self_t lerp(const self_t& other, T_PARAMETER t) const {
      return self_t( as_vector().lerp(other.as_vector(), t) );
    }

  public: // Plane queries and operations

    /// The normal to the plane (not unit)
    inline dual_vector_t normal() const {
      dual_vector_t result;
      for (size_t i=0; i<dimension; ++i) {
	result[i] = storage_[i];
      }
      return result;
    }
  
    /// The offset from origin
    inline value_t offset() const {
      return -storage_[dimension];
    }

    /// Make the normal of unit length
    inline void normalize() {
      SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) nlen = normal().two_norm();
      SL_REQUIRE("Non null", sl::is_positive(nlen));
      for (size_t i=0; i<hdimension; ++i) {
	storage_[i] /= nlen;
      }
    }

    /// An arbitrary point on the plane
    inline point_t point_on_plane() const {
      point_t result; // zero
      size_t i = normal().iamax();
      result[i] = -storage_[dimension] / normal()[i]; 
      return result;
    }

    /// Make the normal point to the opposite direction
    inline void swap_orientation() {
      for (size_t i=0; i<hdimension; ++i) {
	storage_[i] = -storage_[i];
      }
    }

    /// The evaluation of the plane's affine functional for point p
    inline SL_SUMTYPENAME(value_t) value(const point_t& p) const {
      SL_SUMTYPENAME(value_t) result = this->storage_[dimension];
      for (size_t i=0; i<dimension; ++i) {
	result += this->storage_[i] * p[i];
      }
      return result;
    }                     
 
    /// Signed Euclidean distance from this to p, assumed normalized
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) signed_distance(const point_t& p) const {
      return value(p);
    }

    /// Euclidean distance from this to p, assumed normalized
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) distance(const point_t& p) const {
      return sl::abs(signed_distance(p));
    }
      
    /**
     * ray - plane intersection. The routine returns the 
     * location along the ray origin->extremity of the intersection with
     * this which is closer to the origin. -1 is returned if no 
     * intersection exists.
     */
    value_t intersection(const point_t& origin,
			 const point_t& extremity) const {
      value_t result = -sl::scalar_math<value_t>::one();
      value_t denom = as_dual(extremity - origin).dot(normal());
      if (!sl::is_zero(denom)) {
	result = -value(origin) / denom;
      }
      return result;
    }

  }; // class fixed_size_plane_base

}; // namespace sl

// I/O

template <class SELF_T, size_t DIMENSION, class T>
std::ostream& operator <<(std::ostream& s, const sl::fixed_size_plane_base<SELF_T,DIMENSION,T>& a) {
  for (size_t i = 0; i<a.size(); i++) {
      s << a(i) << std::endl;;
  }
  return s;
}
    
template <class SELF_T, size_t DIMENSION, class T>
std::istream& operator >>(std::istream& s, sl::fixed_size_plane_base<SELF_T,DIMENSION,T>& a) {
  for (size_t i = 0; i<a.size(); i++) {
      s >> a(i) ;
  }
  return s;
}


namespace sl {

  /// Planes of fixed dimension
  template <size_t DIMENSION, class T>
  class fixed_size_plane: 
    public fixed_size_plane_base< fixed_size_plane<DIMENSION, T>, DIMENSION, T > {

  public: // Types

    typedef fixed_size_plane_base< fixed_size_plane<DIMENSION, T>, DIMENSION, T > super_t;

    enum { dimension = super_t::dimension, hdimension = super_t::hdimension };

    typedef typename super_t::self_t        self_t;
    typedef typename super_t::value_t       value_t;
    typedef typename super_t::point_t       point_t;
    typedef typename super_t::vector_t      vector_t;
    typedef typename super_t::dual_vector_t dual_vector_t;
    typedef typename super_t::hpoint_t       hpoint_t;
    typedef typename super_t::hvector_t      hvector_t;
    typedef typename super_t::hdual_vector_t hdual_vector_t;

    typedef hdual_vector_t storage_t;

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_plane() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_plane(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }
   
    /**
     *  Initialize from homogeneous vector representation
     */
    inline explicit fixed_size_plane(const hdual_vector_t& hv): super_t(hv) {
    }

    /**
     *  Initialize from normal and point on plane
     */
    inline fixed_size_plane(const dual_vector_t& n,
			    const point_t& p): super_t(tags::not_initialized()) {
      this->storage_[dimension] = 0;
      for (size_t i=0; i<dimension; ++i) {
	this->storage_[i] = n[i];
	this->storage_[dimension] -= p[i]*n[i]; // -p dot n
      }
    }
       
    /**
     *  Initialize from normal and constant, defining the plane (n . P) + d = 0
     */
    inline fixed_size_plane(const dual_vector_t& n, value_t d): super_t(tags::not_initialized()) {
      for (size_t i=0; i<dimension; ++i) {
	this->storage_[i] = n[i];
      }
      this->storage_[dimension] = d;
    }
    
    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }
    
  };

}; // namespace sl


namespace sl {

  /// 2D plane with single precision floating plane components
  typedef fixed_size_plane<2,float> plane2f;
  /// 3D plane with single precision floating plane components
  typedef fixed_size_plane<3,float> plane3f;
  /// 4D plane with single precision floating plane components
  typedef fixed_size_plane<4,float> plane4f;

  /// 2D plane with double precision floating plane components
  typedef fixed_size_plane<2,double> plane2d;
  /// 3D plane with double precision floating plane components
  typedef fixed_size_plane<3,double> plane3d;
  /// 4D plane with double precision floating plane components
  typedef fixed_size_plane<4,double> plane4d;
}

#endif
