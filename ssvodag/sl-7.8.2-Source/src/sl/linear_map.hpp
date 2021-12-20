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
#ifndef LINEAR_MAP_HPP
#define LINEAR_MAP_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/fixed_size_plane.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/utility.hpp>

namespace sl {

  /*
   * Map property tags
   */
  struct identity_tag {
    enum { is_identity = true, is_rigid_body = true, is_affine = true };
  };
  struct rigid_body_tag {
    enum { is_identity = false, is_rigid_body = true, is_affine = true };
  };
  struct affine_tag {
    enum { is_identity = false, is_rigid_body = false, is_affine = true };
  };
  struct general_tag {
    enum { is_identity = false, is_rigid_body = false, is_affine = false };
  };

  template <class T_tag_left, class T_tag_right> 
  struct linear_map_composed_tag {
    typedef general_tag tag;
  };

  template <class T_tag>
  struct linear_map_composed_tag<T_tag, T_tag> { 
    typedef T_tag tag;
  };

  template <>
  struct linear_map_composed_tag<identity_tag, rigid_body_tag> { 
    typedef rigid_body_tag tag;
  };

  template <>
  struct linear_map_composed_tag<identity_tag, affine_tag> { 
    typedef affine_tag tag;
  };

  template <>
  struct linear_map_composed_tag<rigid_body_tag, identity_tag> { 
    typedef rigid_body_tag tag;
  };

  template <>
  struct linear_map_composed_tag<rigid_body_tag, affine_tag> { 
    typedef affine_tag tag;
  };

  template <>
  struct linear_map_composed_tag<affine_tag, identity_tag> { 
    typedef affine_tag tag;
  };

  template <>
  struct linear_map_composed_tag<affine_tag, rigid_body_tag> { 
    typedef affine_tag tag;
  };
  
  /* 
   * Forward declarations
   */
  template <size_t N_dimension, class T_value> 
  class identity_map;

  template <size_t N_dimension, class T_value> 
  class rigid_body_map;

  template <size_t N_dimension, class T_value> 
  class affine_map;

  template <size_t N_dimension, class T_value> 
  class projective_map;

  /**
   *  The actual type corresponding to given properties
   */
  template <class T_tag, size_t N_dimension, class T_value>
  class linear_map_generator {
  public:
    typedef projective_map<N_dimension, T_value> type;
  };

  template <size_t N_dimension, class T_value>
  class linear_map_generator<identity_tag, N_dimension, T_value> {
  public:
    typedef identity_map<N_dimension, T_value> type;
  };

  template <size_t N_dimension, class T_value>
  class linear_map_generator<rigid_body_tag, N_dimension, T_value> {
  public:
    typedef rigid_body_map<N_dimension, T_value> type;
  };

  template <size_t N_dimension, class T_value>
  class linear_map_generator<affine_tag, N_dimension, T_value> {
  public:
    typedef affine_map<N_dimension, T_value> type;
  };

  /**
   *  The actual type corresponding to given properties
   */
  template <
    class T_tag_left,
    class T_tag_right,
    size_t N_dimension, 
    class T_value>
  class linear_map_mul_generator {
  public:
    typedef 
      typename linear_map_generator< 
        typename linear_map_composed_tag<T_tag_left, T_tag_right>::tag, N_dimension, T_value >::type type;
  };

  /**
   *  Generic base class for linear maps, i.e. objects that perform
   *  linear transformation of geometric entities. The internal
   *  representation is an homogeneous matrix.
   */
  template < 
    class T_tag,
    size_t N_dimension, 
    class T_value>
  class linear_map_base {
  public: // Constants and types

    enum { dimension = N_dimension };

    typedef T_tag                                                                                   tag_t;
    typedef typename linear_map_generator<T_tag, N_dimension, T_value>::type                        concrete_t;
    typedef typename linear_map_generator<identity_tag, N_dimension, T_value>::type                 identity_t;
    typedef T_value                                                                                 value_t;

    typedef fixed_size_vector<column_orientation, N_dimension, T_value>                             vector_t;
    typedef fixed_size_vector<row_orientation, N_dimension, T_value>                                dual_vector_t;
    typedef fixed_size_point<N_dimension, T_value>                                                  point_t;
    typedef fixed_size_plane<N_dimension, T_value>                                                  plane_t;

    typedef fixed_size_square_matrix<N_dimension+1, T_value>                                        matrix_t;
    
    enum { is_identity   = tag_t::is_identity };
    enum { is_rigid_body = tag_t::is_rigid_body };
    enum { is_affine     = tag_t::is_affine };

  protected: // Storage

    matrix_t storage_;

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", N_dimension > 0);
    SL_COMPILE_TIME_CHECK("Rigid body implies affine", !is_rigid_body || is_affine);
    SL_COMPILE_TIME_CHECK("Identity implies rigid_body", !is_identity || is_rigid_body);

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << storage_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> storage_;
    }

  public: // Invariant

    /// Is the matrix m representing a compatible transformation?
    static bool is_compatible(const matrix_t& m) {
      bool result = true;
      if (is_identity) {
        result = m.is_identity();
      } else if (is_rigid_body) {
        result = m.is_affine(); // TODO: m.is_scale_preserving();
      } else if (is_affine) {
        result = m.is_affine();
      }
      return result;
    }
    
    /// Is this in a coherent state?
    inline bool invariant() const {
      return is_compatible(storage_);
    }

  protected: // Creation, Copy & Destruction

    /// Default init (identity)
    inline linear_map_base(): storage_(tags::not_initialized()) {
      storage_.to_identity();
    }

    /// Fast init (garbage), handle with care
    inline linear_map_base(tags::not_initialized tag): storage_(tag) {
      // Garbage in, use with care
    }

    /// Init from matrix m, assumed compatible
    inline linear_map_base(const matrix_t& m): storage_(m) {
      SL_REQUIRE("Compatible matrix", is_compatible(m));
      SL_INVARIANT(invariant());
    }

  public: // Element Access (read-only)

    /// Is (i,j) a good element index for the homogeneous matrix representation of this?
    inline static bool good_index(size_t i, size_t j) {
      return matrix_t::good_index(i,j);
    }

    /// The (i,j)-th component of the homogeneous matrix representation of this
    inline value_t operator()(size_t i, size_t j) const {
      SL_REQUIRE("Good index", good_index(i,j));
      return storage_(i,j);
    }
    
    /// Pointer to storage area (homogeneous matrix representation of this)
    inline const value_t *to_pointer() const {
      return storage_.to_pointer();
    }

    /// Storage area (homogeneous matrix representation of this)
    inline const matrix_t& storage() const {
      return storage_;
    } 

    /// The homogeneous matrix representation of this
    inline const matrix_t& as_matrix() const {
      return storage_;
    } 

    /// The i-th axis of the basis
    inline vector_t axis(size_t i) const {
      SL_REQUIRE("Good index", i<dimension);

      vector_t result = tags::not_initialized();
      for (size_t j=0; j<dimension; ++j) {
        result[j] = storage_(i,j);
      }
      return result;
    }
      
  public: // Comparison

    /// -1 if this < t2, +1, if this > t2, 0 otherwise (sequential element comparison)
    template <class T_other_tag>
    inline int compare(const linear_map_base<
		       T_other_tag,
		       dimension,
		       value_t>& t2) const {
      return storage_.compare(t2.storage_);
    }

    /// is this less than t2? (sequential element comparison)
    template <class T_other_tag>
    inline bool operator< (const linear_map_base<
			   T_other_tag,
			   dimension,
			   value_t>& t2) const {
      return storage_ < t2.storage_;
    }

    /// is this equal to t2?
    template <class T_other_tag>
    inline bool operator== (const linear_map_base<
			    T_other_tag,
			    dimension,
			    value_t>& t2) const {
      return storage_ == t2.storage_;
    }

    /// is this greater than t2?
    template <class T_other_tag>
    inline bool operator> (const linear_map_base<
			   T_other_tag,
			   dimension,
			   value_t>& t2)  const { 
      return t2 < *this; 
    }

    /// is this less than or equal to t2?
    template <class T_other_tag>
    inline bool operator<= (const linear_map_base<
			    T_other_tag,
			    dimension,
			    value_t>& t2) const { 
      return !(t2 < *this);
    }
			    
    /// is this greater than or equal to t2?
    template <class T_other_tag>
    inline bool operator>= (const linear_map_base<
			    T_other_tag,
			    dimension,
			    value_t>& t2) const { 
      return !(*this < t2); 
    }

    /// is this not equal to t2?
    template <class T_other_tag>
    inline bool operator!= (const linear_map_base<
			    T_other_tag,
			    dimension,
			    value_t>& t2) const {
      return storage_ != t2.storage_;
    }

    /// Are all elements of the homogenous matrix representation of this and t2 within eps?
    template <class T_other_tag>
    inline bool is_epsilon_equal(const linear_map_base<
				 T_other_tag,
				 dimension,
				 value_t>& t2,
         value_t eps) const {
      SL_USEVAR(eps);
      return storage_.is_epsilon_equal(t2.storage_);
    }

  public: // Initializations
    
    /// Set this to the identity map
    inline void to_identity(void) {
      storage_.to_identity();
    }

    /// The identity map
    static const identity_t& identity();

    /// Set this to matrix m (required compatible)
    inline void set_matrix(const matrix_t & m) {
      SL_REQUIRE("Compatible matrix", is_compatible(m));
      storage_ = m;
    }

  public: /// Inverse

    /// Invert this to mo, and set ok to true iff the operation was successful.
    void invert_to(concrete_t& mo, bool* ok) const {
      if (is_identity) {
	// Nothing to do
	*ok = true;
      } else  if (is_rigid_body) {
	// Transpose rotation matrix
	for (size_t j=0; j<dimension; ++j) {
	  for (size_t i=0; i<dimension; ++i) {
	    mo.storage_(i,j) = storage_(j,i);
	  }  
	  mo.storage_(dimension,j) = sl::zero(value_t());
	}
	mo.storage_(dimension,dimension) = sl::one(value_t());

	// Invert translation
	for (size_t i=0; i<dimension; ++i) {
	  mo.storage_(i,dimension) = sl::zero(value_t());
	  for (size_t j=0; j<dimension; ++j) {
	    mo.storage_(i,dimension) -= storage_(j, i) * storage_(j, dimension);
	  }
	}
	// Signal success
	*ok = true;
      } else if (is_affine) {
	// TODO: Optimize affine case
	storage_.invert_to(mo.storage_, ok);
      } else {
	storage_.invert_to(mo.storage_, ok);
      }
    }
      
    /// The inverse map
    inline concrete_t inverse() const {
      bool ok;
      concrete_t result = tags::not_initialized();
      invert_to(result, &ok);
      SL_CHECK("Invertible", ok);
      return result;
    }

    /// the inverse map
    inline concrete_t operator~() const {
      return inverse();
    }                                                                                                                    

  public: // Multiplicative group operations

    /// Post-multiply this by other
    template <class T_other_tag>
    inline concrete_t& operator*= (const linear_map_base<
				   T_other_tag,
				   dimension,
				   value_t>& other) {
      SL_COMPILE_TIME_CHECK("Other is compatible", 
			    (!is_identity || T_other_tag::is_identity) &&
			    (!is_rigid_body || T_other_tag::is_rigid_body) &&
			    (!is_affine || T_other_tag::is_affine));
      // FIXME *this = (*this) * (other);
      storage_ = as_matrix() * other.as_matrix(); // FIXME!!!!
      return static_cast<concrete_t&>(*this);
    }


    /// Post-multiply this by the inverse of other
    template <class T_other_tag>
    inline concrete_t& operator /= (const linear_map_base<
				    T_other_tag,
				   dimension,
				   value_t>& other) {
      SL_COMPILE_TIME_CHECK("Other is compatible", 
			    (!is_identity || T_other_tag::is_identity) &&
			    (!is_rigid_body || T_other_tag::is_rigid_body) &&
			    (!is_affine || T_other_tag::is_affine));
      *this *= other.inverse();
      return static_cast<concrete_t&>(*this);
    }
    
    /// The result of this multiplied by the inverse of other
    template <class T_other_tag>
    inline typename linear_map_mul_generator<
      T_tag, T_other_tag,
      dimension,
      value_t>::type operator/ (const linear_map_base<
				T_other_tag,
				dimension,
				value_t>& other) const {
      return *this * other.inverse();
    }


  public: // Interpolation

    /// linear interpolation
    concrete_t lerp(const concrete_t& f2, value_t t) const {
      if (is_identity) {
        return *this;
      } else if (is_rigid_body || is_affine) {
        SL_FAIL("Interpolation not implemented for N-dimensional rigid body or affine matrices!");
      } else {
        return concrete_t(storage_ * (sl::one(t)-t) + f2.storage_ * (t));
      }
    }
    

  }; // linear_map_base

}; // namespace sl

// --------------------------------------------------------------------
// Multiplicative operators
// --------------------------------------------------------------------
// TODO: Other special cases

template 
< std::size_t N, class T, 
  class T_tag_left, class T_tag_right>
inline 
typename sl::linear_map_mul_generator<T_tag_left, T_tag_right,N,T>::type 
operator*(const sl::linear_map_base<T_tag_left, N, T>& l,
	  const sl::linear_map_base<T_tag_right, N, T>& r) {
  typedef typename sl::linear_map_mul_generator<T_tag_left, T_tag_right,N,T>::type result_t;
  // Any * any
  return result_t(l.as_matrix() * r.as_matrix());
}

template 
< std::size_t N, class T, 
  class T_tag_right>
inline 
typename sl::linear_map_mul_generator<sl::identity_tag, T_tag_right,N,T>::type 
operator*(const sl::linear_map_base<sl::identity_tag, N, T>&,
	  const sl::linear_map_base<T_tag_right, N, T>& r) {
  return r;
}

template 
< std::size_t N, class T, 
  class T_tag_left >
inline 
typename sl::linear_map_mul_generator<T_tag_left, sl::identity_tag,N,T>::type 
operator*(const sl::linear_map_base<T_tag_left, N, T>& l,
	  const sl::linear_map_base<sl::identity_tag, N, T>&) {
  return l;
}

// --------------------------------------------------------------------
// I/O
// --------------------------------------------------------------------

/// Write linear map to output stream
template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
std::ostream& operator <<(std::ostream& s, 
			  const sl::linear_map_base<
			  T_tag,
			  N_dimension,
			  T_value>& a) {
  s << a.as_matrix(); 
  return s;
}
    
/// Read linear map from input stream
template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
std::istream& operator >>(std::istream& s, 
			  sl::linear_map_base<
			  T_tag,
			  N_dimension,
			  T_value>& a) {
  typename sl::linear_map_base<	
    T_tag,
    N_dimension,T_value>::matrix_t m;
  s >> m;
  a = m;
  return s;
}

// --------------------------------------------------------------------
// Geometric entity transformations - points
// --------------------------------------------------------------------

/// The point p transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  SL_USEVAR(m);
  return p;
}

/// The point p transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  sl::fixed_size_point<N_dimension, T_value> result = sl::tags::not_initialized();
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(i,N_dimension);
    for (size_t j=0; j<N_dimension; ++j) {
      result[i] += m(i,j) * p[j];
    }
  }
  return result;
}

/// The point p transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  sl::fixed_size_point<N_dimension, T_value> result = sl::tags::not_initialized();
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(i,N_dimension);
    for (size_t j=0; j<N_dimension; ++j) {
      result[i] += m(i,j) * p[j];
    }
  }
  return result;
}

/// The point p transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  return from_homogeneous( as_point(m.as_matrix() * as_homogeneous(p).as_vector()) );
}

// --------------------------------------------------------------------
// Geometric entity transformations - vectors
// --------------------------------------------------------------------


/// The vector v transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  SL_USEVAR(m);
  return v;
}

/// The vector v transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value> result = sl::tags::not_initialized();
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(i,0) * v[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] += m(i,j) * v[j];
    }
  }
  return result;
}

/// The vector v transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value> result = sl::tags::not_initialized();
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(i,0) * v[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] += m(i,j) * v[j];
    }
  }
  return result;
}

/// The vector v transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  return from_homogeneous( m.as_matrix() * as_homogeneous(v) );
}

// --------------------------------------------------------------------
// Geometric entity transformations - dual vectors
// --------------------------------------------------------------------


/// The dual vector v transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  SL_USEVAR(m);
  return dual_v;
}

/// The dual vector v transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  // v * inv(R) = v_i * trans(R)
  sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value> result = sl::tags::not_initialized();
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(i,0) * dual_v[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] +=  dual_v[j] * m(i,j);
    }
  }
  return result;
}

/// The dual vector v transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  // v * inv(M)
  // TODO: Optimize this!
  sl::fixed_size_vector<sl::row_orientation, N_dimension+1, T_value> vh = as_homogeneous(dual_v) * m.inverse().as_matrix();
  vh[N_dimension] = sl::scalar_math<T_value>::zero();
  return from_homogeneous( vh );
}

/// The dual vector v transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  // v * inv(M)
  // TODO: Optimize this!
  sl::fixed_size_vector<sl::row_orientation, N_dimension+1, T_value> vh = as_homogeneous(dual_v) * m.inverse().as_matrix();
  vh[N_dimension] = sl::scalar_math<T_value>::zero();
  return from_homogeneous( vh );
}

// --------------------------------------------------------------------
// Geometric entity transformations - hyperplanes
// --------------------------------------------------------------------

/// The hyperplane hp transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  SL_USEVAR(m);
  return hp;
}

/// The hyperplane hp transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  sl::fixed_size_plane<N_dimension, T_value> result = sl::tags::not_initialized();

  // v * inv(TR) = v * inv(R) * inv(T) = v * trans(R) * inv(T)
  
  // Multiply by transpose of rotation
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = hp[0] * m(i,0);
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] +=  hp[j] * m(i,j);
    }
  }
  result[N_dimension] = hp[N_dimension];
  
  // Multiply by inverse of translation
  for (size_t i=0; i<N_dimension; ++i) {
    result[N_dimension] -= result[i] * m(i,N_dimension);
  }

  return result;
}

/// The hyperplane hp transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  return sl::fixed_size_plane<N_dimension, T_value>(hp.as_vector() * m.inverse().as_matrix()); 
}

/// The hyperplane hp transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& dual_v) {
  return sl::fixed_size_plane<N_dimension, T_value>(dual_v.as_vector() * m.inverse().as_matrix()); 
}

// --------------------------------------------------------------------
// Geometric entity inverse transformations - points
// --------------------------------------------------------------------

/// The point p transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
		       const sl::fixed_size_point<N_dimension, T_value>& p) {
  SL_USEVAR(m);
  return p;
}

/// The point p transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  sl::fixed_size_point<N_dimension, T_value> result = sl::tags::not_initialized();
  sl::fixed_size_point<N_dimension, T_value> tmp    = sl::tags::not_initialized();

  for (size_t i=0; i<N_dimension; ++i) {
    tmp[i] = p[i] - m(i,N_dimension);
  }
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(0,i) * tmp[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] += m(j,i) * tmp[j];
    }
  }
  return result;
}

/// The point p transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  return transformation(m.inverse(), p); // no optimization
}

/// The point p transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_point<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_point<N_dimension, T_value>& p) {
  return transformation(m.inverse(), p); // no optimization
}

// --------------------------------------------------------------------
// Geometric entity inverse transformations - vectors
// --------------------------------------------------------------------


/// The vector v transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  SL_USEVAR(m);
  return v;
}

/// The vector v transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value> result = sl::tags::not_initialized();
  // inv(R) * v = trans(R) * v
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(0,i) * v[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] += m(j,i) * v[j];
    }
  }
  return result;
}

/// The vector v transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  return transformation(m.inverse(), v); // no optimization
}

/// The vector v transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::column_orientation, N_dimension, T_value>& v) {
  return transformation(m.inverse(), v); // no optimization
}

// --------------------------------------------------------------------
// Geometric entity inverse transformations - dual vectors
// --------------------------------------------------------------------


/// The dual vector v transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  SL_USEVAR(m);
  return dual_v;
}

/// The dual vector v transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
	
  sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value> result = sl::tags::not_initialized();
  // v * inv(inv(M)) = v * M
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(0,i) * dual_v[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] +=  dual_v[j] * m(j,i);
    }
  }
  return result;
}

/// The dual vector v transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value> result = sl::tags::not_initialized();
  // v * inv(inv(M)) = v * M
  for (size_t i=0; i<N_dimension; ++i) {
    result[i] = m(0,i) * dual_v[0];
    for (size_t j=1; j<N_dimension; ++j) {
      result[i] +=  dual_v[j] * m(j,i);
    }
  }
  return result;
}

/// The dual vector v transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_vector<sl::row_orientation, N_dimension, T_value>& dual_v) {
  // v * inv(inv(M)) = v * M
  // TODO: Optimize this!
  sl::fixed_size_vector<sl::row_orientation, N_dimension+1, T_value> vh = as_homogeneous(dual_v) * m.as_matrix();
  vh[N_dimension] = sl::scalar_math<T_value>::zero();
  return from_homogeneous( vh );
}

// --------------------------------------------------------------------
// Geometric entity inverse transformations - hyperplanes
// --------------------------------------------------------------------

/// The hyperplane hp transformed by the identity map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::identity_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  SL_USEVAR(m);
  return hp;
}

/// The hyperplane hp transformed by the rigid body map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::rigid_body_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  return sl::fixed_size_plane<N_dimension, T_value>(hp.as_vector() * m.as_matrix());
}

/// The hyperplane hp transformed by the affine map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::affine_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  return sl::fixed_size_plane<N_dimension, T_value>(hp.as_vector() * m.as_matrix());
}

/// The hyperplane hp transformed by the general map m
template <size_t N_dimension, class T_value> 
inline sl::fixed_size_plane<N_dimension, T_value>
inverse_transformation(const sl::linear_map_base<sl::general_tag, N_dimension, T_value>& m,
	       const sl::fixed_size_plane<N_dimension, T_value>& hp) {
  return sl::fixed_size_plane<N_dimension, T_value>(hp.as_vector() * m.as_matrix());
}
  
// --------------------------------------------------------------------
// Operator overloading
// --------------------------------------------------------------------

/// ---  point transformation
template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
inline sl::fixed_size_point<N_dimension,T_value> operator * (const sl::linear_map_base<T_tag, N_dimension, T_value>& f, 
							     const sl::fixed_size_point<N_dimension,T_value>& p) {
  return transformation(f,p);
}

/// ---- vector transformation
template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
inline sl::fixed_size_vector<sl::column_orientation,N_dimension,T_value> operator * (const sl::linear_map_base<T_tag, N_dimension, T_value>& f, 
									       const sl::fixed_size_vector<sl::column_orientation,N_dimension,T_value>& v) {
  return transformation(f,v);
}

/// ---  dual vector transformation
template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
inline sl::fixed_size_vector<sl::row_orientation, N_dimension,T_value> operator * (const sl::linear_map_base<T_tag, N_dimension, T_value>& f, 
										   const sl::fixed_size_vector<sl::row_orientation,N_dimension,T_value>& dual_v) {
  return transformation(f,dual_v);
}

/// --- (hyper)plane transformation
template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
inline sl::fixed_size_plane<N_dimension,T_value> operator * (const sl::linear_map_base<T_tag, N_dimension, T_value>& f, 
						       const sl::fixed_size_plane<N_dimension,T_value>& hp) {
  return transformation(f, hp);
}

// --------------------------------------------------------------------
// Identity
// --------------------------------------------------------------------

#include <sl/identity_map.hpp>

template < 
  class T_tag,
  size_t N_dimension, 
  class T_value>
const typename sl::linear_map_base<T_tag, N_dimension, T_value>::identity_t& 
sl::linear_map_base<T_tag, N_dimension, T_value>::identity() {
  typedef typename sl::linear_map_base<T_tag, N_dimension, T_value>::identity_t result_t;
  static const result_t I_;
  return I_;
}

#endif
