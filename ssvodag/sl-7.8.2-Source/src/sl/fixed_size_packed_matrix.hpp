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

#ifndef SL_FIXED_SIZE_PACKED_MATRIX_HPP
#define SL_FIXED_SIZE_PACKED_MATRIX_HPP

#include <sl/fixed_size_matrix.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/fixed_size_vector.hpp>

namespace sl {

  /// Symmetric matrix stored in packed triangular form (to be generalized!)
  template <size_t DIMENSION, class T>
  class fixed_size_packed_symmetric_matrix {
  public:
    enum { dimension = DIMENSION };
    enum { storage_count    = DIMENSION + (DIMENSION*(DIMENSION-1))/2 };
    
    typedef fixed_size_packed_symmetric_matrix<DIMENSION, T>   this_t;
    typedef this_t                                             transposed_t;
    typedef fixed_size_packed_symmetric_matrix<DIMENSION, T>   self_t;
    typedef T                                                  value_t;
    typedef fixed_size_array<storage_count, value_t>           storage_t;
    
    typedef fixed_size_square_matrix<dimension, value_t>        unpacked_matrix_t;

    typedef fixed_size_vector<column_orientation, dimension, T>   vector_t;
    typedef fixed_size_vector<row_orientation, dimension, T>      dual_vector_t;
    
  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", dimension > 0);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);

  protected:
    storage_t storage_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << storage_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> storage_;
    }

  public: // indexing
    
    /// the number of rows
    static inline size_t row_count() {
      return dimension;
    }
    
    /// the number of columns
    static inline size_t column_count() {
      return dimension;
    }
    
    /// is (i,j) a good component index?
    static inline bool good_index(size_t i, size_t j) {
      return 
        i < dimension &&
        j < dimension;
    }
    
    /// convert i,j into the position in a linear array
    static inline size_t array_index(size_t i, size_t j) {
      SL_REQUIRE("Good row index", i < dimension);
      SL_REQUIRE("Good column index", j < dimension);
      
      size_t result = (i<=j) ? (i+j*(j+1)/2) : (j+i*(i+1)/2);
      
      SL_ENSURE("Good storage index", result < storage_count);
      
      return result;
    }

  public: // Creation, Copy & Destruction
    
    /// Default initializer (zero)
    inline fixed_size_packed_symmetric_matrix() {
      // Storage is already 0-filled
    }
    
    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_packed_symmetric_matrix(tags::not_initialized tag): storage_(tag) {
      // Garbage in, use with care
    }
    
    /// Set this to other
    inline fixed_size_packed_symmetric_matrix(const this_t& other): storage_(other.storage_) {
    }

    /// Set this to other, converting from unpacked storage
    inline fixed_size_packed_symmetric_matrix(const unpacked_matrix_t& other): storage_(tags::not_initialized()) {
      for (size_t i=0; i<dimension; ++i) {
        for (size_t j=i; j<dimension; ++j) {
          (*this)(i,j) = other(i,j);
        }
      }
    }

    /// Set this to the outer product of a, b
    inline fixed_size_packed_symmetric_matrix(const vector_t& a, const vector_t& b): storage_(tags::not_initialized()) {
      for(int i=0; i<dimension; ++i) {
        for(int j=i; j<dimension; ++j) {
          (*this)(i,j) = a[i]*b[j];
        }
      }
    }

    /// This as a standard unpacked matrix
    inline unpacked_matrix_t unpacked() const {
      unpacked_matrix_t result;
      for (size_t i=0; i<dimension; ++i) {
        for (size_t j=i; j<dimension; ++j) {
          result(i,j) = result(j,i) = (*this)(i,j);
        }
      }
      return result;
    }

    /// Set this to other, converting from unpacked storage
    inline void from_unpacked(const unpacked_matrix_t& other) {
      for (size_t i=0; i<dimension; ++i) {
        for (size_t j=i; j<dimension; ++j) {
          (*this)(i,j) = other(i,j);
        }
      }
    }
    
  public: // Element Access
    
    /// the (i,j)-th element
    inline value_t& operator()(size_t i, size_t j) {
      SL_REQUIRE("Good index", good_index(i,j));
      return storage_[array_index(i,j)];
    }
      
    /// the (i,j)-th element
    inline value_t operator()(size_t i, size_t j) const {
      SL_REQUIRE("Good index", good_index(i,j));
      return storage_[array_index(i,j)];
    }
    
    /**
     *  Access to storage area.
     *  Storage is Packed Triangular (FORTRAN-LIKE, *NOT* C-LIKE!)
     */
    inline value_t *to_pointer() {
      return( storage_.to_pointer() );
    }
      
    /**
     *  Access to storage area.
     *  Storage is Packed Triangular (FORTRAN-LIKE, *NOT* C-LIKE!)
     *  since this class is used most often with OpenGL, which
     *  uses such a convention.
     */
    inline const value_t *to_pointer() const {
      return( storage_.to_pointer() );
    }
    
    /**
     *  Access to storage area.
     *  Storage is Column-Major (FORTRAN-LIKE, *NOT* C-LIKE!),
     *  since this class is used most often with OpenGL, which
     *  uses such a convention.
     */
    inline const fixed_size_array<storage_count, value_t>& storage() const {
      return storage_;
    } 
    
  public: // STL-style iterators
    
    typedef typename storage_t::iterator          iterator;
    typedef typename storage_t::restrict_iterator restrict_iterator;
    typedef typename storage_t::const_iterator    const_iterator;
    
    /// the number of components
    static inline size_t size()                { return storage_count; }
    /// linear iterator pointing to first component
    inline iterator begin()                    { return storage_.begin(); }
    /// linear restricted iterator pointing to first component
    inline iterator restrict_begin()  { return storage_.restrict_begin(); }
    /// linear const iterator pointing to first component
    inline const_iterator begin() const        { return storage_.begin(); }
    /// linear iterator pointing to end of storage
    inline iterator end()                      { return storage_.end(); }
    /// linear restricted iterator pointing to end of storage
    inline iterator restrict_end()    { return storage_.restrict_end(); }
    /// linear const iterator pointing to end of storage
    inline const_iterator end() const          { return storage_.end(); }
    
  public: // Comparison
    
    /// -1 if this < t2, +1 if this > t2, 0 otherwise (sequential element comparison)
    inline int compare(const self_t& t2) const {
      return storage().compare(t2.storage());
    }
    
    /// is this less than other? (sequential element comparison)
    inline bool operator<(const self_t& t2) const {
      return storage() < t2.storage();
    }
    
    /// is this equal to other?
    inline bool operator== (const self_t& t2) const {
      return storage() == t2.storage();
    }
    
    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);
    
    /// area all elements of this within eps of the corresponding elements of t2?
    inline bool is_epsilon_equal(const self_t& t2,
                                 const T& eps) const {
      return storage().is_epsilon_equal(t2.storage(), eps);
    }
    
  public: // Initializations
    
    /// set all components to s
    inline void fill (value_t s) {
      storage_.fill(s);
    }
      
    /// set all components to zero
    inline void to_zero() {
      fill(scalar_math<value_t>::zero());
    }

    inline void to_identity() {
      fill(scalar_math<value_t>::zero());
      fill_diagonal(scalar_math<value_t>::one());
    }
    
    /// set all diagonal components to s
    void fill_diagonal(value_t s) { 
      // Reimplement with iterators!
      size_t n = dimension;
      for (size_t i=0; i<n; i++) {
        (*this)(i,i) = s;
      }
    }
    
    /// the matrix with all components set to s
    static self_t constant(value_t s) {
      self_t result = tags::not_initialized();
      result.to_constant(s);
      return result;
    }
      
    /// the null matrix
    static self_t zero() {
      self_t result = tags::not_initialized();
      result.to_zero();
      return result;
    }

    
  public: // Arithmetic operations: Linear space (Matrices, Scalars)
    
    /// add other to this
    template <class VALUE_T2>
    inline self_t& operator+= (const fixed_size_packed_symmetric_matrix<DIMENSION, VALUE_T2>& other) {
      const typename fixed_size_packed_symmetric_matrix<DIMENSION, VALUE_T2>::const_iterator other_it = other.begin();
      for (size_t i=0; i<storage_count; ++i) {
        typedef SL_PROMOTENAME(value_t,VALUE_T2) arg_t;
        storage_[i] = value_t(arg_t(storage_[i]) + arg_t(other_it[i]));
      }
      return static_cast<self_t&>(*this);
    }
    
    /// subtract other from this
    template <class VALUE_T2>
    inline self_t& operator-= (const fixed_size_packed_symmetric_matrix<DIMENSION, VALUE_T2>& other) {
      const typename fixed_size_packed_symmetric_matrix<DIMENSION, VALUE_T2>::const_iterator other_it = other.begin();
      for (size_t i=0; i<storage_count; ++i) {
        typedef SL_PROMOTENAME(value_t,VALUE_T2) arg_t;
        storage_[i] = value_t(arg_t(storage_[i]) - arg_t(other_it[i]));
      }
      return static_cast<self_t&>(*this);
    }
      
    /// change sign to all components
    inline self_t operator- () const {
      self_t result = tags::not_initialized();
      restrict_iterator result_it = result.restrict_begin();
      for (size_t i=0; i<storage_count; ++i) {
        result_it[i] = - storage_[i];
      }
      return result;
    }
    
    /// multiply this by scalae
    template <class VALUE_T2>
    inline self_t& operator*= (VALUE_T2 scalar) {
      for (size_t i=0; i<storage_count; ++i) {
        typedef SL_PROMOTENAME(value_t,VALUE_T2) arg_t;
        storage_[i] = value_t(arg_t(storage_[i]) * arg_t(scalar));
      }
      return static_cast<self_t&>(*this);
    }
      
    /// divide this by scalar
    template <class VALUE_T2>
    inline self_t & operator/= (VALUE_T2 scalar) {
      return *this *= reciprocal(scalar);
    }
    
    SL_OP_LINEAR_SPACE(self_t,value_t); // Extend to generic types!
  public: // Interpolation

    /// Linear interpolation
    template <class T_PARAMETER>
    self_t lerp(const self_t& other, T_PARAMETER t) const {
      self_t result =  tags::not_initialized();
      //mbr: element_size not defined//fastest::transform<self_t::element_size>::apply(begin(), other.begin(), result.restrict_begin(), lerp_op<value_t,T_PARAMETER>(t));
      fastest::transform<self_t::storage_count>::apply(begin(), other.begin(), result.restrict_begin(), lerp_op<value_t,T_PARAMETER>(t));
      return result;
    }

  public: // Values

    /// the maximum absolute value of all components
    inline value_t amax() const {
      return sl::abs(storage_[storage_.iamax()]);
    }

    /// the minimum absolute value of all components
    inline value_t amin() const {
      return sl::abs(storage_[storage_.iamin()]);
    }
  public: // Division algebra operations

    /// this multiplied by other
    self_t operator *(const self_t& other) const {
      self_t result = tags::not_initialized();
      
      // General multiply
      for(int i=0; i < dimension; i++) {
        for(int j=i; j < dimension; j++) {
          T sum(0);
          for(int k=0; k < dimension; k++) {
            sum += (*this)(j,k) * other(k,i);
          }
          result(j,i) = sum;
        }           
      }
      return result;
    }

    /// post-multiply this by other
    self_t& operator *= (const self_t& other) {
      *this = *this * other;
      return static_cast<self_t&>(*this);
    }
    
    /// set mo to the inverse of this, and set ok to true if the operation was successfull
    void invert_to(self_t& mo, bool* ok) const {
      // FIXME!!
      unpacked_matrix_t mi2 = unpacked();
      unpacked_matrix_t mo2 = tags::not_initialized();
      mi2.invert_to(mo2, ok);
      if (ok) {
        mo.from_unpacked(mo2);
      }
    }

    /// The inverse of other (requires invertible!)
    self_t inverse() const {
      bool ok;
      self_t result = tags::not_initialized();
      invert_to(result, &ok);
      SL_CHECK("Invertible", ok);
      return result;
    }

    /// This inverse of other (requires invertible!)
    inline self_t operator~() const {
      return inverse();
    }                                                                                                                    
    
    /// Post-multiply by the inverse of other
    inline self_t& operator /= (const self_t& other) {
      *this *= other.inverse();
      return static_cast<self_t&>(*this);
    }

  }; // class fixed_size_packed_symmetric_matrix

} // namespace sl


// Arithmetic operators overloads

template <size_t DIMENSION, class T>
inline sl::fixed_size_packed_symmetric_matrix<DIMENSION,T> operator*(const T& y, const sl::fixed_size_packed_symmetric_matrix<DIMENSION,T>& x) {
  typedef sl::fixed_size_packed_symmetric_matrix<DIMENSION,T> self_t;

  self_t result = sl::tags::not_initialized();
   
  const typename self_t::restrict_iterator result_it = result.restrict_begin();
  const typename self_t::const_iterator    src_it    = x.begin();
  for (size_t i=0; i<self_t::element_size; ++i) {
    result_it[i] = src_it[i] * y;
  }
  return result;
}
// ---------------------------------------------------------------------------
// Matrix/Vector operations overloads
// ---------------------------------------------------------------------------

// TODO: HANDLE TYPE CONVERSIONS!!!

template <size_t DIMENSION, class T>
sl::fixed_size_vector<sl::row_orientation, DIMENSION, T> operator *(const sl::fixed_size_vector<sl::row_orientation, DIMENSION, T>& vec,
                                                                const sl::fixed_size_packed_symmetric_matrix<DIMENSION,T>& mat) {
  sl::fixed_size_vector<sl::row_orientation, DIMENSION, T> result = sl::tags::not_initialized();

  for (size_t i=0; i<DIMENSION; ++i) {
    result[i] = vec[0] * mat(0,i);
    for (size_t j=1; j<DIMENSION; ++j) {
      result[i] += vec[j] * mat(j,i);
    }
  }
  return result;  
}

/// the result of matrix * vector
template <size_t DIMENSION, class T>
sl::fixed_size_vector<sl::column_orientation, DIMENSION, T> operator *(const sl::fixed_size_packed_symmetric_matrix<DIMENSION,T>& mat,
                                                                   const sl::fixed_size_vector<sl::column_orientation, DIMENSION, T>& vec) {
  sl::fixed_size_vector<sl::column_orientation, DIMENSION, T> result = sl::tags::not_initialized();

  for (size_t i=0; i<DIMENSION; ++i) {
    result[i] = mat(i,0) * vec[0];
    for (size_t j=1; j<DIMENSION; ++j) {
      result[i] += mat(i,j) * vec[j];
    }
  }  

  return result;  
}

// ---------------------------------------------------------------------------

namespace sl {

  /// A factory for packed_symmetric matrices
  template <size_t DIMENSION, class T>
  class fixed_size_packed_symmetric_matrix_factory {
    
  public: // Constants and types 

    enum { dimension = DIMENSION };

    typedef fixed_size_packed_symmetric_matrix<DIMENSION, T> matrix_t;
    typedef T                                 value_t;

  public: // Creators

    /// The identity matrix
    static matrix_t identity() {
      matrix_t result;
      result.fill_diagonal(scalar_math<value_t>::one());
      return result;
    }

    /// The null matrix
    static matrix_t zero() {
      value_t result;
      return result;
    }

  }; // class fixed_size_packed_symmetric_matrix_factory
      
}; // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {

  /// 2D packed symmetric matrix with single precision floating point components
  typedef fixed_size_packed_symmetric_matrix<2,float> packed_symmetric_matrix2f;
  /// 3D packed symmetric matrix with single precision floating point components
  typedef fixed_size_packed_symmetric_matrix<3,float> packed_symmetric_matrix3f;
  /// 4D packed symmetric matrix with single precision floating point components
  typedef fixed_size_packed_symmetric_matrix<4,float> packed_symmetric_matrix4f;

  /// 2D packed symmetric matrix with double precision floating point components
  typedef fixed_size_packed_symmetric_matrix<2,double> packed_symmetric_matrix2d;
  /// 3D packed symmetric matrix with double precision floating point components
  typedef fixed_size_packed_symmetric_matrix<3,double> packed_symmetric_matrix3d;
  /// 4D packed symmetric matrix with double precision floating point components
  typedef fixed_size_packed_symmetric_matrix<4,double> packed_symmetric_matrix4d;

};

#endif
