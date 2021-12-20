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
#ifndef SL_FIXED_SIZE_MATRIX_HPP
#define SL_FIXED_SIZE_MATRIX_HPP

#include <sl/conv_to.hpp>
#include <sl/fixed_size_array.hpp>
#include <sl/iterator.hpp>
#include <sl/numeric_traits.hpp>

// --------------------------------------------------------------------
// helpers
// --------------------------------------------------------------------

namespace sl {

  /// Function object, x+y
  template <class arg1_t, class arg2_t, class result_t>
  class plus_op : public std::binary_function<arg1_t,arg2_t,result_t> {
  public:
    inline result_t operator()(arg1_t x, arg2_t y) const { 
      typedef SL_PROMOTENAME(arg1_t,arg2_t) arg_t;
      return result_t(arg_t(x)+arg_t(y)); 
    }
  };

  /// Function object, x-y
  template <class arg1_t, class arg2_t, class result_t>
  class minus_op : public std::binary_function<arg1_t,arg2_t,result_t> {
  public:
    inline result_t operator()(arg1_t x, arg2_t y) const { 
      typedef SL_PROMOTENAME(arg1_t,arg2_t) arg_t;
      return result_t(arg_t(x)-arg_t(y)); 
    }
  };

  /// Function object, y-x
  template <class arg1_t, class arg2_t, class result_t>
  class rev_minus_op : public std::binary_function<arg1_t,arg2_t,result_t> {
  public:
    inline result_t operator()(arg1_t x, arg2_t y) const { 
      typedef SL_PROMOTENAME(arg1_t,arg2_t) arg_t;
      return result_t(arg_t(y)-arg_t(x)); 
    }
  };

  /// Function object, s*x
  template <class arg1_t, class arg2_t, class result_t>
  class scale_op : public std::unary_function<arg1_t,result_t> {
  public:
    typedef SL_PROMOTENAME(arg1_t,arg2_t) arg_t;
    const arg_t s;
    inline scale_op(arg2_t s_arg): s(arg_t(s_arg)) {};
    inline result_t operator()(arg1_t x) const { 
      return result_t(arg_t(x)*s); 
    }
  };

  /// Function object, -x
  template <class arg1_t, class result_t>
  class negate_op : public std::unary_function<arg1_t,result_t> {
  public:
    inline result_t operator()(arg1_t x) const { 
      return result_t(-x); 
    }
  };
 
  /// Function object, x+y^2
  template <class arg1_t, class arg2_t, class result_t>
  class add_sqr_op : public std::binary_function<arg1_t,arg2_t,result_t> {
  public:
    typedef SL_PROMOTENAME(arg1_t,arg2_t) arg_t;
    inline result_t operator()(arg1_t x, arg2_t y) const { 
      return result_t(arg_t(x) + sqr(arg_t(y))); 
    }
  };

} // namespace sl

namespace sl {

  /// Base class for rectangular matrices of fixed dimension
  template <class SELF_T, class TRANSPOSED_T, size_t N_ROW, size_t N_COL, class T>
  class fixed_size_matrix_base {
  public: // Constants and types

    enum { element_size    = N_ROW*N_COL };

    typedef fixed_size_matrix_base<SELF_T, TRANSPOSED_T, N_ROW, N_COL, T> this_t;
    typedef TRANSPOSED_T                  transposed_t;
    typedef SELF_T                        self_t;
    typedef T                             value_t;
    typedef fixed_size_array<element_size, value_t> storage_t;

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Non null rows",    N_ROW > 0);
    SL_COMPILE_TIME_CHECK("Non null columns", N_COL > 0);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);

  protected: // Storage 
    
    storage_t storage_;

  public: // indexing

    /// the number of rows
    static inline size_t row_count() {
      return N_ROW;
    }

    /// the number of columns
    static inline size_t column_count() {
      return N_COL;
    }

    /// is (i,j) a good component index?
    static inline bool good_index(size_t i, size_t j) {
      return 
	i < N_ROW &&
	j < N_COL;
    }
	
    /// convert i,j into the position in a linear array
    static inline size_t array_index(size_t i, size_t j) {
      SL_REQUIRE("Good row index", i < N_ROW);
      SL_REQUIRE("Good column index", j < N_COL);

      size_t result = i + j * N_ROW;

      SL_ENSURE("Good storage index", result < N_ROW*N_COL);
      
      return result;
    }

  protected: // Creation, Copy & Destruction

    /// Default initializer (zero)
    inline fixed_size_matrix_base() {
      // Storage is already 0-filled
    }

    /// Fast initializer, values set to garbage. handle with care!
    inline fixed_size_matrix_base(tags::not_initialized tag): storage_(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_matrix_base(const this_t& other): storage_(other.storage_) {
    }

  public: // Initialization

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as m = 1.0; (fill) and 
     *  m = 1.0, 2.0, 3.0, 4.0; (component init)
     */
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, N_ROW, N_COL, v);
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
     *  Storage is Column-Major (FORTRAN-LIKE, *NOT* C-LIKE!),
     *  since this class is used most often with OpenGL, which
     *  uses such a convention.
     */
    inline value_t *to_pointer() {
      return( storage_.to_pointer() );
    }
    
    /**
     *  Access to storage area.
     *  Storage is Column-Major (FORTRAN-LIKE, *NOT* C-LIKE!),
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
    inline const fixed_size_array<element_size, value_t>& storage() const {
      return storage_;
    } 

  public: // STL-style iterators

    typedef typename storage_t::iterator          iterator;
    typedef typename storage_t::restrict_iterator restrict_iterator;
    typedef typename storage_t::const_iterator    const_iterator;

    /// the number of components
    static inline size_t size()                { return element_size; }
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

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << storage_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> storage_;
    }

  public: // MTL-style iterators

    typedef fixed_strided_2d_accessor< iterator, N_ROW, N_COL, 1, N_ROW> column_accessor_t;
    typedef fixed_strided_2d_accessor< iterator, 1, N_ROW, N_ROW, N_COL> row_accessor_t;

    typedef fixed_strided_2d_accessor< const_iterator, N_ROW, N_COL, 1, N_ROW> const_column_accessor_t;
    typedef fixed_strided_2d_accessor< const_iterator, 1, N_ROW, N_ROW, N_COL> const_row_accessor_t;
    
    typedef fixed_strided_2d_accessor< restrict_iterator, N_ROW, N_COL, 1, N_ROW> restrict_column_accessor_t;
    typedef fixed_strided_2d_accessor< restrict_iterator, 1, N_ROW, N_ROW, N_COL> restrict_row_accessor_t;

    typedef fixed_strided_1d_accessor< iterator, 1, N_ROW > column_t;
    typedef fixed_strided_1d_accessor< iterator, N_ROW, N_COL > row_t;

    typedef fixed_strided_1d_accessor< const_iterator, 1, N_ROW > const_column_t;
    typedef fixed_strided_1d_accessor< const_iterator, N_ROW, N_COL > const_row_t;

    typedef fixed_strided_1d_accessor< restrict_iterator, 1, N_ROW > restrict_column_t;
    typedef fixed_strided_1d_accessor< restrict_iterator, N_ROW, N_COL > restrict_row_t;

    /// Accessor to matrix columns
    inline column_accessor_t columns() {
      return column_accessor_t(begin());
    }

    /// Accessor to matrix columns
    inline const_column_accessor_t columns() const {
      return const_column_accessor_t(begin());
    }

    /// Accessor to matrix columns
    inline restrict_column_accessor_t restrict_columns() {
      return restrict_column_accessor_t(restrict_begin());
    }

    /// The i-th column
    inline column_t column(size_t i) {
      SL_REQUIRE("Good index", i < N_COL);
      return columns()[i];
    }

    /// The i-th column
    inline const_column_t column(size_t i) const {
      SL_REQUIRE("Good index", i < N_COL);
      return columns()[i];
    }

    /// The i-th column
    inline restrict_column_t restrict_column(size_t i) {
      SL_REQUIRE("Good index", i < N_COL);
      return restrict_columns()[i];
    }

    /// Accessor to matrix rows
    inline row_accessor_t rows() {
      return row_accessor_t(begin());
    }

    /// Accessor to matrix rows
    inline const_row_accessor_t rows() const {
      return const_row_accessor_t(begin());
    }

    /// Accessor to matrix rows
    inline restrict_row_accessor_t restrict_rows() {
      return restrict_row_accessor_t(restrict_begin());
    }

    /// The i-th row
    inline row_t row(size_t i) {
      SL_REQUIRE("Good index", i < N_ROW);
      return rows()[i];
    }

    /// The i-th row
    inline const_row_t row(size_t i) const {
      SL_REQUIRE("Good index", i < N_ROW);
      return rows()[i];
    }

    /// The i-th row
    inline restrict_row_t restrict_row(size_t i) {
      SL_REQUIRE("Good index", i < N_ROW);
      return restrict_rows()[i];
    }

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

    /// set all diagonal components to s
    void fill_diagonal(value_t s) { 
      // Reimplement with iterators!
      size_t n = min(row_count(),column_count());
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
    template <class SELF_T2, class TRANSPOSED_T2, class VALUE_T2>
    inline self_t& operator+= (const fixed_size_matrix_base<SELF_T2,TRANSPOSED_T2, N_ROW, N_COL, VALUE_T2>& other) {
#if SL_NO_TEMPLATE_METAPROGRAMS 
      const typename fixed_size_matrix_base<SELF_T2,TRANSPOSED_T2, N_ROW, N_COL, VALUE_T2>::const_iterator other_it = other.begin();
      for (size_t i=0; i<element_size; ++i) {
	typedef SL_PROMOTENAME(value_t,VALUE_T2) arg_t;
	storage_[i] = value_t(arg_t(storage_[i]) + arg_t(other_it[i]));
      }
#else
     fastest::self_transform<element_size>::apply(other.begin(), restrict_begin(), plus_op<value_t,VALUE_T2,value_t>());
#endif
      return static_cast<self_t&>(*this);
    }

    /// subtract other from this
    template <class SELF_T2, class TRANSPOSED_T2, class VALUE_T2>
    inline self_t& operator-= (const fixed_size_matrix_base<SELF_T2,TRANSPOSED_T2, N_ROW, N_COL, VALUE_T2>& other) {
#if SL_NO_TEMPLATE_METAPROGRAMS 
      const typename fixed_size_matrix_base<SELF_T2,TRANSPOSED_T2, N_ROW, N_COL, VALUE_T2>::const_iterator other_it = other.begin();
      for (size_t i=0; i<element_size; ++i) {
	typedef SL_PROMOTENAME(value_t,VALUE_T2) arg_t;
	storage_[i] = value_t(arg_t(storage_[i]) - arg_t(other_it[i]));
      }
#else
      fastest::self_transform<element_size>::apply(other.begin(), restrict_begin(), rev_minus_op<value_t,VALUE_T2,value_t>());
#endif
      return static_cast<self_t&>(*this);
    }

    /// change sign to all components
    inline self_t operator- () const {
      self_t result = tags::not_initialized();
#if SL_NO_TEMPLATE_METAPROGRAMS 
      restrict_iterator result_it = result.restrict_begin();
      for (size_t i=0; i<element_size; ++i) {
	result_it[i] = - storage_[i];
      }
#else
      fastest::transform<element_size>::apply(begin(), result.restrict_begin(), negate_op<value_t,value_t>());
#endif
      return result;
    }
    
    /// multiply this by scalae
    template <class VALUE_T2>
    inline self_t& operator*= (VALUE_T2 scalar) {
#if SL_NO_TEMPLATE_METAPROGRAMS 
      for (size_t i=0; i<element_size; ++i) {
	typedef SL_PROMOTENAME(value_t,VALUE_T2) arg_t;
	storage_[i] = value_t(arg_t(storage_[i]) * arg_t(scalar));
      }
#else
      fastest::self_transform<element_size>::apply(restrict_begin(), scale_op<value_t,VALUE_T2,value_t>(scalar));
#endif
      return static_cast<self_t&>(*this);
    }
      
    /// divide this by scalar
    template <class VALUE_T2>
    inline self_t & operator/= (VALUE_T2 scalar) {
      return *this *= reciprocal(scalar);
    }

    SL_OP_LINEAR_SPACE(self_t,value_t); // Extend to generic types!

    /// the transposed of this
    transposed_t transposed() const {
      transposed_t result = tags::not_initialized();

      for(size_t i=0; i < N_ROW; i++) {
        for(size_t j=0; j < N_COL; j++) {
	  result(j,i) = (*this)(i,j);
	}           
      }

      return result;
    }

  public: // Interpolation

    /// Linear interpolation
    template <class T_PARAMETER>
    self_t lerp(const self_t& other, T_PARAMETER t) const {
      self_t result =  tags::not_initialized();
      fastest::transform<element_size>::apply(begin(), other.begin(), result.restrict_begin(), lerp_op<value_t,T_PARAMETER>(t));
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

  public: // Norm

    /// the square of the Euclidean norm of the matrix (sum of squares)
    SL_SUMTYPENAME(value_t) two_norm_squared() const {
      typedef SL_SUMTYPENAME(value_t)          T_sumtype;
      T_sumtype result = fastest::accumulate<element_size>::apply(begin(), scalar_math<SL_SUMTYPENAME(value_t)>::zero(), add_sqr_op<value_t,value_t,SL_SUMTYPENAME(value_t)>());
      return result;
    }

    /// the Euclidean norm of the matrix (sum of squares)
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) two_norm() const {
      typedef SL_SUMTYPENAME(value_t)          T_sumtype;
      typedef SL_FLOATTYPENAME(T_sumtype)      T_floattype;
      T_floattype result = std::sqrt(T_floattype(two_norm_squared()));
      return result;
    }

    /// the Euclidean norm of the matrix (sum of squares)
    inline SL_FLOATTYPENAME(SL_SUMTYPENAME(value_t)) frobenius_norm() const {
      return two_norm();
    }

    /// the maximum of column abs sum
    SL_SUMTYPENAME(value_t) one_norm() const {
      // rewrite with iterators!
      typedef SL_SUMTYPENAME(value_t) sum_t;
      sum_t result = scalar_math<sum_t>::zero();
      for (size_t j = 0; j < column_count(); j++) {
	sum_t asum = scalar_math<sum_t>::zero();
	for (size_t i = 0; i< row_count(); i++) {
	  asum += sl::abs(sum_t((*this)(i,j)));
	}
	result = sl::max(result, asum);
      }
      return result;
    }

    /// the maximum of row abs sum
    SL_SUMTYPENAME(value_t) infinite_norm() const {
      // rewrite with iterators!
      typedef SL_SUMTYPENAME(value_t) sum_t;
      sum_t result = scalar_math<sum_t>::zero();
      for (size_t i = 0; i< row_count(); i++) {
	sum_t asum = scalar_math<sum_t>::zero();
	for (size_t j = 0; j < column_count(); j++) {
	  asum += sl::abs(sum_t((*this)(i,j)));
	}
	result = sl::max(result, asum);
      }
      return result;
    }

  }; // class fixed_size_matrix_base

}; // namespace sl

// I/O

template <class SELF_T, class TRANSPOSED_T, size_t N_ROW, size_t N_COL, class T>
std::ostream& operator<<(std::ostream& s, const sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,N_ROW,N_COL,T>& a) {
  for (size_t i = 0; i<a.row_count(); i++) {
    for (size_t j = 0; j<a.column_count(); j++) {
      s << a(i,j) << " ";
    }
    s << std::endl;
  }
  return s;
}
    
template <class SELF_T, class TRANSPOSED_T, size_t N_ROW, size_t N_COL, class T>
std::istream& operator>>(std::istream& s, sl::fixed_size_matrix_base<SELF_T,TRANSPOSED_T,N_ROW,N_COL,T>& a) {
  for (size_t i = 0; i<a.row_count(); i++) {
    for (size_t j = 0; j<a.column_count(); j++) {
      s >> a(i,j) ;
    }
  }
  return s;
}

namespace sl {

  /// Rectangular matrixes of fixed size
  template <size_t N_ROW, size_t N_COL, class T>
  class fixed_size_matrix: 
    public fixed_size_matrix_base< fixed_size_matrix<N_ROW, N_COL, T>, 
                              fixed_size_matrix<N_COL, N_ROW, T>,
                              N_ROW, N_COL, T > {

  public: // Types

    typedef fixed_size_matrix_base< fixed_size_matrix<N_ROW, N_COL, T>, 
                               fixed_size_matrix<N_COL, N_ROW, T>,
                               N_ROW, N_COL, T > super_t;

    typedef fixed_size_matrix<N_ROW, N_COL, T>  this_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_matrix() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    inline fixed_size_matrix(const super_t& other): super_t(other) {
    }

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as m = 1.0; (fill) and 
     *  m = 1.0, 2.0, 3.0, 4.0; (component init)
     */
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, N_ROW, N_COL, v);
    }

  public: // Multiplication

    /// the multiplication of this with other
    template <size_t OTHER_N_COL>
    fixed_size_matrix<N_ROW, OTHER_N_COL,value_t> operator*(const fixed_size_matrix<N_COL, OTHER_N_COL,value_t>& other) const {
      fixed_size_matrix<N_ROW, OTHER_N_COL,value_t> result = tags::not_initialized();
      // TODO: REORDER i,k,j TO ENABLE VECTORIZATION
      for(size_t i=0; i < OTHER_N_COL; i++) {
        for(size_t j=0; j < N_ROW; j++) {
	  value_t sum = (*this)(j,0) * other(0,i);
	  for(size_t k=1; k < N_COL; k++) {
	    sum += (*this)(j,k) * other(k,i);
	  }
	  result(j,i) = sum;
	}           
      }

      return result;
    }

#if 0

  public: // Inversion

    /// Pseudo inverse implementation -- WARNING: Really quick'n'dirty
    void pseudo_invert_to(sl::fixed_size_matrix<N_COL,N_ROW,T>& Ainv,
			  bool *ok);

    /// The inverse of other (requires pseudo-invertible!)
    sl::fixed_size_matrix<N_COL,N_ROW,T> pseudo_inverse() const {
      bool ok;
      self_t result = tags::not_initialized();
      pseudo_invert_to(result, &ok);
      SL_CHECK("Invertible", ok);
      return result;
    }
#endif

  }; // class fixed_size_matrix


  template <size_t N_ROW, size_t N_COL, class OUT_ET>
  class conv_to< fixed_size_matrix<N_ROW, N_COL, OUT_ET> > {
  public:
    typedef fixed_size_matrix<N_ROW, N_COL, OUT_ET> result_t;

    // Explicit conversion from matrices of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_matrix_base< fixed_size_matrix<N_ROW, N_COL, IN_ET>, 
							      fixed_size_matrix<N_COL, N_ROW, IN_ET>,
							      N_ROW, N_COL, IN_ET >& in) {
      result_t result = tags::not_initialized();
      for(size_t i=0; i < N_ROW; i++) {
        for(size_t j=0; j < N_COL; j++) {
	  result(i,j) = static_cast<OUT_ET>(in(i,j));
	}
      }
      return result;
    }
      
  }; // class conv_to

} // namespace sl

// Arithmetic operators overloads
template <size_t N_ROW, size_t N_COL, class T>
inline sl::fixed_size_matrix<N_ROW,N_COL,T> operator*(const T& y, const sl::fixed_size_matrix<N_ROW,N_COL,T>& x) {
  typedef sl::fixed_size_matrix<N_ROW,N_COL,T> self_t;

  self_t result = sl::tags::not_initialized();
   
  const typename self_t::restrict_iterator result_it = result.restrict_begin();
  const typename self_t::const_iterator    src_it    = x.begin();
  for (size_t i=0; i<self_t::element_size; ++i) {
    result_it[i] = src_it[i] * y;
  }
  return result;
}

#if 0

#include <sl/fixed_size_square_matrix.hpp>

// Quick and dirty pseudo inverse implementation
template <size_t N_ROW, size_t N_COL, class T>
void sl::fixed_size_matrix<N_ROW,N_COL,T>::pseudo_invert_to(sl::fixed_size_matrix<N_COL,N_ROW,T>& Ainv,
							    bool *ok) {
  assert(ok != 0);
  
  sl::fixed_size_matrix<N_COL,N_ROW,T> AT  = this->transposed();
  sl::fixed_size_matrix<N_COL,N_COL,T> ATA = AT*(*this);
  sl::fixed_size_square_matrix<N_COL,T> ATAs;
  for (std::size_t i=0; i<N_COL; ++i) {
    for (std::size_t j=0; j<N_COL; ++j) {
      ATAs(i,j) = ATA(i,j);
    }
  }
  sl::fixed_size_square_matrix<N_COL,T> ATAinvs; 
  ATAs.invert_to(ATAinvs,ok);
  if (*ok) {
    sl::fixed_size_matrix<N_COL,N_COL,T> ATAinv;
    for (std::size_t i=0; i<N_COL; ++i) {
      for (std::size_t j=0; j<N_COL; ++j) {
	ATAinv(i,j) = ATAinvs(i,j);
      }
    }
    
    Ainv = ATAinv * AT;
  }

}
#endif


#endif 






