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
#ifndef SL_FIXED_SIZE_SQUARE_MATRIX_HPP
#define SL_FIXED_SIZE_SQUARE_MATRIX_HPP

#include <sl/conv_to.hpp>
#include <sl/fixed_size_matrix.hpp>
#include <sl/fixed_size_vector.hpp>

namespace sl {

  /// Base class for fixed size square matrices
  template <class SELF_T, size_t DIMENSION, class T>
  class fixed_size_square_matrix_base: 
    public fixed_size_matrix_base<SELF_T,SELF_T,DIMENSION,DIMENSION,T> {
  public: // Constants and types
    
    typedef fixed_size_matrix_base<SELF_T,SELF_T,DIMENSION,DIMENSION,T> super_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    typedef fixed_size_vector<column_orientation, DIMENSION,T > vector_t;
    typedef fixed_size_vector<row_orientation, DIMENSION,T > dual_vector_t;
    
    enum { dimension = DIMENSION };

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Good dimension", DIMENSION > 0);

  protected: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_square_matrix_base() {
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_square_matrix_base(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_square_matrix_base(const super_t& other): super_t(other) {
    }

    /**
     * Initialize from manifest constant. This allows initializations such
     * as m = 1.0 (fill) and m = 1.0, 2.0, 3.0, 4.0. 
     */ 
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, dimension, dimension, v);
    }

  public: // Initializations

    /// Set this to the identity matrix
    void to_identity() {
      this->to_zero();
      this->fill_diagonal(scalar_math<value_t>::one());
    }

  public: // Queries

    /// Is this an homogeneous matrix representing an affine transform?
    bool is_affine() const {
      for (size_t i=0; i<dimension-1; i++) {
	if ((*this)(dimension-1,i) != scalar_math<value_t>::zero()) return false;
      }
      return (*this)(dimension-1,dimension-1) == scalar_math<value_t>::one();
    }

    /// Is this an identity matrix?
    bool is_identity() const {
      for (size_t i = 0; i<dimension; i++) {
	for (size_t j = 0; j< dimension; j++) {
	  if ((*this)(i,j) != ((i==j) ? scalar_math<value_t>::one() : scalar_math<value_t>::zero())) return false;
	}
      }
      return true;
    }

    /// Is this an homogeneous matrix representing a scale-preserving transform?
    inline bool is_scale_preserving() const {
      return 
        sl::abs(determinant() - scalar_math<value_t>::one()) < 
	static_cast<value_t>(100 * dimension * dimension * dimension) * std::numeric_limits<value_t>::epsilon(); 
    }
    
    /// Is the matrix a right handed basis?
    inline bool is_right_handed() const {
      return sl::is_positive(determinant());
    }

    /// Is the matrix representing a left handed basis?
    inline bool is_left_handed() const {
      return sl::is_negative(determinant());
    }

  public: // Handedness change

    /// Change the handedness of the matrix
    inline void toggle_handedness() {
      // Invert direction of first basis vector 
      for (size_t i=0; i< dimension; i++) {
	(*this)(i,0) = -(*this)(i,0);
      }
    }

    /// Make the matrix right handed (if not singular)
    inline void make_right_handed() {
      SL_REQUIRE("Not singular", is_right_handed() || is_left_handed());
      if (!is_right_handed()) {
	toggle_handedness();
      }
      SL_ENSURE("is_right_handed", is_right_handed());
    }

    /// Make the matrix left handed (if not singular)
    inline void make_left_handed() {
      SL_REQUIRE("Not singular", is_right_handed() || is_left_handed());
      if (!is_left_handed()) {
	toggle_handedness();
      }
      SL_ENSURE("is_right_handed", is_left_handed());
    }

  public: // Division algebra operations

    /// this multiplied by other
    self_t operator *(const self_t& other) const {
      self_t result = tags::not_initialized();

      // General multiply
      for(int i=0; i < dimension; i++) {
	for(int j=0; j < dimension; j++) {
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
      value_t det;
      general_invert(static_cast<const self_t&>(*this), mo, &det);
      *ok = !is_zero(det);
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

  public: // G++ Workarounds

    /// scale this by scalar
    inline self_t& operator *= (const value_t& scalar) { 
      return super_t::operator *= (scalar); 
    }

    /// scale this by the inverse of scalar
    inline self_t& operator /= (const value_t& scalar) { 
      return super_t::operator /= (scalar); 
    }

  public: // Transposition

    /// transpose all components
    inline void transpose() {
      *this = transposed(); // REWRITE!
    }                                                                                                                          
    
    /// the transposed matrix
    self_t transposed() const {
      self_t result = tags::not_initialized();

      for(int i=0; i < dimension; i++) {
        for(int j=0; j < dimension; j++) {
	  result(j,i) = (*this)(i,j);
	}           
      }

      return result;
    }

  public: // Determinant

    /// the determinant of this (currently very slow!)
    inline value_t determinant() const {
      // HACK, SHOULD OPTIMIZE...
      value_t result;
      self_t tmp = tags::not_initialized();
      general_invert(static_cast<const self_t&>(*this), tmp, &result);
      return result;
    }

  public: // Eigenvalues / eigenvectors

    /** 
     * Computes the eigenvectors and eigenvalues of symmetric square matrices
     * by the Jacobi method. this must be a symmetric matrix, d will
     * be filled with its eigenvalues and the columns of v will be
     * the eigenvectors.
     * Based on Numerical Recipes.
     */
    void symmetric_eigen_in(vector_t& d, self_t& v, bool *ok) const {
#define SL_SE_ROT(a, i, j, k, l) { value_t g = a(i,j); value_t h = a(k,l); a(i,j) = g-s*(h+g*tau); a(k,l) = h+s*(g-h*tau); }

      SL_REQUIRE("OK exists", ok != NULL);
      SL_REQUIRE("Symmetric", true);

      const std::size_t MAX_ITER = 50;
      const value_t Zero = scalar_math<value_t>::zero();
      const value_t One  = scalar_math<value_t>::one();

      *ok = false;
      self_t a = (static_cast<const self_t&>(*this));

      v.to_identity();
      
      vector_t b, z;
      for (std::size_t ip=0; ip<dimension; ++ip) {
	b[ip] = d[ip] = a(ip,ip);
	z[ip] = Zero;
      }

      for (std::size_t i=0; i<MAX_ITER && !(*ok); i++) {
        // Sum off-diagonal elements
	value_t sm = Zero;
        for (std::size_t ip=0; ip<dimension-1; ++ip) {
          for (std::size_t iq=ip+1; iq<dimension; ++iq) {
            sm += sl::abs(a(ip,iq));
          }
	}
        // Convergence to machine precision
	// This assumes undeflows -> zero
	*ok = (sm == Zero);
	if (!*ok) {
	  const value_t tresh = (i<4) ? static_cast<value_t>(0.2/sqr(double(dimension)))*sm : Zero;

	  for (std::size_t ip=0; ip<dimension-1; ++ip) {
            for (std::size_t iq =ip+1; iq<dimension; ++iq) {
	      // After four sweeps, skip rotation if the off-diagonal element is small
              value_t g = static_cast<value_t>(100)*sl::abs(a(ip,iq));
              if ((i>4) &&
		  (float(sl::abs(d[ip])+g) == float(sl::abs(d[ip]))) &&
		  (float(sl::abs(d[iq])+g) == float(sl::abs(d[iq])))) {
                a(ip,iq)=Zero;
              } else if (sl::abs(a(ip,iq))>tresh) {
                value_t h =d[iq]-d[ip];
		value_t t;
                if (float(sl::abs(h)+g) == float(sl::abs(h))) {
                  t = a(ip,iq)/h;
                } else {
                  value_t theta = static_cast<value_t>(0.5)*h/a(ip,iq);
                  t = One/(sl::abs(theta)+std::sqrt(One+sqr(theta)));
                  if (theta<Zero) { t = -t; }
                }
                value_t c = One/std::sqrt(One+sqr(t)); 
		value_t s = t*c; 
		value_t tau = s/(One+c);
                h =t*a(ip,iq);
		z[ip] -= h; z[iq] += h;
		d[ip] -= h; d[iq] += h;
		a(ip,iq)=Zero;

                for (std::size_t j=0; j+1<=ip; ++j) {
		  SL_SE_ROT(a, j, ip, j, iq);
                } 
                for (std::size_t j=ip+1; j+1<=iq; ++j) {
		  SL_SE_ROT(a, ip, j, j, iq);
                }
                for (std::size_t j=iq+1; j<dimension; j++) {
		  SL_SE_ROT(a, ip, j, iq, j);
		}
                for (std::size_t j =0; j<dimension; j++) {
		  SL_SE_ROT(v, j, ip, j, iq);
                }
              }
            }
          }
          b += z; d = b; z.to_zero();
	}
      }

#undef SL_SE_ROT
    }
    
    /** 
     * Computes the eigenvectors and eigenvalues of symmetric square matrices
     * by the Jacobi method. this must be a symmetric matrix, d will
     * be filled with its eigenvalues and the columns of v will be
     * the eigenvectors. Eigenvalues and eigenvectors are sorted by
     * decreasing eigenvalue magnitude. 
     * Based on Numerical Recipes.
     */
    void symmetric_sorted_eigen_in(vector_t& d, self_t& v, bool *ok) const {
      SL_REQUIRE("OK exists", ok != NULL);
      SL_REQUIRE("Symmetric", true);\
      symmetric_eigen_in(d,v,ok);
      if (*ok) {
	// Sort by eigenvalue magnitude
	for (size_t i = 0; i<dimension; i++) {
	  size_t k = i; value_t p = d[k];
	  for (size_t j=i+1; j<dimension; j++) {
	    if (d[j]>p) {
	      p = d[j];
	      k = j;
	    }
	  }
	  if (k!= i) {
	    d[k] = d[i]; d[i] = p;
	    for (size_t j=0; j<dimension; j++) {
	      p = v(j,i); 
	      v(j,i) = v(j,k);
	      v(j,k) = p;
	    }
	  }
	}
      }
    }

    /// The i-th axis of the basis
    inline vector_t axis(size_t i) const {
      SL_REQUIRE("Good index", i<dimension);

      vector_t result = tags::not_initialized();
      for (size_t j=0; j<dimension; ++j) {
	result[j] = (*this)(i,j);
      }
      return result;
    }

  public: // sqrt

    /// The square root of the symmetric matrix A
    void symmetric_sqrt_in(self_t& A_sqrt, bool *ok) const {
      self_t V = tags::not_initialized();
      vector_t d = tags::not_initialized();
      symmetric_eigen_in(d, V, ok);
      if (!(*ok)) { std::cerr << "SVD FAILED" << std::endl; }
      if (*ok) {
        self_t D = tags::not_initialized();
        for (size_t i=0; i<dimension; ++i) {
          D(i,i) = sl::is_positive(d[i]) ? std::sqrt(d[i]) : sl::scalar_math<T>::zero();
          *ok = *ok && !(d[i]<=0.0);
          if (!*ok) { std::cerr << "NEGATIVE d" << std::endl << d << std::endl; }
        }
        A_sqrt = V * D * (V.transposed());
      }
    }

    // The rotational part of Apq=this, defined as R = Apq * sqrt(Apq^T*Apq);
    void rotational_part_in(self_t& R, bool *ok) const {
      self_t A = transposed()*(*this);
      self_t S = tags::not_initialized();
      A.symmetric_sqrt_in(S, ok);
      if (*ok) {
        self_t S_inv = tags::not_initialized();
        S.invert_to(S_inv, ok);
        //if (!*ok) { std::cerr << "SINGULAR S" << std::endl << S << std::endl;  }
        R = (*this) * S_inv;
      }
    }

  public: // rank 1-2 updat

    /// this += alpha*v*v^T
    void rank_one_update(const vector_t& v,
                         const value_t& alpha = sl::scalar_math<value_t>::one()) {
      for (std::size_t i=0; i<dimension; ++i) {
        for (std::size_t j=0; j<dimension; ++j) {
          (*this)(i,j) += alpha * v[i] * v[j];
        }
      }
    }

    /// this += alpha*v*w^T
    void rank_two_update(const vector_t& v,
                         const vector_t& w,
                         const value_t& alpha = sl::scalar_math<value_t>::one()) {
      for (std::size_t i=0; i<dimension; ++i) {
        for (std::size_t j=0; j<dimension; ++j) {
          (*this)(i,j) += alpha * v[i] * w[j];
        }
      }
    }
    
    void to_orthogonal_projection_onto_line(const vector_t& v) {
      this->to_zero();
      this->rank_two_update(v, v, sl::reciprocal(v.dot(v))); 
    }

    static self_t orthogonal_projection_onto_line(const vector_t& v) {
      self_t result = sl::tags::not_initialized();
      result.to_orthogonal_projection_onto_line(v);
      return result;
    }

  protected: // Helpers

    /// General square matrix inversion
    static void general_invert(const self_t& mat, 
			       self_t& inv,
			       value_t* determ) {

      inv = mat;

      int i,j,k;
      int pvt_i[dimension], pvt_j[dimension]; // Locations of pivot elements 
      value_t pvt_val;                     // Value of current pivot element
      value_t hold;                        // Temporary storage

      *determ = scalar_math<value_t>::one();

      for (k = 0; k < (int)dimension; k++) {
        // Locate k'th pivot element
        pvt_val = inv(k,k);      
        pvt_i[k] = k;
        pvt_j[k] = k;
        for (i = k; i < (int)dimension; i++) {
          for (j = k; j < (int)dimension; j++) {
            if (sl::abs(inv(i,j)) > sl::abs(pvt_val)) {
	      pvt_i[k] = i;
	      pvt_j[k] = j;
	      pvt_val = inv(i,j);
            }
	  }
	}

        // Product of pivots, gives determinant when finished 
        *determ *= pvt_val;
        if (is_zero(*determ)) return;

	// "Interchange" rows (with sign change stuff)
        i = pvt_i[k];
        if (i != k)  { 
          for (j = 0; j < (int)dimension; j++) {
            hold = -inv(k,j);
            inv(k,j) = inv(i,j);
            inv(i,j) = hold;
          }
	}

        // "Interchange" columns
        j = pvt_j[k];
        if (j != k) {
          for (i = 0; i < (int)dimension; i++) {
            hold = -inv(i,k);
            inv(i,k) = inv(i,j);
            inv(i,j) = hold;
          }
	}
	
        // Divide column by minus pivot value 
        for (i = 0; i < (int)dimension; i++) {
          if (i != k)                   // Don't touch the pivot entry 
            inv(i,k) /= (-pvt_val) ;  
	}


        // Reduce the invrix 
        for (i = 0; i < (int)dimension; i++) {
	  hold = inv(i,k);
	  for (j = 0; j < (int)dimension; j++) {
	    if ( i != k && j != k )     
	      inv(i,j) += hold * inv(k,j);
	  }
        }

        // Divide row by pivot
        for (j = 0; j < (int)dimension; j++) {
          if (j != k)                   
            inv(k,j) /= pvt_val;
	}

        inv(k,k) = scalar_math<value_t>::one();
	inv(k,k) /= pvt_val;
      }

      // That was most of the work, one final pass of row/column interchange 
      // to finish
      for (k = (int)dimension-2; k >= 0; k--) { 
        // Rows to swap correspond to pivot COLUMN 
	i = pvt_j[k];            
        if (i != k) { 
          for(j = 0; j < (int)dimension; j++) {
            hold = inv(k,j);
            inv(k,j) = -inv(i,j);
            inv(i,j) = hold;
          }
	}
	// Columns to swap correspond to pivot ROW */
        j = pvt_i[k];           
        if (j != k) {                    
          for (i = 0; i < (int)dimension; i++) {
            hold = inv(i,k);
            inv(i,k) = -inv(i,j);
            inv(i,j) = hold;
          }
	}
      }
    }

  }; // fixed_size_square_matrix_base

}; // namespace sl

namespace sl {

  /// Fixed size square matrices
  template <size_t DIMENSION, class T>
  class fixed_size_square_matrix: 
    public fixed_size_square_matrix_base< fixed_size_square_matrix<DIMENSION, T>, 
                                     DIMENSION, T > {

  public: // Types

    typedef fixed_size_square_matrix_base< fixed_size_square_matrix<DIMENSION, T>, 
                                     DIMENSION, T > super_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    enum { dimension = DIMENSION };

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_square_matrix() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_square_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_square_matrix(const super_t& other): super_t(other) {
    }

    /**
     * Initialize from manifest constant. This allows initializations such
     * as m = 1.0 (fill) and m = 1.0, 2.0, 3.0, 4.0. 
     */ 
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, dimension, dimension, v);
    }

  };

  template <size_t DIMENSION, class OUT_ET>
  class conv_to< fixed_size_square_matrix<DIMENSION, OUT_ET> > {
  public:
    typedef fixed_size_square_matrix<DIMENSION, OUT_ET> result_t;

    // Explicit conversion from matrices of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_matrix<DIMENSION, DIMENSION, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for(size_t i=0; i < DIMENSION; i++) {
        for(size_t j=0; j < DIMENSION; j++) {
	  result(i,j) = static_cast<OUT_ET>(in(i,j));
	}
      }
      return result;
    }
      
    // Explicit conversion from matrices of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_square_matrix<DIMENSION, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for(size_t i=0; i < DIMENSION; i++) {
        for(size_t j=0; j < DIMENSION; j++) {
	  result(i,j) = static_cast<OUT_ET>(in(i,j));
	}
      }
      return result;
    }
      
  }; // class conv_to

}; // namespace sl

// Arithmetic operators overloads

template <size_t DIMENSION, class T>
inline sl::fixed_size_square_matrix<DIMENSION,T> operator*(const T& y, const sl::fixed_size_square_matrix<DIMENSION,T>& x) { 
  typedef sl::fixed_size_square_matrix<DIMENSION,T> self_t;

  self_t result = sl::tags::not_initialized();
  
  const typename self_t::restrict_iterator result_it = result.restrict_begin();
  const typename self_t::const_iterator    src_it    = x.begin();
  for (size_t i=0; i<self_t::element_size; ++i) {
    result_it[i] = src_it[i] * y;
  }
  return result;
}

namespace sl {

  /// 4x4 matrices
  template <class T>
  class fixed_size_square_matrix<4,T>: 
    public fixed_size_square_matrix_base< fixed_size_square_matrix<4, T>, 
                                     4, T > {

  public: // Types
    typedef fixed_size_square_matrix_base< fixed_size_square_matrix<4, T>, 
                                      4, T > super_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;

    enum { dimension = super_t::dimension };

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_square_matrix() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_square_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_square_matrix(const super_t& other): super_t(other) {
    }

    /// Component by component init
    explicit inline fixed_size_square_matrix(const value_t& m00, const value_t& m01, const value_t& m02, const value_t& m03,
					     const value_t& m10, const value_t& m11, const value_t& m12, const value_t& m13,
					     const value_t& m20, const value_t& m21, const value_t& m22, const value_t& m23,
					     const value_t& m30, const value_t& m31, const value_t& m32, const value_t& m33)  
      : 
      super_t(tags::not_initialized()) {
      (*this)(0,0) = m00; (*this)(0,1) = m01; (*this)(0,2) = m02; (*this)(0,3) = m03;
      (*this)(1,0) = m10; (*this)(1,1) = m11; (*this)(1,2) = m12; (*this)(1,3) = m13;
      (*this)(2,0) = m20; (*this)(2,1) = m21; (*this)(2,2) = m22; (*this)(2,3) = m23;
      (*this)(3,0) = m30; (*this)(3,1) = m31; (*this)(3,2) = m32; (*this)(3,3) = m33;
    }

    /**
     * Initialize from manifest constant. This allows initializations such
     * as m = 1.0 (fill) and m = 1.0, 2.0, 3.0, 4.0. 
     */ 
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, dimension, dimension, v);
    }

  protected: // Special cases

    static void affine4x4_invert(const self_t& mat, 
				 self_t& inv,
				 value_t* determ) {
      SL_REQUIRE("4x4", dimension == 4);

      inv(0,0) = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2);
      inv(0,1) = mat(0,1)*mat(2,2) - mat(2,1)*mat(0,2);
      inv(0,2) = mat(0,1)*mat(1,2) - mat(1,1)*mat(0,2);

      inv(1,0) = mat(1,0)*mat(2,2) - mat(2,0)*mat(1,2);
      inv(1,1) = mat(0,0)*mat(2,2) - mat(2,0)*mat(0,2);
      inv(1,2) = mat(0,0)*mat(1,2) - mat(1,0)*mat(0,2);

      inv(2,0) = mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1);
      inv(2,1) = mat(0,0)*mat(2,1) - mat(2,0)*mat(0,1);
      inv(2,2) = mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);

      inv(3,0) = 0;
      inv(3,1) = 0;
      inv(3,2) = 0;
      inv(3,3) = 1;

      *determ = mat(0,0)*inv(0,0) - mat(1,0)*inv(0,1) + mat(2,0)*inv(0,2);

      if (is_zero(*determ)) return;

      inv(0,0) /= *determ;
      inv(2,0) /= *determ;
      inv(1,1) /= *determ;
      inv(0,2) /= *determ;
      inv(2,2) /= *determ;

      *determ *= -1.;

      inv(1,0) /= *determ;
      inv(0,1) /= *determ;
      inv(2,1) /= *determ;
      inv(1,2) /= *determ;
 
      inv(0,3) = -(inv(0,0)*mat(0,3)+inv(0,1)*mat(1,3)+inv(0,2)*mat(2,3));
      inv(1,3) = -(inv(1,0)*mat(0,3)+inv(1,1)*mat(1,3)+inv(1,2)*mat(2,3));
      inv(2,3) = -(inv(2,0)*mat(0,3)+inv(2,1)*mat(1,3)+inv(2,2)*mat(2,3));
    }

    static self_t general4x4_mult_general4x4(const self_t& m1, const self_t& m2) {
      SL_REQUIRE("4x4", dimension == 4);
      return 
	self_t(
	     m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0)+m1(0,2)*m2(2,0)+m1(0,3)*m2(3,0),
	     m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1)+m1(0,2)*m2(2,1)+m1(0,3)*m2(3,1),
	     m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2)+m1(0,2)*m2(2,2)+m1(0,3)*m2(3,2),
	     m1(0,0)*m2(0,3)+m1(0,1)*m2(1,3)+m1(0,2)*m2(2,3)+m1(0,3)*m2(3,3),
	     
	     m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0)+m1(1,2)*m2(2,0)+m1(1,3)*m2(3,0),
	     m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1)+m1(1,2)*m2(2,1)+m1(1,3)*m2(3,1),
	     m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)+m1(1,2)*m2(2,2)+m1(1,3)*m2(3,2),
	     m1(1,0)*m2(0,3)+m1(1,1)*m2(1,3)+m1(1,2)*m2(2,3)+m1(1,3)*m2(3,3),
	     
	     m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0)+m1(2,2)*m2(2,0)+m1(2,3)*m2(3,0),
	     m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1)+m1(2,2)*m2(2,1)+m1(2,3)*m2(3,1),
	     m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)+m1(2,2)*m2(2,2)+m1(2,3)*m2(3,2),
	     m1(2,0)*m2(0,3)+m1(2,1)*m2(1,3)+m1(2,2)*m2(2,3)+m1(2,3)*m2(3,3),
	     
	     m1(3,0)*m2(0,0)+m1(3,1)*m2(1,0)+m1(3,2)*m2(2,0)+m1(3,3)*m2(3,0),
	     m1(3,0)*m2(0,1)+m1(3,1)*m2(1,1)+m1(3,2)*m2(2,1)+m1(3,3)*m2(3,1),
	     m1(3,0)*m2(0,2)+m1(3,1)*m2(1,2)+m1(3,2)*m2(2,2)+m1(3,3)*m2(3,2),
	     m1(3,0)*m2(0,3)+m1(3,1)*m2(1,3)+m1(3,2)*m2(2,3)+m1(3,3)*m2(3,3)
	     );
    }
    
    static self_t general4x4_mult_affine4x4(const self_t& m1, const self_t& m2) {
      SL_REQUIRE("4x4", dimension == 4);
      SL_REQUIRE("m2 affine", m2.is_affine());
      return 
	self_t(
	     m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0)+m1(0,2)*m2(2,0),
	     m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1)+m1(0,2)*m2(2,1),
	     m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2)+m1(0,2)*m2(2,2),
	     m1(0,0)*m2(0,3)+m1(0,1)*m2(1,3)+m1(0,2)*m2(2,3)+m1(0,3),
	     
	     m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0)+m1(1,2)*m2(2,0),
	     m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1)+m1(1,2)*m2(2,1),
	     m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)+m1(1,2)*m2(2,2),
	     m1(1,0)*m2(0,3)+m1(1,1)*m2(1,3)+m1(1,2)*m2(2,3)+m1(1,3),
	     
	     m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0)+m1(2,2)*m2(2,0),
	     m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1)+m1(2,2)*m2(2,1),
	     m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)+m1(2,2)*m2(2,2),
	     m1(2,0)*m2(0,3)+m1(2,1)*m2(1,3)+m1(2,2)*m2(2,3)+m1(2,3),
	     
	     m1(3,0)*m2(0,0)+m1(3,1)*m2(1,0)+m1(3,2)*m2(2,0),
	     m1(3,0)*m2(0,1)+m1(3,1)*m2(1,1)+m1(3,2)*m2(2,1),
	     m1(3,0)*m2(0,2)+m1(3,1)*m2(1,2)+m1(3,2)*m2(2,2),
	     m1(3,0)*m2(0,3)+m1(3,1)*m2(1,3)+m1(3,2)*m2(2,3)+m1(3,3)
	     );
    }

    static self_t affine4x4_mult_general4x4(const self_t& m1, const self_t& m2) {
      SL_REQUIRE("4x4", dimension == 4);
      SL_REQUIRE("m1 affine", m1.is_affine());
      return 
	self_t(
	     m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0)+m1(0,2)*m2(2,0)+m1(0,3)*m2(3,0),
	     m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1)+m1(0,2)*m2(2,1)+m1(0,3)*m2(3,1),
	     m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2)+m1(0,2)*m2(2,2)+m1(0,3)*m2(3,2),
	     m1(0,0)*m2(0,3)+m1(0,1)*m2(1,3)+m1(0,2)*m2(2,3)+m1(0,3)*m2(3,3),
	     
	     m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0)+m1(1,2)*m2(2,0)+m1(1,3)*m2(3,0),
	     m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1)+m1(1,2)*m2(2,1)+m1(1,3)*m2(3,1),
	     m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)+m1(1,2)*m2(2,2)+m1(1,3)*m2(3,2),
	     m1(1,0)*m2(0,3)+m1(1,1)*m2(1,3)+m1(1,2)*m2(2,3)+m1(1,3)*m2(3,3),
	     
	     m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0)+m1(2,2)*m2(2,0)+m1(2,3)*m2(3,0),
	     m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1)+m1(2,2)*m2(2,1)+m1(2,3)*m2(3,1),
	     m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)+m1(2,2)*m2(2,2)+m1(2,3)*m2(3,2),
	     m1(2,0)*m2(0,3)+m1(2,1)*m2(1,3)+m1(2,2)*m2(2,3)+m1(2,3)*m2(3,3),
	     
	     m2(3,0),
	     m2(3,1),
	     m2(3,2),
	     m2(3,3)
	     );
    }

    static self_t affine4x4_mult_affine4x4(const self_t& m1, const self_t& m2) {
      SL_REQUIRE("4x4", dimension == 4);
      SL_REQUIRE("m1 affine", m1.is_affine());
      SL_REQUIRE("m2 affine", m2.is_affine());
      return 
	self_t(
	     m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0)+m1(0,2)*m2(2,0),
	     m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1)+m1(0,2)*m2(2,1),
	     m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2)+m1(0,2)*m2(2,2),
	     m1(0,0)*m2(0,3)+m1(0,1)*m2(1,3)+m1(0,2)*m2(2,3)+m1(0,3),
	     
	     m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0)+m1(1,2)*m2(2,0),
	     m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1)+m1(1,2)*m2(2,1),
	     m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)+m1(1,2)*m2(2,2),
	     m1(1,0)*m2(0,3)+m1(1,1)*m2(1,3)+m1(1,2)*m2(2,3)+m1(1,3),
	     
	     m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0)+m1(2,2)*m2(2,0),
	     m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1)+m1(2,2)*m2(2,1),
	     m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)+m1(2,2)*m2(2,2),
	     m1(2,0)*m2(0,3)+m1(2,1)*m2(1,3)+m1(2,2)*m2(2,3)+m1(2,3),
	     
	     sl::zero(value_t()),
	     sl::zero(value_t()),
	     sl::zero(value_t()),
	     sl::one(value_t())
	     );
    }

  public: // Optimizations 

    inline bool is_affine() const {
      return 
	is_zero((*this)(3,0)) &&
	is_zero((*this)(3,1)) &&
	is_zero((*this)(3,2)) &&
	is_one ((*this)(3,3));
    }

    inline self_t operator *(const self_t& other) const {
      if (is_affine()) {
	if (other.is_affine()) {
	  return affine4x4_mult_affine4x4(*this, other);
	} else {
	  return affine4x4_mult_general4x4(*this, other);
	}
      } else {
	if (other.is_affine()) {
	  return general4x4_mult_affine4x4(*this, other);
	} else {
	  return general4x4_mult_general4x4(*this, other);
	}
      }
    }

    inline void invert_to(self_t& mo, bool* ok) const {
      if (is_affine()) {
	value_t det;
	affine4x4_invert(*this, mo, &det);
	*ok = !is_zero(det);
      } else {
	value_t det;
	this->general_invert(*this, mo, &det);
	*ok = !is_zero(det);
      }
    }

  }; // class fixed_size_square_matrix<4,T>

}; // namespace sl  

namespace sl {

  template <class T>
  class fixed_size_square_matrix<3,T>: 
    public fixed_size_square_matrix_base< fixed_size_square_matrix<3, T>, 
                                     3, T > {

  public: // Types

    typedef fixed_size_square_matrix_base< fixed_size_square_matrix<3, T>, 
                                      3, T > super_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;


    enum { dimension = super_t::dimension };

  public: // Creation, Copy & Destruction

    
    /// Default init (zero)
    inline fixed_size_square_matrix() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_square_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }

    /// Set this to other
    inline fixed_size_square_matrix(const super_t& other): super_t(other) {
    }

    /// Component by component init
    inline explicit fixed_size_square_matrix(const value_t& m00, const value_t& m01, const value_t& m02,
					     const value_t& m10, const value_t& m11, const value_t& m12,
					     const value_t& m20, const value_t& m21, const value_t& m22)
      : 
      super_t(tags::not_initialized()) {
      (*this)(0,0) = m00; (*this)(0,1) = m01; (*this)(0,2) = m02;
      (*this)(1,0) = m10; (*this)(1,1) = m11; (*this)(1,2) = m12;
      (*this)(2,0) = m20; (*this)(2,1) = m21; (*this)(2,2) = m22;
    }
    
    /**
     * Initialize from manifest constant. This allows initializations such
     * as m = 1.0 (fill) and m = 1.0, 2.0, 3.0, 4.0. 
     */ 
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, dimension, dimension, v);
    }

  }; // class fixed_size_square_matrix<3,T>

}; // namespace sl

namespace sl {

  template <class T>
  class fixed_size_square_matrix<2,T>: 
    public fixed_size_square_matrix_base< fixed_size_square_matrix<2, T>, 
                                     2, T > {

  public: // Types

    typedef fixed_size_square_matrix_base< fixed_size_square_matrix<2, T>, 
                                      2, T > super_t;
    typedef typename super_t::self_t       self_t;
    typedef typename super_t::transposed_t transposed_t;
    typedef typename super_t::value_t      value_t;


    enum { dimension = super_t::dimension };

  public: // Creation, Copy & Destruction

    /// Default init (zero)
    inline fixed_size_square_matrix() {
      // Storage is already 0-filled
    }

    /// Fast init (garbage), handle with care!
    inline fixed_size_square_matrix(tags::not_initialized tag): super_t(tag) {
      // Garbage in, use with care
    }
    
    /// Set this to other
    inline fixed_size_square_matrix(const super_t& other): super_t(other) {
    }

    /// Component by component init
    explicit fixed_size_square_matrix(const value_t& m00, const value_t& m01,
				 const value_t& m10, const value_t& m11)    
      : 
      super_t(tags::not_initialized()) {

      SL_REQUIRE("Good dimension", dimension == 2);
      (*this)(0,0) = m00; (*this)(0,1) = m01;
      (*this)(1,0) = m10; (*this)(1,1) = m11;
    }
 
    /**
     * Initialize from manifest constant. This allows initializations such
     * as m = 1.0 (fill) and m = 1.0, 2.0, 3.0, 4.0. 
     */ 
    inline manifest_array2d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array2d_initializer<self_t&,value_t>(*this, dimension, dimension, v);
    }

  }; // class fixed_size_square_matrix<2,T>

}; // namespace sl

namespace sl {

  /// A factory for square matrices
  template <size_t DIMENSION, class T>
  class fixed_size_square_matrix_factory {
    
  public: // Constants and types 

    enum { dimension = DIMENSION };

    typedef fixed_size_square_matrix<DIMENSION, T> matrix_t;
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

  }; // class fixed_size_square_matrix_factory
      
}; // namespace sl

// ---------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------

namespace sl {

  /// A 2x2 matrix with single-precision floating point components
  typedef fixed_size_square_matrix<2,float> matrix2f;
  /// A 3x3 matrix with single-precision floating point components
  typedef fixed_size_square_matrix<3,float> matrix3f;
  /// A 4x4 matrix with single-precision floating point components
  typedef fixed_size_square_matrix<4,float> matrix4f;
  
  /// A 2x2 matrix factory with single-precision floating point components
  typedef fixed_size_square_matrix_factory<2,float> matrix_factory2f;
  /// A 3x3 matrix factory with single-precision floating point components
  typedef fixed_size_square_matrix_factory<3,float> matrix_factory3f;
  /// A 4x4 matrix factory with single-precision floating point components
  typedef fixed_size_square_matrix_factory<4,float> matrix_factory4f;

  /// A 2x2 matrix with double-precision floating point components
  typedef fixed_size_square_matrix<2,double> matrix2d;
  /// A 3x3 matrix with double-precision floating point components
  typedef fixed_size_square_matrix<3,double> matrix3d;
  /// A 4x4 matrix with double-precision floating point components
  typedef fixed_size_square_matrix<4,double> matrix4d;
  
  /// A 2x2 matrix factory with double-precision floating point components
  typedef fixed_size_square_matrix_factory<2,double> matrix_factory2d;
  /// A 3x3 matrix factory with double-precision floating point components
  typedef fixed_size_square_matrix_factory<3,double> matrix_factory3d;
  /// A 4x4 matrix factory with double-precision floating point components
  typedef fixed_size_square_matrix_factory<4,double> matrix_factory4d;
};


#endif




