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
#ifndef SL_FIXED_SIZE_ARRAY_HPP
#define SL_FIXED_SIZE_ARRAY_HPP

#include <sl/assert.hpp>
#include <sl/conv_to.hpp>
#include <sl/operators.hpp>
#include <sl/math.hpp>
#include <sl/interpolation.hpp>
#include <sl/fastest.hpp>
#include <sl/manifest_array_initializer.hpp>
#include <sl/serializer.hpp>
#include <sl/hash.hpp>
#include <functional>

namespace sl {

  namespace detail {

    template <size_t SZ, class PTR, bool INIT_REQUIRED> 
    class numeric_initializer {
    public:
    };

    template <size_t SZ, class PTR>
    class numeric_initializer<SZ, PTR, true> {
    public:
      static inline void apply(PTR ptr) {
        fastest::fill<SZ>::apply(ptr, sl::zero(*ptr));
      }
    };

    template <size_t SZ, class PTR>
    class numeric_initializer<SZ, PTR, false> {
    public:
      static inline void apply(PTR) {
        // No op
      }
    };
  }
 
  /// Arrays with size determined at compile time
  template <size_t SZ, class T>
  class fixed_size_array {
  public: // Constants and types

    enum { dimension = SZ };

    typedef fixed_size_array<dimension,T> self_t;
    typedef T                             value_t;

  public: // Type constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);
    //SL_COMPILE_TIME_CHECK("No extra baggage", sizeof(self_t) == dimension * sizeof(value_t));

  protected: // Implementation

    value_t storage_ [dimension];

  public: // STL-style iterators

    typedef value_t* iterator;
    typedef const value_t* const_iterator;
    typedef value_t* restrict restrict_iterator;

    static inline size_t size()                { return dimension; }
    inline iterator begin()                    { return &(storage_[0]); }
    inline iterator restrict_begin()           { return &(storage_[0]); }
    inline const_iterator begin() const        { return &(storage_[0]); }
    inline iterator end()                      { return &(storage_[0]) + dimension; }
    inline iterator restrict_end()             { return &(storage_[0]) + dimension; }
    inline const_iterator end() const          { return &(storage_[0]) + dimension; }

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      sl::fixed_size_vector_container_serialize_helper<iterator,const_iterator,sl::serialization_traits<value_t>::is_array_base_type>::store(s,begin(),end());
    }
    
    void retrieve_from(input_serializer& s) {
      sl::fixed_size_vector_container_serialize_helper<iterator,const_iterator,sl::serialization_traits<value_t>::is_array_base_type>::retrieve(s,begin(),end());
    }
    
  public: // Creation
    
    /// Default init (zero for numeric values, default constructor for others)
    inline fixed_size_array() {
      // Zero filled
      detail::numeric_initializer<dimension,restrict_iterator,std::numeric_limits<value_t>::is_specialized>::
        apply(restrict_begin());
    }
    
    /// Fast init (garbage), handle with care!
    inline fixed_size_array(const tags::not_initialized) {
      // Garbage memory contents!
    }

    /// Set this to other
    inline fixed_size_array(const self_t& other) {
      fastest::copy<dimension>::apply(other.begin(), restrict_begin());
    }

    /// Set this to other (hoping pointer points  to the right area!)
    inline fixed_size_array(const value_t* other) { 
      // Hoping src dimension is correct 
      SL_REQUIRE("Other exists", other);
      SL_REQUIRE("No sharing", other < begin() && other > end());
      fastest::copy<dimension>::apply(other, restrict_begin());
    }

    /// Explicit init from components (1D)
    inline fixed_size_array(const value_t& i0) {
      SL_REQUIRE("Good dimension", dimension == 1);
      storage_[0] = i0;
    }

    /// Explicit init from components (2D)
    inline fixed_size_array(const value_t& i0, 
                            const value_t& i1) {
      SL_REQUIRE("Good dimension", dimension == 2);
      storage_[0] = i0;
      storage_[1] = i1;
    }

    /// Explicit init from components (3D)
    inline fixed_size_array(const value_t& i0, 
                            const value_t& i1,
                            const value_t& i2) {
      SL_REQUIRE("Good dimension", dimension == 3);
      storage_[0] = i0;
      storage_[1] = i1;
      storage_[2] = i2;
    }

    /// Explicit init from components (4D)
    inline fixed_size_array(const value_t& i0, 
                            const value_t& i1,
                            const value_t& i2,
                            const value_t& i3) {
      SL_REQUIRE("Good dimension", dimension == 4);
      storage_[0] = i0;
      storage_[1] = i1;
      storage_[2] = i2;
      storage_[3] = i3;
    }

    /** Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this,  dimension, v);
    }

  public:  // Indexing
    
    /// is i a good element index?
    static inline bool good_index(const size_t i) {
      return (i< size());
    }

    /// the i-th element
    inline value_t& operator[](size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// the i-th element
    inline const value_t& operator[](size_t i) const {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// the i-th element
    inline value_t& operator()(size_t i) {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// the i-th element
    inline const value_t& operator()(size_t i) const {
      SL_REQUIRE("Good index", good_index(i));
      return( storage_[i] );
    }

    /// pointer to the storage area
    inline const value_t* to_pointer() const {
      return( &(storage_[0]) );
    }

    /// pointer to the storage area
    inline value_t* to_pointer() {
      return( &(storage_[0]) );
    }

  public: // Initializing

    /// fill with c0
    inline void fill(value_t c0) {
      fastest::fill<dimension>::apply(restrict_begin(), c0);
    }

  public: // Searching

    /// index of component with largest value
    inline size_t imax() const {
      const_iterator imax_iter = fastest::max_element<dimension>::apply(begin());
      const size_t result = imax_iter - begin();
      return result;
    }

    /// index of component with largest absolute value
    inline size_t iamax() const { // Loop : avoids multiple calls to abs
      size_t      result = 0;
      value_t  a_res(abs(storage_[result]));
      for (size_t i = 1; i<dimension; i++) {
        value_t a_i(abs(storage_[i]));
        if (a_res < a_i) { 
          result = i; 
          a_res = a_i; 
        }
      }
      return result;
    }

    /// index of component with smallest value
    inline size_t imin() const {
      const_iterator imax_iter = fastest::min_element<dimension>::apply(begin());
      size_t result = imax_iter - begin();
      return result;
    }

    /// index of component with smallest absolute value
    inline size_t iamin() const { // Loop : avoids multiple calls to abs
      size_t      result = 0;
      value_t  a_res(abs(storage_[result]));
      for (size_t i = 1; i<dimension; i++) {
        value_t a_i(abs(storage_[i]));
        if (a_i < a_res) { 
          result = i; 
          a_res = a_i; 
        }
      }
      return result;
    }

  public: // Comparison

    /// -1 if this < t2, +1 if this > t2, 0 otherwise (sequential element comparison)
    int compare(const self_t& t2) const {
      for (size_t i= 0; i< dimension; i++) {
        if (storage_[i] < t2.storage_[i]) {
          return -1;
        } else if (storage_[i] > t2.storage_[i]) {
          return 1;
        }
      }
      return 0;
    }

    /// is this less than other? (sequential element comparison)
    inline bool operator< (const self_t& t2) const {
      return compare(t2)== -1;
    }

    /// is this equal to other?
    inline bool operator== (const self_t& t2) const {
      return compare(t2) == 0;
    }

    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

    /// are all components of this with eps from the corresponding components of t2?
    bool is_epsilon_equal(const self_t& t2,
                          const value_t& eps) const {
      for (size_t i = 0; i<dimension; i++) {
        if (abs(storage_[i]- t2.storage_[i]) > eps) {
          return false;
        }
      }
      return true;
    }

  public: // Search

    /// Dois this contain x?
    bool has(const value_t& x) const {
      bool result = false;
      for (size_t i = 0; !result && i<dimension; ++i) {
        result = (storage_[i] == x);
      }
      return true;
    }

    /// A pair with first equal to has(x) and second equal to the index of x if found 
    std::pair<bool,size_t> find(const value_t& x) const {
      std::pair<bool,size_t> result(false,0);
      for (size_t i = 0; !result.first && i<dimension; ++i) {
        if (storage_[i] == x) {
          result.first = true;
          result.second = i;
        }
      }
      return result;
    }
      
  public: // Interpolation

    /// linear interpolation with other
    template <class T_PARAMETER>
    inline self_t lerp(const self_t& other, T_PARAMETER t) const {
      self_t result = tags::not_initialized();
      fastest::transform<dimension>::apply(begin(), other.begin(), result.restrict_begin(), lerp_op<value_t,T_PARAMETER>(t));
      return result;
    }

  }; // class fixed_size_array

  template <size_t N, class T>
  inline std::size_t hash_value(const fixed_size_array<N,T>& a) {
    std::size_t seed = 0;
    for (std::size_t i=0; i<N; ++i) {
      hash_combine(seed, a[i]);
    }
    return seed;
  }

  template <std::size_t SZ, typename OUT_ET>
  class conv_to< fixed_size_array<SZ, OUT_ET> > {
  public:
    typedef fixed_size_array<SZ, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const fixed_size_array<SZ, IN_ET>& in) {
      result_t result = tags::not_initialized();
      for (std::size_t i=0; i<SZ; ++i) {
	result[i] = static_cast<OUT_ET>(in[i]);
      }
      return result;
    }    
  }; // class conv_to
  
} // namespace sl

// I/O

template <size_t N, class T>
std::ostream& operator <<(std::ostream& s, const sl::fixed_size_array<N,T>& a) {
  for (size_t i = 0; i<a.size(); i++) {
    if (i!=0){
      s << " ";
    }
    s << a[i];
  }
  return s;
}
    
template <size_t N, class T>
std::istream & operator >>(std::istream & s, sl::fixed_size_array<N,T>& a) {
  for (size_t i = 0; i<a.size(); i++) {
    s >> a[i];
  }
  return s;
}


namespace sl {

  /// array of 2 ints
  typedef fixed_size_array<2,int> tuple2i;
  /// array of 3 ints
  typedef fixed_size_array<3,int> tuple3i;
  /// array of 4 ints
  typedef fixed_size_array<4,int> tuple4i;

  /// array of 2 floats
  typedef fixed_size_array<2,float> tuple2f;
  /// array of 3 floats
  typedef fixed_size_array<3,float> tuple3f;
  /// array of 4 floats
  typedef fixed_size_array<4,float> tuple4f;

  /// array of 2 doubles
  typedef fixed_size_array<2,double> tuple2d;
  /// array of 3 doubles
  typedef fixed_size_array<3,double> tuple3d;
  /// array of 4 doubles
  typedef fixed_size_array<4,double> tuple4d;

}

#endif
