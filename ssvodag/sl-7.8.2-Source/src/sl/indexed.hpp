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
#ifndef SL_INDEXED_HPP
#define SL_INDEXED_HPP

#include <sl/math.hpp>
#include <sl/numeric_traits.hpp>
#include <sl/xpr.hpp>
#include <sl/indexed_iterator.hpp>
#include <sl/serializer.hpp>
#include <sl/std_serializer.hpp>

//#include <sl/tensor_indexing.hpp>
//#include <sl/xpr_scalar.hpp>

namespace sl {

  // ----- Init values - move elsewhere?

  namespace detail {

    template <class T, bool is_numeric>
    class init_value_helper {
    public: 
      static const bool init_value_is_default = true;

      static T value() { return T(); }
    };
    
    template <class T>
    class init_value_helper<T, true> {
    public:
      static const bool init_value_is_default = false;
      static T value() { return scalar_math<T>::zero(); }
    };      
    
    template <class T>
    class init_value {
    public:
      static const bool init_value_is_default = init_value_helper<T, std::numeric_limits<T>::is_specialized>::init_value_is_default;
      static T value() { return init_value_helper<T, std::numeric_limits<T>::is_specialized>::value(); }
    };
 
  }

  // ----- Sparsity tags - move elsewhere?
    
  /// Common base class for sparsity types
  class sparsity_tag_base {};

  /// Sparse vectors or expressions (iterating on non-zero elements is faster)
  class sparse_tag: public sparsity_tag_base {};

  /// Dense vectors or expressions (iterating on all elements is faster)
  class dense_tag : public sparsity_tag_base {};

  // ----- Forward declarations
  template <
    class  G_indexed, 
    size_t N0,
    size_t N1 = 0,
    size_t N2 = 0,
    size_t N3 = 0,
    size_t N4 = 0,
    size_t N5 = 0,
    size_t N6 = 0,
    size_t N7 = 0,
    size_t N8 = 0,
    size_t N9 = 0,
    size_t N10 = 0
  >
  class indexed_subscript_remapping;
  
  template <
    class  G_numtype, 
    size_t G_rank, 
    class  G_derived, 
    class  G_userdefined,
    class  G_discriminant
  >
  class indexed;
}

  template <
    class  G_numtype1, 
    class  G_numtype2, 
    size_t G_rank, 
    class  G_derived1, 
    class  G_derived2, 
    class  G_userdefined1,
    class  G_userdefined2,
    class  G_discriminant
    >
  int compare(const sl::indexed<G_numtype1, G_rank, G_derived1, G_userdefined1, G_discriminant>& arg1,
	      const sl::indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& arg2,
	      const SL_PROMOTENAME(G_numtype1,G_numtype2)& eps = sl::scalar_math<SL_PROMOTENAME(G_numtype1,G_numtype2)>::zero());

namespace sl {

  // ----- Class definition

  /**
   *  N-dimensional indexed entities
   */
  template <
    class  G_numtype, 
    size_t G_rank, 
    class  G_derived, 
    class  G_userdefined,
    class  G_discriminant
  >
  class indexed {
    
  public:
    
    typedef indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>   this_t;
    typedef G_derived                                                              derived_t;
    typedef G_userdefined                                                          user_defined_data_t;
    typedef G_numtype                                                              value_t;
    typedef G_discriminant                                                         discriminant_t;

    typedef ::sl::index<G_rank>                                                    subscript_t;
    
    typedef indexed_const_iterator_wrapper<derived_t>                              const_iterator;
    typedef indexed_const_sparse_iterator_wrapper<derived_t>                       const_sparse_iterator;

    /// The rank of the object
    enum { rank_c = G_rank };

  public:

    /// True iff the dimensions are adaptable inside an expression
    enum { xpr_dynamic_dimensions = false };

  protected:

    /// User defined data. Declared here to remove the space overhead associated to empty base classes.
    SL_DECLARE_GENERIC_SUPERCLASS_DATA(user_defined_data_t,data_);
    
  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << data_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> data_;
    }

  public:
    
    /// Default init (no operation)
    explicit inline indexed() {
    }

    /// Data init
    explicit inline indexed(const user_defined_data_t& data): data_(data) {
    }

  public:

    /// Features to access this as a derived_t pointer
    SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);

  public:
    
    /// The rank of this entity (0: scalars, 1: vectors, 2: matrices, ...)
    inline size_t rank() const {
      return rank_c;
    }

    /// Is i a valid rank index for this?
    inline bool good_rank_index(size_t i) const {
      return i<rank_c;
    }

    /// The number of indexed elements for each of the ranks
    //  SUBCLASS RESPONSIBILITY
    inline subscript_t extent() const {
      return derived_ref().extent();
    }

    inline size_t extent(size_t i) const {
      SL_REQUIRE("Good subscript", this->good_subscript(i));
      return extent()[i];
    }
    
    /// Is idx a good subscript?
    inline bool good_subscript(const subscript_t& idx) const {
      return idx.all_lt(extent());
    }

    /// The total number of elements
    inline size_t count() const {
      return extent().element_count();
    }

  public: // Dynamic expression bounds

    /**
     *  The upper bound for rank r when this is part of an 
     *  expresssion. This is extent(r) if r has a fixed dimension,
     *  XPR_DYNAMIC_BOUND otherwise.
     */
    // SUBCLASS RESPONSIBILITY
    inline size_t xpr_ubound(size_t r) const {
      return derived_ref().xpr_ubound(r);
    }

    /// Is the upper bound for rank i dynamic?
    inline size_t xpr_has_dynamic_ubound(size_t r) const {
      return xpr_is_dynamic_bound(xpr_ubound(r));
    }

    /** 
     *  Remapping of a generalized subscript to the corresponding subscript of this.
     *  All missing indices and all indices corresponding to a dynamically sized dimension
     *  are set to zero.
     */
    template <size_t G_rank2>
    inline subscript_t xpr_from_generalized_subscript(const ::sl::index<G_rank2>& idx) const {
      SL_REQUIRE("Dynamic", xpr_dynamic_dimensions || G_rank2 <= G_rank);
      if (G_rank2 == G_rank && !xpr_dynamic_dimensions) { // Help optimizer with compile-time constant
	return idx;
      } else {
	subscript_t result; 
	for (size_t i=0; i<min(G_rank2,G_rank); i++) {
	  result[i] = xpr_has_dynamic_ubound(idx[i]) ? 0 : idx[i];
	}
	return result;
      }
    }

    /// The element referenced by the generalized subscript idx
    template <size_t G_rank2>
    inline value_t xpr_generalized_item(const ::sl::index<G_rank2>& idx) const {
      SL_REQUIRE("Dynamic", xpr_dynamic_dimensions || G_rank2 <= G_rank);
      return derived_ref().item(xpr_from_generalized_subscript(idx));
    }
      
  public: // Element access
    
    /// The element referenced by subscript idx
    // SUBCLASS RESPONSIBILITY
    inline value_t item(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return derived_ref().item(idx);
    }

    /// The element referenced by subscript idx
    inline value_t operator()(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return derived_ref().item(idx);
    }

    /// The element referenced by subscript idx (alias of operator())
    inline value_t operator[](const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return derived_ref().item(idx);
    }

    /// Iterator pointing to the beginning of the element sequence
    //  SUBCLASS RESPONSIBILITY
    inline const_iterator begin() const {
      return const_iterator(derived_ref().begin());
    }

    /// Iterator pointing to the beginning of the non zero element sequence
    //  SUBCLASS RESPONSIBILITY
    inline const_sparse_iterator sparse_begin() const {
      return const_sparse_iterator(derived_ref().sparse_begin());
    }

  public: // Helpers for indexing

    // TODO: SL_DECLARE_EXPLICIT_CONST_SUBSCRIPTING

    inline value_t operator()(size_t i0) const  {
      return operator()(subscript_t(i0));
    }

    inline value_t operator()(size_t i0, size_t i1) const {
      return operator()(subscript_t(i0, i1));
    }

    inline value_t operator()(size_t i0, size_t i1, size_t i2) const {
      return operator()(subscript_t(i0, i1, i2));
    }

    inline value_t operator()(size_t i0, size_t i1, size_t i2, size_t i3) const {
      return operator()(subscript_t(i0, i1, i2, i3));
    }

  public: // Subscript remapping

#define SL_DECLARE_SUBSCRIPT_REMAPPING(derived_t) \
    template<size_t N0> \
    inline indexed_subscript_remapping<derived_t, N0> \
    operator() (index_placeholder<N0>) const { \
      return indexed_subscript_remapping<derived_t,N0>(derived_ref()); \
    }\
    template<size_t N0, size_t N1> \
    inline indexed_subscript_remapping<derived_t, N0, N1> \
    operator() (index_placeholder<N0>, index_placeholder<N1>) const { \
      return indexed_subscript_remapping<derived_t,N0,N1>(derived_ref()); \
    }\
    template<size_t N0, size_t N1, size_t N2> \
    inline indexed_subscript_remapping<derived_t, N0, N1, N2> \
    operator() (index_placeholder<N0>, index_placeholder<N1>, index_placeholder<N2>) const { \
      return indexed_subscript_remapping<derived_t,N0,N1,N2>(derived_ref()); \
    }\
    template<size_t N0, size_t N1, size_t N2, size_t N3> \
    inline indexed_subscript_remapping<derived_t, N0, N1, N2, N3> \
    operator() (index_placeholder<N0>, index_placeholder<N1>, index_placeholder<N2>, index_placeholder<N3>) const { \
      return indexed_subscript_remapping<derived_t,N0,N1,N2,N3>(derived_ref()); \
    }\
    template<size_t N0, size_t N1, size_t N2, size_t N3, size_t N4> \
    inline indexed_subscript_remapping<derived_t, N0, N1, N2, N3, N4> \
    operator() (index_placeholder<N0>, index_placeholder<N1>, index_placeholder<N2>, index_placeholder<N3>, index_placeholder<N4>) const { \
      return indexed_subscript_remapping<derived_t,N0,N1,N2,N3,N4>(derived_ref()); \
    }\
    template<size_t N0, size_t N1, size_t N2, size_t N3, size_t N4, size_t N5> \
    inline indexed_subscript_remapping<derived_t, N0, N1, N2, N3, N4, N5> \
    operator() (index_placeholder<N0>, index_placeholder<N1>, index_placeholder<N2>, index_placeholder<N3>, index_placeholder<N4>, index_placeholder<N5>) const { \
      return indexed_subscript_remapping<derived_t,N0,N1,N2,N3,N4,N5>(derived_ref()); \
    }\
    template<size_t N0, size_t N1, size_t N2, size_t N3, size_t N4, size_t N5, size_t N6> \
    inline indexed_subscript_remapping<derived_t, N0, N1, N2, N3, N4, N5, N6> \
    operator() (index_placeholder<N0>, index_placeholder<N1>, index_placeholder<N2>, index_placeholder<N3>, index_placeholder<N4>, index_placeholder<N5>, index_placeholder<N6>) const { \
      return indexed_subscript_remapping<derived_t,N0,N1,N2,N3,N4,N5,N6>(derived_ref()); \
    }\

    //    SL_DECLARE_SUBSCRIPT_REMAPPING(derived_t);

  public: // Comparison
    
    // TODO: Replace with comparison expressions

    /// -1 if this < other, +1 if this > other, 0 otherwise (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    int compare(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other, 
		const SL_PROMOTENAME(G_numtype, G_numtype2)& eps) const {
      SL_REQUIRE("Good size", this->extent() == other.extent());
      return ::compare(*this, other, eps);
    }
    
    /// -1 if this < other, +1 if this > other, 0 otherwise (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline int compare(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      SL_REQUIRE("Good size", extent() == other.extent());
      SL_PROMOTENAME(G_numtype, G_numtype2) zero_eps = sl::scalar_math< SL_PROMOTENAME(G_numtype, G_numtype2) >::zero();
      return this->compare(other, zero_eps);
    }
     
    /// is this less than other? (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline bool operator<(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      return this->compare(other) < 0;
    }
    
    /// is this greater than other? (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline bool operator>(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      return this->compare(other) > 0;
    }

    /// is this equal to other? (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline bool operator==(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      return this->compare(other) == 0;
    }
    
    /// is this not equal to other? (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline bool operator!=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      return this->compare(other) != 0;
    }
      
    /// is this less than or equal to other? (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline bool operator<=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      return this->compare(other) <= 0;
    }

    /// is this greater than or equal to other? (sequential element comparison)
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline bool operator>=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& other) const {
      return this->compare(other) >= 0;
    }

  };
}
 
//---------------------------------------------------------------------------
// Input/output

#include <sl/indexed_io.hpp>

//---------------------------------------------------------------------------
// Comparison

#include <sl/indexed_compare.hpp>

//---------------------------------------------------------------------------
// Math. expressions

//#include <sl/indexed_unary.hpp>
//#include <sl/indexed_binary.hpp>

//---------------------------------------------------------------------------
// Subscript remapping

#include <sl/indexed_subscript_remapping.hpp>


#endif
