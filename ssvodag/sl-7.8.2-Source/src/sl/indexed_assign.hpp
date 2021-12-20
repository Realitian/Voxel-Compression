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
#ifndef SL_INDEXED_ASSIGN_HPP
#define SL_INDEXED_ASSIGN_HPP

#include <sl/indexed_container.hpp>

namespace sl {

  /// Common base class for indexed_assigners
  class indexed_assigner_base {};

  /// Default assigner interface
  template <class T_LHS, class T_LHS_sparsity_t, class T_RHS, class T_RHS_sparsity_t>
  class indexed_assigner: public indexed_assigner_base {
  public:

    static inline void apply(T_LHS& lhs, const T_RHS& rhs) {
      SL_REQUIRE("Good size", lhs.extent() == rhs.extent());

      for (typename T_RHS::const_iterator rhs_iter= rhs.begin();
	   !rhs_iter.off();
	   ++rhs_iter) {
	lhs(rhs_iter.index()) = rhs_iter.value();
      }
    }
  };

  /// Dense-dense assigner interface
  template <class T_LHS, class T_RHS>
  class indexed_assigner<T_LHS, dense_tag, T_RHS, dense_tag>: public indexed_assigner_base {
  public:

    static inline void apply(T_LHS& lhs, const T_RHS& rhs) {
      SL_REQUIRE("Good size", lhs.extent() == rhs.extent());
      for (typename T_LHS::const_iterator rhs_iter =  rhs.begin();
	   !rhs_iter.off();
	   ++rhs_iter) {
	lhs(rhs_iter.index()) = rhs_iter.value();
      }
    }
  };

  /// Dense-sparse assigner interface
  template <class T_LHS, class T_RHS>
  class indexed_assigner<T_LHS, dense_tag, T_RHS, sparse_tag>: public indexed_assigner_base {
  public:

    static inline void apply(T_LHS& lhs, const T_RHS& rhs) {
      SL_REQUIRE("Good size", lhs.extent() == rhs.extent());

      typedef typename T_LHS::value_t value_t;
      const value_t Zero = scalar_math<value_t>::zero();

      typename T_LHS::const_iterator lhs_iter = lhs.begin();

      typename T_RHS::const_sparse_iterator rhs_iter = rhs.sparse_begin();

      while (lhs_iter.index() < rhs_iter.index()) {
	lhs(lhs_iter().index()) = Zero;
	++lhs_iter;
      }
      while (!lhs_iter.off()) {
	lhs(lhs_iter().index()) = rhs_iter.value();
	++rhs_iter;
	while (lhs_iter.index() < rhs_iter.index()) {
	  lhs(lhs_iter().index()) = Zero;
	  ++lhs_iter;
	}
      }
    }
  };

  /// Sparse-dense assigner interface
  template <class T_LHS, class T_RHS>
  class indexed_assigner<T_LHS, sparse_tag, T_RHS, dense_tag>: public indexed_assigner_base {
  public:

    static inline void apply(T_LHS& lhs, const T_RHS& rhs) {
      SL_REQUIRE("Good size", lhs.extent() == rhs.extent());
     
      // No need to clear, since all values will be set!
      for (typename T_RHS::const_iterator rhs_iter= rhs.begin();
	   !rhs_iter.off();
	   ++rhs_iter) {
	lhs.put(rhs_iter.value(), rhs_iter.index());
      }
    }
  };

  /// Sparse-sparse assigner interface
  template <class T_LHS, class T_RHS>
  class indexed_assigner<T_LHS, sparse_tag, T_RHS, sparse_tag>: public indexed_assigner_base {
  public:

    static inline void apply(T_LHS& lhs, const T_RHS& rhs) {
      SL_REQUIRE("Good size", lhs.extent() == rhs.extent());

      lhs.clear(); // Assume fast for sparse

      for (typename T_RHS::const_sparse_iterator rhs_iter= rhs.sparse_begin();
	   !rhs_iter.off();
	   ++rhs_iter) {
	const typename T_LHS::value_t v = rhs_iter.value();
	if (!is_zero(v)) {
	  lhs.put_non_zero(v, rhs_iter.index());
	}
      }
    }
  };

  template <class G_numtype, size_t G_rank, class G_derived, class G_userdefined, class G_discriminant>
  template <class G_numtype2, class G_derived2, class G_userdefined2> 
  inline 
  typename indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant >::derived_t&
  indexed_container<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant >::
  assign_from(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
        
    typedef typename G_derived::sparsity_t  lhs_sparsity_t;
    typedef typename G_derived2::sparsity_t rhs_sparsity_t;

    indexed_assigner<G_derived, lhs_sparsity_t, G_derived2, rhs_sparsity_t>::
      apply(this->derived_ref(), rhs.derived_ref());

    return this->derived_ref();
  };

}


#endif
