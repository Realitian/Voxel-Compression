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
#ifndef SL_INDEXED_COMPARE_HPP
#define SL_INDEXED_COMPARE_HPP

#include <sl/indexed.hpp>

/// Indexed comparison: -1 if arg1+eps < arg2, +1 if arg1-eps > arg2, 0 otherwise (sequential element comparison)
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
	    const SL_PROMOTENAME(G_numtype1,G_numtype2)& eps) {
  SL_REQUIRE("Good size", arg1.extent() == arg2.extent());

  typedef SL_PROMOTENAME(G_numtype1,G_numtype2) value_t;

  typedef typename G_derived1::sparsity_t arg1_sparsity_t;
  typedef typename G_derived2::sparsity_t arg2_sparsity_t;
  
  if (sl::is_same<arg1_sparsity_t,sl::sparse_tag>::value &&
      sl::is_same<arg2_sparsity_t,sl::sparse_tag>::value) {
    typename G_derived1::const_sparse_iterator it1      = arg1.derived_ref().sparse_begin();
    typename G_derived2::const_sparse_iterator it2      = arg2.derived_ref().sparse_begin();

    for (;!it1.off() && !it2.off(); ++it1, ++it2) {
      if (it1.index() == it2.index()) {
	const value_t x_i = static_cast<value_t>(it1.value());
	const value_t y_i = static_cast<value_t>(it2.value());
	if ((x_i > y_i) && (x_i - y_i > eps)) {
	  return 1;
	} else if ((x_i < y_i) && (y_i - x_i > eps)) {
	  return -1;
	}
	++it1; ++it2;
      } else if (it1.index() < it2.index()) {
	const value_t x_i = static_cast<value_t>(it1.value());
	const value_t y_i = sl::scalar_math<value_t>::zero();
	if ((x_i > y_i) && (x_i - y_i > eps)) {
	  return 1;
	} else if ((x_i < y_i) && (y_i - x_i > eps)) {
	  return -1;
	}
	++it1;
      } else {
	const value_t x_i = sl::scalar_math<value_t>::zero();
	const value_t y_i = static_cast<value_t>(it2.value());
	if ((x_i > y_i) && (x_i - y_i > eps)) {
	  return 1;
	} else if ((x_i < y_i) && (y_i - x_i > eps)) {
	  return -1;
	}
	++it2;
      }	
    }
  } else {
    typename G_derived1::const_iterator it1 = arg1.derived_ref().begin();
    typename G_derived2::const_iterator it2 = arg2.derived_ref().begin();
    for (;!it1.off(); ++it1,++it2) {
      const value_t x_i = static_cast<value_t>(it1.value());
      const value_t y_i = static_cast<value_t>(it2.value());
      if ((x_i > y_i) && (x_i - y_i > eps)) {
	return 1;
      } else if ((x_i < y_i) && (y_i - x_i > eps)) {
	return -1;
      }
    }
  }
  return 0;
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
	    const sl::indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& arg2) {
  return compare(arg1, arg2, sl::scalar_math<SL_PROMOTENAME(G_numtype1,G_numtype2)>::zero());
}

#endif
