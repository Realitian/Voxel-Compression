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
#ifndef SL_INDEXED_FUNCTIONS_HPP
#define SL_INDEXED_FUNCTIONS_HPP

#include <sl/indexed.hpp>
#include <sl/math.hpp>

#ifdef _WIN32
#undef max
#undef min
#endif

namespace sl {

  template <
    class G_numtype,
    size_t G_rank,
    class G_derived,
    class G_userdefined,
    class  G_discriminant
    >
  inline G_numtype max(const indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& arg1) {
    typename G_derived::const_iterator it1 = arg1.derived_ref().begin();
    G_numtype result = scalar_math<G_numtype>::finite_lower_bound();
    for (;!it1.off(); ++it1) {
      G_numtype a_ij = arg1(it1.index());
      if (a_ij > result) result = a_ij;
    }
    return result;
  }
  
  template <
    class G_numtype,
    size_t G_rank,
    class G_derived,
    class G_userdefined,
    class  G_discriminant
    >
  inline G_numtype min(const indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& arg1) {
    typename G_derived::const_iterator it1 = arg1.derived_ref().begin();
    G_numtype result = scalar_math<G_numtype>::finite_upper_bound();
    for (;!it1.off(); ++it1) {
      G_numtype a_ij = arg1(it1.index());
      if (a_ij < result) result = a_ij;
    }
    return result;
  }
  
  template <
    class G_numtype,
    size_t G_rank,
    class G_derived,
    class G_userdefined,
    class  G_discriminant
    >
  inline G_numtype amax(const indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& arg1) {
    typename G_derived::const_sparse_iterator it1 = arg1.derived_ref().sparse_begin();
    G_numtype result = scalar_math<G_numtype>::zero();
    for (;!it1.off(); ++it1) {
      G_numtype a_ij = absolute_value(arg1(it1.index()));
      if (a_ij > result) result = a_ij;
    }
    return result;
  }
  
  template <
    class G_numtype,
    size_t G_rank,
    class G_derived,
    class G_userdefined,
    class  G_discriminant
    >
  inline G_numtype amin(const indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& arg1) {
    typename G_derived::const_sparse_iterator it1 = arg1.derived_ref().sparse_begin();
    G_numtype result = scalar_math<G_numtype>::zero();
    for (;!it1.off(); ++it1) {
      G_numtype a_ij = absolute_value(arg1(it1.index()));
      if (a_ij < result) result = a_ij;
    }
    return result;
  }

  template <
    class G_numtype,
    size_t G_rank,
    class G_derived,
    class G_userdefined,
    class  G_discriminant
    >
  inline SL_SUMTYPENAME(G_numtype) energy(const indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& arg1) {
    typedef SL_SUMTYPENAME(G_numtype) result_t;
    typename G_derived::const_sparse_iterator it1 = arg1.derived_ref().sparse_begin();
    result_t result = scalar_math<result_t>::zero();
    for (;!it1.off(); ++it1) {
      result_t a_ij = result_t(arg1(it1.index()));
      result += a_ij*a_ij;
    }
    return result;
  }
  
  template <
    class G_numtype,
    size_t G_rank,
    class G_derived,
    class G_userdefined,
    class  G_discriminant
    >
  inline SL_FLOATTYPENAME(SL_SUMTYPENAME(G_numtype)) mean(const indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& arg1) {
    typedef SL_FLOATTYPENAME(SL_SUMTYPENAME(G_numtype)) result_t;
    typename G_derived::const_sparse_iterator it1 = arg1.derived_ref().sparse_begin();
    result_t result = scalar_math<result_t>::zero();
    std::size_t N = arg1.count();
    if (N) {
      for (;!it1.off(); ++it1) {
      result_t a_ij = result_t(arg1(it1.index()));
      result += a_ij;
      }
      return result/= result_t(N);
    }
    return result;
  }
  
  template <
    class G_numtype1,
    class G_numtype2,
    size_t G_rank,
    class G_derived1,
    class G_derived2,
    class G_userdefined1,
    class G_userdefined2,
    class G_discriminant
    >
  inline SL_FLOATTYPENAME(SL_SUMTYPENAME(SL_PROMOTENAME(G_numtype1,G_numtype2))) rms(const indexed<G_numtype1, G_rank, G_derived1, G_userdefined1, G_discriminant>& arg1,
                                                                                     const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& arg2) {
    typedef SL_FLOATTYPENAME(SL_SUMTYPENAME(SL_PROMOTENAME(G_numtype1,G_numtype2))) result_t;
    typename G_derived1::const_iterator it1 = arg1.derived_ref().begin();
    result_t result = scalar_math<result_t>::zero();
    std::size_t N = arg1.count();
    if (N) {
      for (;!it1.off(); ++it1) {
        result_t a1_ij = result_t(arg1(it1.index()));
        result_t a2_ij = result_t(arg2(it1.index()));
        result += sqr(a1_ij-a2_ij);
      }
      return std::sqrt(result/result_t(N));
    }
    return result;
  }
  
  template <
    class G_numtype1,
    class G_numtype2,
    size_t G_rank,
    class G_derived1,
    class G_derived2,
    class G_userdefined1,
    class G_userdefined2,
    class G_discriminant
    >
  inline SL_FLOATTYPENAME(SL_SUMTYPENAME(SL_PROMOTENAME(G_numtype1,G_numtype2))) amax_diff(const indexed<G_numtype1, G_rank, G_derived1, G_userdefined1, G_discriminant>& arg1,
                                                                                           const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& arg2) {
    typedef SL_FLOATTYPENAME(SL_SUMTYPENAME(SL_PROMOTENAME(G_numtype1,G_numtype2))) result_t;
    typename G_derived1::const_iterator it1 = arg1.derived_ref().begin();
    result_t result = scalar_math<result_t>::zero();
    for (;!it1.off(); ++it1) {
      result_t a1_ij = result_t(arg1(it1.index()));
      result_t a2_ij = result_t(arg2(it1.index()));
      result_t d_ij = a1_ij>a2_ij ? a1_ij-a2_ij : a2_ij-a1_ij;
      result = std::max(result, d_ij);
    }
    return result;
  }
  
  template <
    class G_numtype1,
    class G_numtype2,
    size_t G_rank,
    class G_derived1,
    class G_derived2,
    class G_userdefined1,
    class G_userdefined2,
    class G_discriminant
    >
  inline SL_FLOATTYPENAME(SL_SUMTYPENAME(SL_PROMOTENAME(G_numtype1,G_numtype2))) psnr(const indexed<G_numtype1, G_rank, G_derived1, G_userdefined1, G_discriminant>& arg1,
                                                                                      const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& arg2) {
    typedef SL_FLOATTYPENAME(SL_SUMTYPENAME(SL_PROMOTENAME(G_numtype1,G_numtype2))) result_t;
    typename G_derived1::const_iterator it1 = arg1.derived_ref().begin();
    result_t result      = scalar_math<result_t>::finite_upper_bound();
    std::size_t N = arg1.count();
    if (N) {
      result_t sum_sqr_err = scalar_math<result_t>::zero();
      result_t amax        = scalar_math<result_t>::zero();
      for (;!it1.off(); ++it1) {
        result_t a1_ij = result_t(arg1(it1.index()));
        result_t a2_ij = result_t(arg2(it1.index()));
        sum_sqr_err += sqr(a1_ij-a2_ij);
        amax = std::max(abs(a1_ij), amax);
      }
      result_t rms = std::sqrt(sum_sqr_err/result_t(N));
      result = 
        is_positive(rms)
        ?
        (is_positive(amax) ? result_t(20.0 * log(amax/rms) / log(10.0)) : scalar_math<result_t>::finite_lower_bound())
        :
        (is_positive(amax) ? scalar_math<result_t>::finite_upper_bound() : scalar_math<result_t>::one());
      
    }
    return result;
  }

}

#endif
