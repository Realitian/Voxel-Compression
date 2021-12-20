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
#ifndef NUMERIC_TRAITS_HPP
#define NUMERIC_TRAITS_HPP

#include <sl/cstdint.hpp> // for int8_t, ...
#include <sl/generative_types.hpp> 
#include <sl/type_traits.hpp> 
#include <string> // to implement feature "what"

#define SL_SUMTYPE(X)    sl::numeric_traits<X>::T_sumtype
#define SL_DIFFTYPE(X)   sl::numeric_traits<X>::T_difftype
#define SL_FLOATTYPE(X)  sl::numeric_traits<X>::T_floattype
#define SL_SIGNEDTYPE(X) sl::numeric_traits<X>::T_signedtype

#define SL_SUMTYPENAME(X)    typename sl::numeric_traits<X>::T_sumtype
#define SL_DIFFTYPENAME(X)   typename sl::numeric_traits<X>::T_difftype
#define SL_FLOATTYPENAME(X)  typename sl::numeric_traits<X>::T_floattype
#define SL_SIGNEDTYPENAME(X) typename sl::numeric_traits<X>::T_signedtype

#define SL_PROMOTE(X,Y)      sl::promote_traits<X,Y>::type
#define SL_PROMOTENAME(X,Y)  typename sl::promote_traits<X,Y>::type

namespace sl {

  /// Basic numeric traits for type P_numtype
  template<class P_numtype>
  class numeric_traits {
  public:
    static std::string what() { return "unknown_type"; }
    typedef P_numtype T_sumtype;     /**< Type to be used for summing */
    typedef P_numtype T_difftype;    /**< Type to be used for difference */
    typedef P_numtype T_floattype;   /**< Type to be used for real number calculations */
    typedef P_numtype T_signedtype;  /**< Type to be used for signed calculations */
    enum { has_trivial_constructor = 0 };    // Assume the worst
    enum { is_specialized = 0 };     // Default to not specialized!
    enum { precision_rank = 30000 }; // "Weight" for type promotion: by default: HUGE!
  };
} // namespace sl

#define SL_DECL_NUMERIC_TRAITS(STR,X,Y,Z,W,U,C,P)                   \
  namespace sl {                                                    \
    template<>                                                      \
    class numeric_traits< X > {                                     \
    public:                                                         \
        static std::string what() { return STR; }                   \
        typedef Y T_sumtype;                                        \
        typedef Z T_difftype;                                       \
        typedef W T_floattype;                                      \
        typedef U T_signedtype;                                     \
        enum { has_trivial_constructor = C };                       \
        enum { is_specialized = 1 };                                \
        enum { precision_rank = P };                                \
        SL_COMPILE_TIME_CHECK("Good rank", (precision_rank < 30000)); \
    };                                                               \
  } // namespace sl


// Standard integer types

//                       STR       TP          SUM    DIFF    FLOAT   SIGN  TC

SL_DECL_NUMERIC_TRAITS("int8_t",   int8_t,  int16_t, int16_t, float,  int8_t, 1,   100);
SL_DECL_NUMERIC_TRAITS("uint8_t", uint8_t, uint16_t, int16_t, float, int16_t, 1,   200);

SL_DECL_NUMERIC_TRAITS("int16_t",  int16_t,  int32_t, int32_t, float, int16_t, 1,  300);
SL_DECL_NUMERIC_TRAITS("uint16_t",uint16_t, uint32_t, int32_t, float, int32_t, 1,  400);

# if HAVE_INTMAX_BITS == 32
SL_DECL_NUMERIC_TRAITS("int32_t",  int32_t,  int32_t, int32_t, float, int32_t, 1,  500);
SL_DECL_NUMERIC_TRAITS("uint32_t",uint32_t, uint32_t, int32_t, float, int32_t, 1,  600);
# else
SL_DECL_NUMERIC_TRAITS("int32_t",  int32_t,  int64_t, int64_t, float, int32_t, 1,  500);
SL_DECL_NUMERIC_TRAITS("uint32_t",uint32_t, uint64_t, int64_t, float, int64_t, 1,  600);
SL_DECL_NUMERIC_TRAITS("int64_t",  int64_t,  int64_t, int64_t, double, int64_t, 1, 700);
SL_DECL_NUMERIC_TRAITS("uint64_t",uint64_t, uint64_t, int64_t, double, int64_t, 1, 800);
# endif
    
// Standard floating point types
//                       STR       TP      SUM    DIFF    FLOAT   SIGN  TC
SL_DECL_NUMERIC_TRAITS("float",   float,  float,  float,  float,  float, 1, 1000);
SL_DECL_NUMERIC_TRAITS("double", double, double, double, double, double, 1, 2000);

// Type promotion

namespace sl {

  /// Trait class for type promotion.
  template<class T1_orig, class T2_orig>
  struct promote_traits {
    // 1. Make both signed if at least one is signed
    typedef typename
      gen_if<
      is_same<T2_orig, SL_SIGNEDTYPENAME(T2_orig)>::value,
      SL_SIGNEDTYPENAME(T1_orig),
      T1_orig>::type T1; 
	
    typedef typename
      gen_if<
      is_same<T1_orig, SL_SIGNEDTYPENAME(T1_orig)>::value,
      SL_SIGNEDTYPENAME(T2_orig),
      T2_orig>::type T2; 
					      
    // 2. Promote to the larger type
    typedef typename
    gen_if< (((int)numeric_traits<T1>::precision_rank) > ((int)numeric_traits<T2>::precision_rank)),
	      T1,
	      typename gen_if<(((int)numeric_traits<T1>::precision_rank) < ((int)numeric_traits<T2>::precision_rank)),
		       T2,
		       typename gen_if< (sizeof(T1) > sizeof(T2)), T1, T2>::type
		    >::type
    >::type  type;
  }; // promote_traits

} // namespace sl
    

#endif



