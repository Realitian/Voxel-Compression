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
#ifndef SL_FLOAT_CAST_HPP
#define SL_FLOAT_CAST_HPP

#include <sl/config.hpp>

#ifndef _ISOC9X_SOURCE
#error "<sl/float_cast.hpp> must be compiled with -D_ISOC9X_SOURCE=1 -D_ISOC99_SOURCE=1 -D__USE_ISOC9X=1 -D__USE_ISOC99=1"
#endif

// Include math.h while ensuring we have access to ISO C99 features
#define	_ISOC9X_SOURCE	1
#define _ISOC99_SOURCE	1
#define	__USE_ISOC9X	1
#define	__USE_ISOC99	1

#include <cmath>

namespace sl {

#if (HAVE_LRINT && HAVE_LRINTF)
  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(double x) {
    return lrint(x);
  }
  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(float x) {
    return lrintf(x);
  }

#elif (defined (WIN32) || defined (_WIN32))
// Win32 doesn't seem to have these functions. 
// Therefore implement inline versions of these functions here.
// The following code comes from  Erik de Castro Lopo <erikd AT mega-nerd DOT com>
// Extended by Samuel Iacolina for 64bit targets

#ifdef _M_X64

#include <emmintrin.h>

  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(double x) {

    return _mm_cvtsd_si32(_mm_load_sd(&x));
  }
  
  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(float x) {

    return _mm_cvtss_si32(_mm_load_ss(&x));
  }

#elif defined(_M_IX86)

  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(double x) {
    long int result;
    _asm
      {fld x
         fistp result
         };
    return result;
  }
  
  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(float x) {
    long int result;
    _asm
      {	fld x
          fistp result
          };
			
    return result;
  }
  
#else
#warning "We don't have the functions lrint() and lrintf (), and the targe is not X64 or IX86."
#warning "Replacing these functions with a standard C cast."

  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(double x) {
    return ((long int)(x));
  }
  
  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(float x) {
    return ((long int)(x));
  }

#endif

#else
#warning "Don't have the functions lrint() and lrintf ()."
#warning "Replacing these functions with a standard C cast."

  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(double x) {
    return ((long int)(x));
  }
  
  /// Round x to the  nearest integer value, using  the current rounding direction.
  static inline long int fast_round_to_integer(float x) {
    return ((long int)(x));
  }

#endif

} // namespace sl


#endif
