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
#ifndef SL_CSTDINT_HPP
#define SL_CSTDINT_HPP

#include <sl/integer.hpp>

SL_COMPILE_TIME_CHECK("int is at least 16 bits", sizeof(int)*8 >= 16);
SL_COMPILE_TIME_CHECK("long is at least 32 bits", sizeof(long)*8 >= 32);

namespace sl {

  typedef int_t<8>::least int8_t;             
  typedef int_t<8>::least int_least8_t;             
  typedef int_t<8>::fast  int_fast8_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(int8_t) * 8 == 8); 

  typedef uint_t<8>::least uint8_t;             
  typedef uint_t<8>::least uint_least8_t;             
  typedef uint_t<8>::fast  uint_fast8_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(uint8_t) * 8 == 8); 
  
  typedef int_t<16>::least int16_t;             
  typedef int_t<16>::least int_least16_t;             
  typedef int_t<16>::fast  int_fast16_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(int16_t) * 8 == 16); 

  typedef uint_t<16>::least uint16_t;             
  typedef uint_t<16>::least uint_least16_t;             
  typedef uint_t<16>::fast  uint_fast16_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(uint16_t) * 8 == 16); 

  typedef int_t<32>::least int32_t;             
  typedef int_t<32>::least int_least32_t;             
  typedef int_t<32>::fast  int_fast32_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(int32_t) * 8 == 32); 

  typedef uint_t<32>::least uint32_t;             
  typedef uint_t<32>::least uint_least32_t;             
  typedef uint_t<32>::fast  uint_fast32_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(uint32_t) * 8 == 32); 

# if HAVE_LONG_LONG
  typedef long long          intmax_t;
  typedef unsigned long long uintmax_t;

#   if (HAVE_LONG_LONG_BITS >= 64)
      typedef int_t<64>::least int64_t;             
      typedef int_t<64>::least int_least64_t;             
      typedef int_t<64>::fast  int_fast64_t;
      SL_COMPILE_TIME_CHECK("Correct size", sizeof(int64_t) * 8 == 64); 

      typedef uint_t<64>::least uint64_t;             
      typedef uint_t<64>::least uint_least64_t;             
      typedef uint_t<64>::fast  uint_fast64_t;
      SL_COMPILE_TIME_CHECK("Correct size", sizeof(uint64_t) * 8 == 64); 
#   endif

# elif HAVE_MSVC_INT64
  typedef __int64          intmax_t;
  typedef unsigned __int64 uintmax_t;

  typedef int_t<64>::least int64_t;             
  typedef int_t<64>::least int_least64_t;             
  typedef int_t<64>::fast  int_fast64_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(int64_t) * 8 == 64); 

  typedef uint_t<64>::least uint64_t;             
  typedef uint_t<64>::least uint_least64_t;             
  typedef uint_t<64>::fast  uint_fast64_t;
  SL_COMPILE_TIME_CHECK("Correct size", sizeof(uint64_t) * 8 == 64); 

# else
  typedef long              intmax_t;
  typedef unsigned long     uintmax_t;

#   if (HAVE_LONG_BITS >= 64)
      typedef int_t<64>::least int64_t;             
      typedef int_t<64>::least int_least64_t;             
      typedef int_t<64>::fast  int_fast64_t;
      SL_COMPILE_TIME_CHECK("Correct size", sizeof(int64_t) * 8 == 64); 

      typedef uint_t<64>::least uint64_t;             
      typedef uint_t<64>::least uint_least64_t;             
      typedef uint_t<64>::fast  uint_fast64_t;
      SL_COMPILE_TIME_CHECK("Correct size", sizeof(uint64_t) * 8 == 64); 
#   endif

# endif  

} // namespace sl

// -------------------------------------------
// Macros for literal constant definition

# define SL_INT8_C(literal_) static_cast<int8_t>(literal_)
# define SL_UINT8_C(literal_) static_cast<uint8_t>(literal_##u)
# define SL_INT16_C(literal_) static_cast<int16_t>(literal_)
# define SL_UINT16_C(literal_) static_cast<uint16_t>(literal_##u)
# if HAVE_INT_BITS >= 32
#   define SL_INT32_C(literal_) static_cast<int32_t>(literal_)
#   define SL_UINT32_C(literal_) static_cast<uint32_t>(literal_##u)
# else 
#   define SL_INT32_C(literal_) static_cast<int32_t>(literal_##L)
#   define SL_UINT32_C(literal_) static_cast<uint32_t>(literal_##uL)
# endif
# if HAVE_LONG_BITS >= 64
#   define SL_INT64_C(literal_) static_cast<int64_t>(literal_##L)
#   define SL_UINT64_C(literal_) static_cast<uint64_t>(literal_##uL)
# elif HAVE_LONG_LONG
#   if HAVE_LONG_LONG_BITS >= 64
#     define SL_INT64_C(literal_) static_cast<int64_t>(literal_##LL)
#     define SL_UINT64_C(literal_) static_cast<uint64_t>(literal_##uLL)
#   endif
# endif

// -------------------------------------------
// Math functions on integers

namespace sl {

  static inline int8_t abs(int8_t x) {
    return (x>=0) ? x : -x;
  }

  static inline int16_t abs(int16_t x) {
    return (x>=0) ? x : -x;
  }

  static inline int32_t abs(int32_t x) {
    return (x>=0) ? x : -x;
  }

# if (HAVE_LONG_BITS >= 64) || (HAVE_LONG_LONG && (HAVE_LONG_LONG_BITS >= 64))
  static inline int64_t abs(int64_t x) {
    return (x>=0) ? x : -x;
  }
# endif
  
  /**
   *  The largest positive integer y such as y*y <= x
   */
  extern uint32_t isqrt(uint32_t x);

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53
   */
  extern int idet2x2_sign(double a, double b, double c, double d);

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Slower version of  integer_det2x2_sign that does not use
   * IEEE arithmetic.
   */
  extern int not_lazy_idet2x2_sign(double a, double b , double c, double d);

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Returns 2 if if one entry is too big or 3 f one entry is not an integer
   */
  extern int checked_idet2x2_sign(double a, double b , double c, double d);

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53
   */
  inline int idet2x2_sign(int32_t a, int32_t b, int32_t c, int32_t d) {
    return idet2x2_sign(double(a), double(b), double(c), double(d));
  }

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Slower version of  integer_det2x2_sign that does not use
   * IEEE arithmetic.
   */
  inline int not_lazy_idet2x2_sign(int32_t a, int32_t b , int32_t c, int32_t d) {
    return not_lazy_idet2x2_sign(double(a), double(b), double(c), double(d));
  }

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Returns 2 if if one entry is too big or 3 f one entry is not an integer
   */
  inline int checked_idet2x2_sign(int32_t a, int32_t b , int32_t c, int32_t d) {
    return checked_idet2x2_sign(double(a), double(b), double(c), double(d));
  }

# if HAVE_LONG_BITS >= 64  

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53
   */
  inline int idet2x2_sign(int64_t a, int64_t b, int64_t c, int64_t d) {
    return idet2x2_sign(double(a), double(b), double(c), double(d));
  }

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Slower version of  integer_det2x2_sign that does not use
   * IEEE arithmetic.
   */
  inline int not_lazy_idet2x2_sign(int64_t a, int64_t b , int64_t c, int64_t d) {
    return not_lazy_idet2x2_sign(double(a), double(b), double(c), double(d));
  }

  /**
   * The sign (-1,0,1) of the 2x2 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Returns 2 if if one entry is too big or 3 f one entry is not an integer
   */
  inline int checked_idet2x2_sign(int64_t a, int64_t b , int64_t c, int64_t d) {
    return checked_idet2x2_sign(double(a), double(b), double(c), double(d));
  }

#endif
  
  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53
   */
  extern int idet3x3_sign(double a, double b, double c, double d, double e, double f, double g, double h, double i);

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Slower version of  integer_det3x3_sign that does not use
   * IEEE arithmetic.
   */
  extern int not_lazy_idet3x3_sign(double a, double b , double c, double d, double e, double f, double g, double h, double i);

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Returns 2 if if one entry is too big or 3 f one entry is not an integer
   */
  extern int checked_idet3x3_sign(double a, double b , double c, double d, double e, double f, double g, double h, double i);

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53
   */
  inline int idet3x3_sign(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e, int32_t f, int32_t g, int32_t h, int32_t i) {
    return idet3x3_sign(double(a), double(b), double(c), double(d), double(e), double(f), double(g), double(h), double(i));
  }

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Slower version of  integer_det3x3_sign that does not use
   * IEEE arithmetic.
   */
  inline int not_lazy_idet3x3_sign(int32_t a, int32_t b , int32_t c, int32_t d, int32_t e, int32_t f, int32_t g, int32_t h, int32_t i) {
    return not_lazy_idet3x3_sign(double(a), double(b), double(c), double(d), double(e), double(f), double(g), double(h), double(i));
  }

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Returns 2 if if one entry is too big or 3 f one entry is not an integer
   */
  inline int checked_idet3x3_sign(int32_t a, int32_t b , int32_t c, int32_t d, int32_t e, int32_t f, int32_t g, int32_t h, int32_t i) {
    return checked_idet3x3_sign(double(a), double(b), double(c), double(d), double(e), double(f), double(g), double(h), double(i));
  }

# if HAVE_LONG_BITS >= 64  

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53
   */
  inline int idet3x3_sign(int64_t a, int64_t b, int64_t c, int64_t d, int64_t e, int64_t f, int64_t g, int64_t h, int64_t i) {
    return idet3x3_sign(double(a), double(b), double(c), double(d), double(e), double(f), double(g), double(h), double(i));
  }

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Slower version of  integer_det3x3_sign that does not use
   * IEEE arithmetic.
   */
  inline int not_lazy_idet3x3_sign(int64_t a, int64_t b , int64_t c, int64_t d, int64_t e, int64_t f, int64_t g, int64_t h, int64_t i) {
    return not_lazy_idet3x3_sign(double(a), double(b), double(c), double(d), double(e), double(f), double(g), double(h), double(i));
  }

  /**
   * The sign (-1,0,1) of the 3x3 determinant of the matrix a,b,c,d.
   * The matrix must contain integer values between -2^53 and 2^53.
   * Returns 2 if if one entry is too big or 3 f one entry is not an integer
   */
  inline int checked_idet3x3_sign(int64_t a, int64_t b , int64_t c, int64_t d, int64_t e, int64_t f, int64_t g, int64_t h, int64_t i) {
    return checked_idet3x3_sign(double(a), double(b), double(c), double(d), double(e), double(f), double(g), double(h), double(i));
  }

#endif
}

#endif // SL_CSTDINT_HPP
