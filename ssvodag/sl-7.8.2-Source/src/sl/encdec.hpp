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
#ifndef SL_ENCDEC_HPP
#define SL_ENCDEC_HPP

#include <sl/cstdint.hpp>

namespace sl {

  namespace encdec {
    union convert_float {
      float f;
      uint32_t i;
    };
#if HAVE_INTMAX_BITS >= 64
    union convert_double {
      double d;
      uint64_t i;
    };
#endif
  }

  /// Encode s into dst using a big endian convetion
  template <class T>
  inline size_t be_encode(T s, void* ptr) {
    SL_USEVAR(s); SL_USEVAR(ptr);
    SL_FAIL("NOT SPECIALIZED"); return 0;
  }
  
  /// Encode s into dst using a little endian convetion
  template <class T>
  inline size_t le_encode(T s, void* ptr) {
    SL_USEVAR(s); SL_USEVAR(ptr);
    SL_FAIL("NOT SPECIALIZED"); return 0;
  }

  
  /// Decode s from src using a big endian convetion
  template <class T>
  inline T be_decode(const void* ptr) {
    SL_USEVAR(ptr);
    SL_FAIL("NOT SPECIALIZED"); return 0;
  }

  /// Decode s from src using a little endian convetion
  template <class T>
  inline T le_decode(const void* ptr) {
    SL_USEVAR(ptr);
    SL_FAIL("NOT SPECIALIZED"); return 0;
  }

  // -- UINT8

  template<>
  inline size_t be_encode(uint8_t s, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    dst[0] = s;
    return 1;
  }
  
  template<>
  inline size_t le_encode(uint8_t s, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    dst[0] = s;
    return 1;
  }
  
  template<>
  inline uint8_t be_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    return src[0];
  }

  template<>
  inline uint8_t le_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    return src[0];
  }

  // -- INT8

  template<>
  inline size_t be_encode(int8_t s, void *ptr) {
    return be_encode(uint8_t(s), ptr);
  }
  
  template<>
  inline size_t le_encode(int8_t s, void *ptr) {
    return le_encode(uint8_t(s), ptr);
  }
  
  template<>
  inline int8_t be_decode(const void *ptr) {
    return int8_t(be_decode<uint8_t>(ptr));
  }

  template<>
  inline int8_t le_decode(const void *ptr) {
    return int8_t(le_decode<uint8_t>(ptr));
  }
  
  // -- UINT16
  
  template<>
  inline size_t be_encode(uint16_t s, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    dst[0] = (s >> 8) & 0xFF;
    dst[1] = s & 0xFF;
    return 2;
  }

  template<>
  inline size_t le_encode(uint16_t s, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    dst[0] = s & 0xFF;
    dst[1] = (s >> 8) & 0xFF;
    return 2;
  }

  template<>
  inline uint16_t be_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    return ((unsigned)src[0] << 8) | src[1];
  }
  
  template<>
  inline uint16_t le_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    return src[0] | ((unsigned)src[1] << 8);
  }
  
  // -- INT16

  template<>
  inline size_t be_encode(int16_t s, void *ptr) {
    return be_encode(uint16_t(s), ptr);
  }
  
  template<>
  inline size_t le_encode(int16_t s, void *ptr) {
    return le_encode(uint16_t(s), ptr);
  }
  
  template<>
  inline int16_t be_decode(const void *ptr) {
    return int16_t(be_decode<uint16_t>(ptr));
  }

  template<>
  inline int16_t le_decode(const void *ptr) {
    return int16_t(le_decode<uint16_t>(ptr));
  }

  // -- UINT32
  
  template<>
  inline size_t be_encode(uint32_t i, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    dst[0] = (i >> 24) & 0xFF;
    dst[1] = (i >> 16) & 0xFF;
    dst[2] = (i >> 8) & 0xFF;
    dst[3] = i & 0xFF;
    return 4;
  }

  template<>
  inline size_t le_encode(uint32_t i, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    dst[0] = i & 0xFF;
    dst[1] = (i >> 8) & 0xFF;
    dst[2] = (i >> 16) & 0xFF;
    dst[3] = (i >> 24) & 0xFF;
    return 4;
  }

  template<>
  inline uint32_t be_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    return ((unsigned)src[0] << 24) | ((unsigned)src[1] << 16) |
      ((unsigned)src[2] << 8) | src[3];
  }

  template<>
  inline uint32_t le_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    return src[0] | ((unsigned)src[1] << 8) |
      ((unsigned)src[2] << 16) | ((unsigned)src[3] << 24);
  }

  // -- INT32

  template<>
  inline size_t be_encode(int32_t s, void *ptr) {
    return be_encode(uint32_t(s), ptr);
  }
  
  template<>
  inline size_t le_encode(int32_t s, void *ptr) {
    return le_encode(uint32_t(s), ptr);
  }
  
  template<>
  inline int32_t be_decode(const void *ptr) {
    return int32_t(be_decode<uint32_t>(ptr));
  }

  template<>
  inline int32_t le_decode(const void *ptr) {
    return int32_t(le_decode<uint32_t>(ptr));
  }


#if HAVE_INTMAX_BITS >= 64
  // -- UINT64

  template<>
  inline size_t be_encode(uint64_t i, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    be_encode(uint32_t(i & SL_UINT64_C(0xFFFFFFFF)), dst + 4);
    be_encode(uint32_t((i >> SL_UINT64_C(32)) & SL_UINT64_C(0xFFFFFFFF)), ptr);
    return 8;
  }

  template<>
  inline size_t le_encode(uint64_t i, void *ptr) {
    unsigned char* dst = (unsigned char*)ptr; 
    le_encode(uint32_t(i & SL_UINT64_C(0xFFFFFFFF)), ptr);
    le_encode(uint32_t((i >> SL_UINT64_C(32)) & SL_UINT64_C(0xFFFFFFFF)), dst + 4);
    return 8;
  }

  template<>
  inline uint64_t be_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    uint64_t i;
    i = be_decode<uint32_t>(ptr);
    i <<= SL_UINT64_C(32);
    i |= be_decode<uint32_t>(src + 4);
    return i;
  }

  template<>
  inline uint64_t le_decode(const void *ptr) {
    const unsigned char* src = (const unsigned char*)ptr; 
    uint64_t i;
    i = le_decode<uint32_t>(src + 4);
    i <<= SL_UINT64_C(32);
    i |= le_decode<uint32_t>(ptr);
    return i;
  }

  // -- INT64

  template<>
  inline size_t be_encode(int64_t s, void *ptr) {
    return be_encode(uint64_t(s), ptr);
  }
  
  template<>
  inline size_t le_encode(int64_t s, void *ptr) {
    return le_encode(uint64_t(s), ptr);
  }
  
  template<>
  inline int64_t be_decode(const void *ptr) {
    return int64_t(be_decode<uint64_t>(ptr));
  }

  template<>
  inline int64_t le_decode(const void *ptr) {
    return int64_t(le_decode<uint64_t>(ptr));
  }

  
#endif
  
  // -- FLOAT

  template<>
  inline size_t le_encode(const float f, void *ptr) {
    union encdec::convert_float c;
    c.f = f;
    return le_encode(c.i, ptr);
  }

  template<>
  inline size_t be_encode(const float f, void *ptr) {
    union encdec::convert_float c;
    c.f = f;
    return be_encode(c.i, ptr);
  }

  template<>
  inline float le_decode(const void *ptr) {
    union encdec::convert_float c;
    c.i = le_decode<uint32_t>(ptr);
    return c.f;
  }
  
  template<>
  inline float be_decode(const void *ptr) {
    union encdec::convert_float c;
    c.i = be_decode<uint32_t>(ptr);
    return c.f;
  }

  // -- DOUBLE

#if HAVE_INTMAX_BITS>=64
  template<>
  inline size_t le_encode(const double d, void *ptr) {
    union encdec::convert_double c;
    c.d = d;
    return le_encode(c.i, ptr);
  }

  template<>
  inline size_t be_encode(const double d, void *ptr) {
    union encdec::convert_double c;
    c.d = d;
    return be_encode(c.i, ptr);
  }

  template<>
  inline double le_decode(const void *ptr) {
    union encdec::convert_double c;
    c.i = le_decode<uint64_t>(ptr);
    return c.d;
  }
  
  template<>
  inline double be_decode(const void *ptr) {
    union encdec::convert_double c;
    c.i = be_decode<uint64_t>(ptr);
    return c.d;
  }

#endif
  
}

#endif
