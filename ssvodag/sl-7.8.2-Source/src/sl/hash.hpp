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
#ifndef SL_HASH_HPP
#define SL_HASH_HPP

#include <sl/utility.hpp>

namespace sl {

  template <class T> struct hash;
  
  /// A good hash value for a single byte (by random bit shuffling)
  extern unsigned char hash_byte(const unsigned char byte);

  /// A general purpose hashing function for a byte array
  extern std::size_t hash_bytes(const unsigned char* byte,
				std::size_t n);

  /// A general purpose hashing function for a constant size byte array
  template <std::size_t N>
  struct buffer_hasher {
    static inline std::size_t hash(const unsigned char* byte) {
      const unsigned char P_N = (unsigned char)small_prime(53-(N%54));
      return
        buffer_hasher<N-1>::hash(byte) + std::size_t(P_N) * hash_byte(byte[N-1] ^ P_N);
    }
  };

  template <>
  struct buffer_hasher<0> {
    static inline std::size_t hash(const unsigned char*) {
      return 0;
    }
  };
  
  /// A general purpose hashing function for a constant size byte array
  template <std::size_t N>
  static inline std::size_t hash_bytes(const unsigned char* byte) {
    return buffer_hasher<N>::hash(byte);
  }

  /* ==== Hash values ==== */
  inline std::size_t hash_value(bool x)    { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(int8_t x)  { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(uint8_t x) { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(int16_t x) { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(uint16_t x) { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(int32_t x) { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(uint32_t x) { return static_cast<std::size_t>(x); }
  inline std::size_t hash_value(int64_t x) { return hash_bytes<sizeof(x)>((const unsigned char*)&x); }
  inline std::size_t hash_value(uint64_t x) { return hash_bytes<sizeof(x)>((const unsigned char*)&x); }
  inline std::size_t hash_value(float x) { return hash_bytes<sizeof(x)>((const unsigned char*)&x); }
  inline std::size_t hash_value(double x) { return hash_bytes<sizeof(x)>((const unsigned char*)&x); }

  template <class T>
  inline void hash_combine(std::size_t& seed, T const& v) {
    sl::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }
  
  template <class A, class B>
  inline std::size_t hash_value(std::pair<A, B> const& v) {
    std::size_t seed = 0;
    hash_combine(seed, v.first);
    hash_combine(seed, v.second);
    return seed;
  }

#define SL_HASH_SPECIALIZE(type)			 \
    template <> struct hash<type>			 \
      : public std::unary_function<type, std::size_t>	 \
    { \
      std::size_t operator()(type v) const	\
      {						\
	return sl::hash_value(v);		\
      }						\
    };

#define SL_HASH_SPECIALIZE_REF(type)			 \
  template <> struct hash<type>				 \
    : public std::unary_function<type, std::size_t>	 \
  {							 \
    std::size_t operator()(type const& v) const		 \
    {							 \
      return sl::hash_value(v);				 \
    }							 \
  };

  SL_HASH_SPECIALIZE(bool)
  SL_HASH_SPECIALIZE(int8_t)
  SL_HASH_SPECIALIZE(uint8_t)
  SL_HASH_SPECIALIZE(int16_t)
  SL_HASH_SPECIALIZE(uint16_t)
  SL_HASH_SPECIALIZE(int32_t)
  SL_HASH_SPECIALIZE(uint32_t)
  SL_HASH_SPECIALIZE(int64_t)
  SL_HASH_SPECIALIZE(uint64_t)
  SL_HASH_SPECIALIZE(float)
  SL_HASH_SPECIALIZE(double)

  template <class T>
  struct hash<T*>: public std::unary_function<T*, std::size_t> {
    std::size_t operator()(T* v) const {
      std::size_t x = static_cast<std::size_t>(reinterpret_cast<std::ptrdiff_t>(v));
      return x + (x >> 3);
    }
  };

  template <class T> struct hash: std::unary_function<T, std::size_t> {
    std::size_t operator()(T const& val) const {
      return hash_value(val);
    }
  };

}

#endif

