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
#ifndef SL_BITOPS_HPP
#define SL_BITOPS_HPP

#include <sl/config.hpp>
#include <sl/utility.hpp>
#undef log2
// FIXME: Windows hack

namespace sl {

  namespace detail {

    template <class G_int, std::size_t G_pattern, std::size_t G_width, std::size_t N>
    struct bitpattern {
      static const G_int value = (bitpattern<G_int,G_pattern,G_width,N/2+N%2>::value<<G_int(G_width*(N/2))
				  | bitpattern<G_int,G_pattern,G_width,N/2>::value);
    };

    template <class G_int, std::size_t G_pattern, std::size_t G_width>
    struct bitpattern<G_int,G_pattern,G_width,0> {
      static const G_int value = 0;
    };

    template <class G_int, std::size_t G_pattern, std::size_t G_width>
    struct bitpattern<G_int,G_pattern,G_width,1> {
      static const G_int value = G_pattern;
    };
    
    /**
     *  Implementation of bit operations on integers of type G_int,
     *  using manipulations on G_pattern_bits bits at a time.
     *  See Freed, Edwin E.  1983.  "Binary Magic Numbers --- Some
     *  Applications and Algorithms", Dr Dobb's Journal, 8(4: April), pp
     *  24--37.
     */
    template <class G_int, std::size_t G_pattern_bits>
    class bitops {
    public:
      /// The integer type
      typedef G_int value_t;
      
      /// The number of bits in value_t
      enum { bits  = 8 * sizeof(value_t) };
      
      /// The number of bytes in value_t
      enum { bytes = sizeof(value_t) };

      /// The number of consecutive bits with the same value in the pattern
      enum { pattern_bits = G_pattern_bits };

      /// the binary magic number bit pattern
      static const G_int pattern      = bitpattern<G_int,
						   bitpattern<G_int,1,1,pattern_bits>::value,
						   pattern_bits*2,
						   bits/(pattern_bits*2)>::value;

      /// The binary magic number bit mirror pattern
      static const G_int mirror_pattern = ~(bitpattern<G_int,
					               bitpattern<G_int,1,1,pattern_bits>::value,
					               pattern_bits*2,
					               bits/(pattern_bits*2)>::value);

    public:
      
      static value_t one_count(value_t n);
      static value_t reverse(value_t n);
      static value_t low_order_zero_count(value_t n, std::size_t c);
      static value_t log2(value_t n, value_t c);
      static value_t parity(value_t n);
      static value_t interleaved(value_t value);
      static value_t deinterleaved(value_t code);
      static value_t binary_from_gray(value_t g);
   };

    /// Specialization of bitops for G_pattern_bits = 1
    template <class G_int> 
    class bitops<G_int, 1> {
    public:
      /// The integer type
      typedef G_int value_t;
      
      /// The number of bits in value_t
      enum { bits  = 8 * sizeof(value_t) };
      
      /// The number of bytes in value_t
      enum { bytes = sizeof(value_t) };

      /// The number of consecutive bits with the same value in the pattern
      enum { pattern_bits = 1 };

      /// the binary magic number bit pattern
      static const G_int pattern = bitpattern<G_int,1,2,bits/2>::value;

      /// The binary magic number bit mirror pattern
      static const G_int mirror_pattern   = ~(bitpattern<G_int,1,2,bits/2>::value);

    public:
      
      static value_t one_count(value_t n);
      static value_t reverse(value_t n);
      static value_t low_order_zero_count(value_t n, std::size_t c);
      static value_t log2(value_t n, value_t c);
      static value_t parity(value_t n);
      static value_t interleaved(value_t value);
      static value_t deinterleaved(value_t code);
      static value_t binary_from_gray(value_t g);
    };

    template <class value_t> 
    inline value_t bitops<value_t,1>::one_count(value_t n) {
      return ((n >> 1) & pattern) + (n & pattern);
    }

    template <class value_t, std::size_t pattern_bits> 
    inline value_t bitops<value_t,pattern_bits>::one_count(value_t n) {
      n = bitops<value_t,pattern_bits/2>::one_count(n);
      return ((n >> pattern_bits) & pattern) + (n & pattern);
    }

    template <class value_t> 
    inline value_t bitops<value_t,1>::reverse(value_t n) {
      return ((n >> 1) & pattern) | ((n << 1) & mirror_pattern); 
    }

    template <class value_t, std::size_t pattern_bits> 
    inline value_t bitops<value_t,pattern_bits>::reverse(value_t n) {
      n = bitops<value_t,pattern_bits/2>::reverse(n);
      return ((n >> pattern_bits) & pattern) | ((n << pattern_bits) & mirror_pattern);
    }
    
    template <class value_t> 
    inline value_t bitops<value_t,1>::low_order_zero_count(value_t n, std::size_t c) {
      if (n & pattern) { n <<= 1; c--; }
      return n ? c-1 : c;
    }

    template <class value_t, std::size_t pattern_bits> 
    inline value_t bitops<value_t,pattern_bits>::low_order_zero_count(value_t n, std::size_t c) {
      return (n & pattern)
	? bitops<value_t,pattern_bits/2>::low_order_zero_count(n << pattern_bits, c - pattern_bits)
	: bitops<value_t,pattern_bits/2>::low_order_zero_count(n, c);
    }

    template <class value_t> 
    inline value_t bitops<value_t,1>::log2(value_t n, value_t c) {
      return c | (n & mirror_pattern ? 1 : 0); 
    }

    template <class value_t, std::size_t pattern_bits> 
    inline value_t bitops<value_t,pattern_bits>::log2(value_t n, value_t c) {
      return (n & mirror_pattern)
	? bitops<value_t,pattern_bits/2>::log2 (n >> pattern_bits, c | pattern_bits)
	: bitops<value_t,pattern_bits/2>::log2 (n, c);
    }

    template <class value_t> 
    inline value_t bitops<value_t,1>::parity(value_t n) {
      return  n ^ (n >> 1);
    }

    template <class value_t, std::size_t pattern_bits> 
    inline value_t bitops<value_t,pattern_bits>::parity(value_t n) {
      return bitops<value_t,pattern_bits/2>::parity (n ^ (n >> pattern_bits));
    }

    template <class value_t>
    inline value_t bitops<value_t,1>::interleaved(value_t value) {
      return (value | (value << 1)) & pattern;
    }

    template <class value_t, std::size_t pattern_bits>
    inline value_t bitops<value_t,pattern_bits>::interleaved(value_t value) {
      return bitops<value_t,pattern_bits/2>::interleaved((value | (value << pattern_bits)) & pattern);
    }

    template <class value_t>
    inline value_t bitops<value_t,1>::deinterleaved(value_t value) {
      value &= pattern; 
      return value | (value >> 1);
    }

    template <class value_t, std::size_t pattern_bits>
    inline value_t bitops<value_t,pattern_bits>::deinterleaved(value_t value) {
      value = bitops<value_t,pattern_bits/2>::deinterleaved(value) & pattern;
      return value | (value >> pattern_bits);
    }

    template <class value_t>
    inline value_t bitops<value_t,1>::binary_from_gray(value_t value) {
      value^= (value >> 1);
      return value;
    }

    template <class value_t, std::size_t pattern_bits>
    inline value_t bitops<value_t,pattern_bits>::binary_from_gray(value_t value) {
      value ^= (value >> pattern_bits);
      return bitops<value_t,pattern_bits/2>::binary_from_gray(value);
    }

  } // namespace detail


  /**
   *  Bit operations on integers of type G_int
   */
  template <class G_int>
  class bitops {
  public:
    /// The integer type
    typedef G_int value_t;
    
    /// The number of bits in value_t
    enum { bits  = 8 * sizeof(value_t) };

    /// The number of bytes in value_t
    enum { bytes = sizeof(value_t) };

  public:
    
    /// The number of bits set in n
    static inline value_t one_count(value_t n) {
      return detail::bitops<value_t, bits/2>::one_count(n);
    }

    /// The number of bits unset in n
    static inline value_t zero_count(value_t n) {
      return bits - one_count(n);
    }

    /// The reverse of the bits in n
    static inline value_t reverse(value_t n) {
      return detail::bitops<value_t, bits/2>::reverse(n);
    }

    /// The number of low order bits set to zero
    static inline value_t low_order_zero_count(value_t n) {
      return detail::bitops<value_t, bits/2>::low_order_zero_count(n, bits);
    }      

    /// The log2 of n
    static inline value_t log2 (value_t n) {
      SL_REQUIRE("Positive", n > 0);
      return detail::bitops<value_t, bits/2>::log2(n,0);
    }

    /// Is n a power of 2?
    static inline bool is_power2(value_t n) {
      return (n & (n - 1)) == 0;
    }

    /// The log_base of n, requires power of 2 base
    static inline value_t log(value_t base, value_t n) {
      SL_REQUIRE("Good base", is_power2(base)); 
      SL_REQUIRE("Positive", n > 0);
      return log2(n)/log2(base);
    }

    /// The largest power of 2 <= n
    static inline value_t prev_power2(value_t n) {
      return value_t(1<<log2(n));
    }

    /// The smallest power of 2 >= n
    static inline value_t next_power2(value_t n) {
      return (n == 0) ? value_t(1) : (is_power2(n) ? n : prev_power2(2*n));
    }
    
    /// The parity of n
    static inline bool parity(value_t n) {
      return detail::bitops<value_t,bits/2>::parity(n) & 1;
    }

    /// The gray code of n
    static inline value_t gray_from_binary(value_t n) {
      return n ^ (n >> 1);
    }

    /// The binary number corresponding to gray code g
    static inline value_t binary_from_gray(value_t g) {
      return detail::bitops<value_t,bits/2>::binary_from_gray(g);
    }

    /**
     *  Extract count bits from n, starting from bit from. 
     */
    static inline value_t extract(value_t n, 
				  std::size_t from,
				  std::size_t count) {
      static const value_t allones = detail::bitpattern<value_t,0xff,8,bytes>::value;
      return (n >> (bits - from - count)) & ~(allones << count);
    }

    /**
     *  Extract count bits from the bit string starting at p, 
     *  starting from bit from. The routine requires that 
     *  the bit pattern does not straddle element boundaries.
     */
    static inline value_t extract_unstraddled_safe(const value_t *restrict p,
						   std::size_t from,
						   std::size_t count) {
      return extract (* (p + from/bits),
		      from % bits,
		      count);
    }

    /**
     *  Extract count bits from the bit string starting at p, 
     *  starting from bit from.
     */
    static inline value_t extract_straddled_safe(const value_t *restrict p,
						  std::size_t from,
						  std::size_t count) {
      std::size_t count_word_1 = min(count,bits-from%bits);
      std::size_t count_word_2 = count - count_word_1;
           
      return 
	(extract_unstraddled_safe (p, from, count_word_1) << count_word_2) | 
	(extract_unstraddled_safe (p, from+count_word_1, count_word_2));
    }

    /**
     *  Extract count bits from bit vector p, starting at bit position from, 
     *  but not exceeding maxbit. The integer value vector starting at p is 
     *  treated as a long bit vector. This function returns the requested bit 
     *  from this vector. The caller guarantees that [from, from+count) does 
     *  not straddle boundaries of elements of p. from can overflow maxbit, 
     *  in which case nothing is returned. However if from does not exceed 
     *  maxbit, it is assumed that neither does from+count. 
     */
    static inline value_t extract_unstraddled(const value_t *restrict p,
					      std::size_t maxbit,
					      std::size_t from,
					      std::size_t count) {
       return from > maxbit ? 0 : extract_unstraddled_safe (p, from, count);      
    }

    /**
     *  Extract count bits from bit vector p, starting at bit position from, 
     *  but not exceeding maxbit, possibly straddling
     *  element boundaries.
     */
    static inline value_t extract_straddled(const value_t *restrict p,
					    std::size_t maxbit,
					    std::size_t from,
					    std::size_t count) {
      return from > maxbit ? 0 : extract_straddled_safe (p, from, count);
    }
    
  };

  /// Morton bit operations
  template <class G_int, std::size_t G_dimension>
  class morton_bitops: public bitops<G_int> {
  public:
    /// The integer type
    typedef G_int value_t;

    /// The dimension
    enum { dimension = G_dimension };

    /// The number of bits in value_t
    enum { bits  = 8 * sizeof(value_t) };

    /// The number of bytes in value_t
    enum { bytes = sizeof(value_t) };

  public:

    /// The morton code of the d-th coordinated
    static value_t encoded(value_t value, std::size_t d);

    /// The d-th coordinate value from the morton code
    static value_t decoded(value_t code, std::size_t d);
  };

  template <class value_t, std::size_t dimension>
  inline value_t morton_bitops<value_t,dimension>::encoded(value_t value, std::size_t d) {
    value_t result = 0;
    value_t dec_bit = value_t(1);
    value_t enc_bit = value_t(1)<<value_t(d);
    for (std::size_t i = 0; i < bits/dimension; ++i) {
      result |= (value&dec_bit) ? value_t(0) : enc_bit;
      dec_bit <<= value_t(1);
      enc_bit <<= value_t(dimension);
    }
    return result;
  }

  template <class value_t, std::size_t dimension>
  inline value_t morton_bitops<value_t,dimension>::decoded(value_t code, std::size_t d) {
    value_t result = 0;
    value_t dec_bit = value_t(1);
    value_t enc_bit = value_t(1)<<value_t(d);
    for (std::size_t i = 0; i < bits/dimension; ++i) {
      result |= (code&enc_bit) ? value_t(0) : dec_bit;
      dec_bit <<= value_t(1);
      enc_bit <<= value_t(dimension);
    }
    return result;
  }

  /// Morton bit operations, dimension = 2
  template <class G_int>
  class morton_bitops<G_int,2>: public bitops<G_int> {
  public:
    /// The integer type
    typedef G_int value_t;

    /// The dimension
    enum { dimension = 2 };

    /// The number of bits in value_t
    enum { bits  = 8 * sizeof(value_t) };

    /// The number of bytes in value_t
    enum { bytes = sizeof(value_t) };

  public:

    /// The morton code of the d-th coordinated
    static value_t encoded(value_t value, std::size_t d);

    /// The d-th coordinate value from the morton code
    static value_t decoded(value_t code, std::size_t d);
  };

  template <class value_t>
  inline value_t morton_bitops<value_t,2>::encoded(value_t value, std::size_t d) {
    return detail::bitops<value_t,bits/4>::interleaved(value) << d;
  }

  template <class value_t>
  inline value_t morton_bitops<value_t,2>::decoded(value_t code, std::size_t d) {
    return detail::bitops<value_t,bits/4>::deinterleaved(code >> d);
  }

} // namespace sl

#endif


