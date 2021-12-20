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
#ifndef SL_FIXED_RC_INT_CODEC_HPP
#define SL_FIXED_RC_INT_CODEC_HPP

#include <sl/cstdint.hpp>
#include <sl/range_codec.hpp>
#include <sl/bitops.hpp>


namespace sl {

  class rc_adaptive_bit_model {
  protected:
    uint8_t c_;
  public:
    inline rc_adaptive_bit_model(uint8_t c = 128) {
      c_ = c;
    }
    
    inline void reset(uint8_t c = 128) {
      c_ = c;
    }

  public:
    
    inline void encode_bit(range_codec& rc_codec, bool x) {
      rc_codec.encode_bit(&c_, (x?1:0));
    }

    inline bool decode_bit(range_codec& rc_codec) {
      return rc_codec.decode_bit(&c_);
    }
  };

  class rc_adaptive_int_model {
  protected:
    rc_adaptive_bit_model zero_tag_context_;
    rc_adaptive_bit_model sign_context_[32];
    rc_adaptive_bit_model mantissa_context_[32];
    rc_adaptive_bit_model exponent_context_[32];
  public:
    inline rc_adaptive_int_model() {
    }
    inline void reset() {
      zero_tag_context_.reset();
      for (std::size_t i=0; i<32; ++i) {
        sign_context_[i].reset();
        exponent_context_[i].reset();
        mantissa_context_[i].reset();
      }     
    }

    inline void  encode_uint(range_codec& rc_codec, int32_t x) {
      assert(x>=0);
      bool is_zero = (x == 0);
      zero_tag_context_.encode_bit(rc_codec, is_zero);
      if (!is_zero) {
        uint32_t a = x;
        int32_t e = 0; while (a>=(uint32_t(1)<<uint32_t(e+1))) ++e;
        assert(e<31);
        
        // Exponent
        for (int i=0; i<e; i++){
          exponent_context_[i].encode_bit(rc_codec, 1);
        }
        exponent_context_[e].encode_bit(rc_codec, 0);

        // Mantissa
        for (int i=e-1; i>=0; i--) {
          mantissa_context_[i].encode_bit(rc_codec, (a>>i)&1);
        }
      }
    }

    inline int32_t decode_uint(range_codec& rc_codec) {
      bool is_zero = zero_tag_context_.decode_bit(rc_codec);
      if (is_zero) {
        return 0;
      } else {
        // Exponent
        int32_t e = 0;
        while (exponent_context_[e].decode_bit(rc_codec)) {
          ++e;
        }
        // Mantissa
        uint32_t a = uint32_t(1)<<uint32_t(e);
        for (int i=e-1; i>=0; i--) {
          uint32_t bit_i = mantissa_context_[i].decode_bit(rc_codec);
          a |= (bit_i<<uint32_t(i));
        }
        return a;
      }
    }
    
    inline void  encode_int(range_codec& rc_codec, int32_t x) {
#if 0
      // Positive/Negative interleaving
      if (x==0) {
        encode_uint(rc_codec, 0);
      } else if (x>0) {
        encode_uint(rc_codec, x*2-1);
      } else {
        encode_uint(rc_codec, -x*2);
      }
#else
      // Sign bit
      bool is_zero = (x == 0);
      zero_tag_context_.encode_bit(rc_codec, is_zero);
      if (!is_zero) {
        uint32_t a = x > 0 ? x : -x;
        int32_t e = 0; while (a>=(uint32_t(1)<<uint32_t(e+1))) ++e;
        assert(e<32);
        assert(e == bitops<int32_t>::log2(a));
        
        // Exponent
        for (int32_t i=0; i<e; i++){
          exponent_context_[i].encode_bit(rc_codec, 1);
        }
        exponent_context_[e].encode_bit(rc_codec, 0);

        // Mantissa
        for (int32_t i=e-1; i>=0; i--) {
          mantissa_context_[i].encode_bit(rc_codec, (a>>i)&1);
        }

        // Sign
        sign_context_[e].encode_bit(rc_codec, x<0);
      }
#endif
    }

    inline int32_t decode_int(range_codec& rc_codec) {
#if 0
      // Positive/Negative interleaving
      int32_t result = decode_uint(rc_codec);
      if (result==0) {
        // OK
      } else if (result&1) {
        ++result; result>>=1;
      } else {
        result>>=1; result= -result;
      }
      return result;
#else
      // Sign bit
      bool is_zero = zero_tag_context_.decode_bit(rc_codec);
      if (is_zero) {
        return 0;
      } else {
        // Exponent
        int32_t e = 0;
        while (exponent_context_[e].decode_bit(rc_codec)) {
          ++e;
        }
        assert(e<32);

        // Mantissa
        uint32_t a = 1<<e;
        for (int i=e-1; i>=0; i--) {
          uint32_t bit_i = mantissa_context_[i].decode_bit(rc_codec);
          a |= (bit_i<<uint32_t(i));
        }
        // Sign
        bool sign = sign_context_[e].decode_bit(rc_codec);
        return (sign) ? -int32_t(a) : int32_t(a);
      }
#endif
    }
  };

  class rc_adaptive_symbol_model {
  protected:
    uint32_t symbol_count_;
    uint32_t symbol_bit_count_;
    rc_adaptive_bit_model bit_context_[16];
  public:
    inline rc_adaptive_symbol_model() {
      symbol_count_ = 0;
      symbol_bit_count_ = 0;
    }

    inline uint32_t model_symbols() const {
      return symbol_count_;
    }

    inline void set_alphabet(uint32_t number_of_symbols) {
      symbol_count_ = number_of_symbols;
      symbol_bit_count_ = 0;
      if (number_of_symbols > 1) {
	symbol_bit_count_ = 1;
	uint32_t n = 1;
	while (n<number_of_symbols-1) {
	  ++symbol_bit_count_;
	  n*=2;
	}
      }
      assert(symbol_bit_count_ < 16);
    }

    inline void reset() {
      for (std::size_t i=0; i<16; ++i) {
        bit_context_[i].reset();
      }     
    } 
   
    inline void  encode_symbol(range_codec& rc_codec, uint32_t x) {
      assert(x<model_symbols());
      for (uint32_t i=0; i<symbol_bit_count_; ++i) {
	bit_context_[i].encode_bit(rc_codec, (x>>i)&1);
      }
    }

    inline uint32_t decode_symbol(range_codec& rc_codec) {
      uint32_t x = 0;
      for (uint32_t i=0; i<symbol_bit_count_; ++i) {
	uint32_t bit_i = bit_context_[i].decode_bit(rc_codec);
	x |= (bit_i<<i);
      }
      return x;
    }
  };
  
  /**
   *  Simple fixed range coder for 32 bit integers. Integers
   *  are arithmetic coded simply using a single symbol to encode
   *  if a number is 0 and if not by encoding the number using its
   *  exponent, mantissa and sign, with different context for
   *  each encoded bit. This is a variation of the
   *  symbol encoding method used in the FFV1 Video Codec.
   *  It is slower but more effective than the fixed
   *  Huffman rle codec.
   */
  class fixed_rc_int_codec {
  public:
    typedef rc_adaptive_symbol_model symbol_context_t;
    typedef rc_adaptive_bit_model    bit_context_t;
    typedef rc_adaptive_int_model    int_context_t;
  protected:
    bool               is_encoding_;
    bool               is_decoding_;

    range_codec   rc_codec_;
    int_context_t       default_int_context_;
    bit_context_t      default_bit_context_;

  public: // constructor

    // Construct codec
    inline fixed_rc_int_codec() {
      is_encoding_ = 0;
      is_decoding_ = 0;
    }

    // Destruct codec
    inline ~fixed_rc_int_codec() {
    }
    
  public: // buffer

    /// Are we currently encoding a buffer?
    inline bool is_encoding() const {
      return is_encoding_;
    }

    /// Are we currently decoding a buffer?
    inline bool is_decoding() const {
      return is_decoding_;
    }

    /// Set current buffer for encoding/decoding.
    inline void set_buffer(uint8_t* b, std::size_t n) {
      assert(!is_encoding());
      assert(!is_decoding());
      rc_codec_.set_buffer(b, n);
    }

    /// The current encoding/decoding buffer
    inline const uint8_t *buffer() const {
      return rc_codec_.buffer();
    }

    /// The maximum capacity of current buffer
    inline std::size_t buffer_capacity() const {
      return rc_codec_.buffer_capacity();
    }

    /// The number of bits encoded or decoded so far
    inline std::size_t current_bit_count() const {
      return rc_codec_.current_bit_count();
    }

    /// The number of bytes encoded or decoded so far
    inline std::size_t current_byte_count() const {
      return rc_codec_.current_byte_count();
    }
    
  protected: // helpers

    inline void rc_reset_probabilities() {
      default_bit_context_.reset();
      default_int_context_.reset();
    }
    
  public: // encoding

    /**
     *  Start encoding on current buffer. If the buffer is NULL,
     *  the codec just performs a dry run, counting the number
     *  of emitted bits.
     */
    inline void start_encoder() {
      assert(!is_encoding());
      assert(!is_decoding());
      is_encoding_ = true;
      rc_codec_.start_encoder();
      rc_reset_probabilities();
      assert(is_encoding());
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_symbol(symbol_context_t& ctx, uint32_t x) {
      assert(is_encoding());
      ctx.encode_symbol(rc_codec_, x);
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_int(int_context_t& ctx, int32_t x) {
      assert(is_encoding());
      ctx.encode_int(rc_codec_, x);
    }
    
    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_int(int32_t x) {
      assert(is_encoding());
      default_int_context_.encode_int(rc_codec_, x);
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_bit(bit_context_t& ctx, bool x) {
      assert(is_encoding());
      ctx.encode_bit(rc_codec_, x);
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_bit(bool x) {
      assert(is_encoding());
      default_bit_context_.encode_bit(rc_codec_, x);
    }
    
    /**
     *  Stop the encoder, flushing output so as to end on a byte
     *  boundary
     */
    inline void stop_encoder() {
      assert(is_encoding());
      is_encoding_ = false;
      rc_codec_.stop_encoder();
      assert(!is_encoding());
    }
    
  public: // decoding

    /**
     *  Start decoding on current buffer. The buffert cannot be NULL.
     */
    inline void start_decoder() {
      assert(!is_encoding());
      assert(!is_decoding());
      assert(buffer());
      is_decoding_ = true;
      rc_codec_.start_decoder();
      rc_reset_probabilities();
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline uint32_t decode_symbol(symbol_context_t& ctx) {
      assert(is_decoding());
      return ctx.decode_symbol(rc_codec_);
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline int32_t decode_int(int_context_t& ctx) {
      assert(is_decoding());
      return ctx.decode_int(rc_codec_);
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline int32_t decode_int() {
      assert(is_decoding());
      return default_int_context_.decode_int(rc_codec_);
    }
    
    /**
     *  Return next bit in the encoded buffer and advance.
     */
    inline bool decode_bit(bit_context_t& ctx) {
      assert(is_decoding());
      return ctx.decode_bit(rc_codec_);
    }

    /**
     *  Return next bit in the encoded buffer and advance.
     */
    inline bool decode_bit() {
      assert(is_decoding());
      return default_bit_context_.decode_bit(rc_codec_);
    }

    /**
     *  Stop decoding decoding on current buffer. 
     */
    inline void stop_decoder() {
      assert(is_decoding());
      is_decoding_ = false;
      rc_codec_.stop_decoder();
      assert(!is_decoding());
    }
      
  };
} // namespace sl

#endif
