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
#ifndef SL_FIXED_AC_INT_CODEC_HPP
#define SL_FIXED_AC_INT_CODEC_HPP

#include <sl/cstdint.hpp>
#include <sl/arithmetic_codec.hpp>
#include <sl/bitops.hpp>

namespace sl {

  /**
   *  Simple arithmetic coding context based on
   *  a fixed huffmann coding of ints. Each bit
   *  is modeled using an adaptive bit model.
   */
  class ac_adaptive_int_model {
  protected:
    adaptive_bit_model zero_tag_context_;
    adaptive_bit_model sign_context_[32];
    adaptive_bit_model mantissa_context_[32];
    adaptive_bit_model exponent_context_[32];
  public:
    inline ac_adaptive_int_model() {
    }
    inline void reset() {
      zero_tag_context_.reset();
      for (std::size_t i=0; i<32; ++i) {
        sign_context_[i].reset();
        mantissa_context_[i].reset();
        exponent_context_[i].reset();
      }     
    }

  public:
    
    inline void encode_int(arithmetic_codec& ac_codec, int32_t x) {
      bool is_zero = (x == 0);
      ac_codec.encode(is_zero, zero_tag_context_);
      if (!is_zero) {
        int32_t a = x > 0 ? x : -x;
        int32_t e = bitops<int32_t>::log2(a);

        // Exponent
        for (int i=0; i<e; i++){
            ac_codec.encode(1, exponent_context_[i]);
        }
        ac_codec.encode(0, exponent_context_[e]);

        // Mantissa
        for (int i=e-1; i>=0; i--) {
          ac_codec.encode((a>>i)&1, mantissa_context_[i]);
        }

        // Sign
        ac_codec.encode(x<0, sign_context_[e]);
      }
    }

    inline int32_t decode_int(arithmetic_codec& ac_codec) {
      bool is_zero = ac_codec.decode(zero_tag_context_)!=0;
      if (is_zero) {
        return 0;
      } else {
        // Exponent
        int32_t e = 0;
        while (ac_codec.decode(exponent_context_[e])) {
          ++e;
        }
        // Mantissa
        int32_t a = 1<<e;
        for (int i=e-1; i>=0; i--) {
          int32_t bit_i = ac_codec.decode(mantissa_context_[i]);
          a |= (bit_i<<i);
        }
        // Sign
        bool sign = ac_codec.decode(sign_context_[e])!=0;
        if (sign) {
          return -a;
        } else {
          return a;
        }
      }
    }
    
    inline void encode_uint(arithmetic_codec& ac_codec, int32_t a) {
      assert(a>=0);
      bool is_zero = (a == 0);
      ac_codec.encode(is_zero, zero_tag_context_);
      if (!is_zero) {
        int32_t e = bitops<int32_t>::log2(a);

        // Exponent
        for (int i=0; i<e; i++){
            ac_codec.encode(1, exponent_context_[i]);
        }
        ac_codec.encode(0, exponent_context_[e]);

        // Mantissa
        for (int i=e-1; i>=0; i--) {
          ac_codec.encode((a>>i)&1, mantissa_context_[i]);
        }
      }
    }

    inline int32_t decode_uint(arithmetic_codec& ac_codec) {
      bool is_zero = ac_codec.decode(zero_tag_context_)!=0;
      if (is_zero) {
        return 0;
      } else {
        // Exponent
        int32_t e = 0;
        while (ac_codec.decode(exponent_context_[e])) {
          ++e;
        }
        // Mantissa
        int32_t a = 1<<e;
        for (int i=e-1; i>=0; i--) {
          int32_t bit_i = ac_codec.decode(mantissa_context_[i]);
          a |= (bit_i<<i);
        }
        return a;
      }
    }
    
  };

  
  /**
   *  Simple fixed arithmetic coder for 32 bit integers. Integers
   *  are arithmetic coded simply using a single symbol to encode
   *  if a number is 0 and if not by encoding the number using its
   *  exponent, mantissa and sign, with different context for
   *  each encoded bit. This is a variation of the
   *  symbol encoding method used in the FFV1 Video Codec.
   *  It is slower but more effective than the fixed
   *  Huffman rle codec.
   */
  class fixed_ac_int_codec {
  public:
    typedef adaptive_data_model    symbol_context_t;
    typedef adaptive_bit_model     bit_context_t;
    typedef ac_adaptive_int_model  int_context_t;
  protected:
    bool               is_encoding_;
    bool               is_decoding_;

    arithmetic_codec   ac_codec_;
    int_context_t      default_int_context_;
    bit_context_t      default_bit_context_;
    
  public: // constructor

    // Construct codec
    inline fixed_ac_int_codec() {
      is_encoding_ = 0;
      is_decoding_ = 0;
    }

    // Destruct codec
    inline ~fixed_ac_int_codec() {
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
      ac_codec_.set_buffer(static_cast<unsigned>(n), b);
    }

    /// The current encoding/decoding buffer
    inline const uint8_t *buffer() const {
      return ac_codec_.buffer();
    }

    /// The maximum capacity of current buffer
    inline std::size_t buffer_capacity() const {
      return ac_codec_.buffer_size();
    }

    /// The number of bits encoded or decoded so far
    inline std::size_t current_bit_count() const {
      return ac_codec_.current_encode_byte_count()*8;
    }

    /// The number of bytes encoded or decoded so far
    inline std::size_t current_byte_count() const {
      return ac_codec_.current_encode_byte_count();
    }
    
  protected: // helpers

    inline void ac_reset_probabilities() {
      default_int_context_.reset();
      default_bit_context_.reset();
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
      ac_codec_.start_encoder();
      ac_reset_probabilities();
      assert(is_encoding());
    }

    /**
     *  Encode x in current buffer using the given context
     */
    inline void encode_symbol(symbol_context_t& ctx, uint32_t x) {
      assert(is_encoding());
      ac_codec_.encode(x, ctx);
    }

    /**
     *  Encode x in current buffer using the given context context.
     */
    inline void encode_int(int_context_t& ctx, int32_t x) {
      assert(is_encoding());
      ctx.encode_int(ac_codec_, x);
    }

      /**
     *  Encode x in current buffer using the 'default' int context.
     */
    inline void encode_int(int32_t x) {
      assert(is_encoding());
      default_int_context_.encode_int(ac_codec_, x);
    }

    /**
     *  Encode x in current buffer using the given context context.
     */
    inline void encode_bit(bit_context_t& ctx, bool x) {
      assert(is_encoding());
      ac_codec_.encode(x, ctx);
    }

      /**
     *  Encode x in current buffer using the 'default' int context.
     */
    inline void encode_bit(bool x) {
      assert(is_encoding());
      ac_codec_.encode(x, default_bit_context_);
    }
    
    /**
     *  Stop the encoder, flushing output so as to end on a byte
     *  boundary
     */
    inline void stop_encoder() {
      assert(is_encoding());
      is_encoding_ = false;
      ac_codec_.stop_encoder();
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
      ac_codec_.start_decoder();
      ac_reset_probabilities();
    }

    /**
     *  Encode x in current buffer using the given context
     */
    inline uint32_t decode_symbol(symbol_context_t& ctx) {
      assert(is_decoding());
      return ac_codec_.decode(ctx);
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline int32_t decode_int(int_context_t& ctx) {
      assert(is_decoding());
      return ctx.decode_int(ac_codec_);
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline int32_t decode_int() {
      assert(is_decoding());
      return default_int_context_.decode_int(ac_codec_);
    }
    
    /**
     *  Return next bit in the encoded buffer and advance.
     */
    inline bool decode_bit(bit_context_t& ctx) {
      assert(is_decoding());
      return ac_codec_.decode(ctx)!=0;
    }

    /**
     *  Return next bit in the encoded buffer and advance.
     */
    inline bool decode_bit() {
      assert(is_decoding());
      return ac_codec_.decode(default_bit_context_)!=0;
    }

    /**
     *  Stop decoding decoding on current buffer. 
     */
    inline void stop_decoder() {
      assert(is_decoding());
      is_decoding_ = false;
      ac_codec_.stop_decoder();
      assert(!is_decoding());
    }
      
  };
} // namespace sl

#endif
