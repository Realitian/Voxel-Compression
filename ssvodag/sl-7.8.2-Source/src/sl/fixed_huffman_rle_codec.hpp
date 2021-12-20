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
#ifndef SL_FIXED_HUFFMAN_RLE_CODEC_HPP
#define SL_FIXED_HUFFMAN_RLE_CODEC_HPP

#include <sl/cstdint.hpp>
#include <cassert>

namespace sl {

  class null_context {
  public:
    inline null_context() {}
    inline void reset() {}
  };

  class huf_adaptive_symbol_context {
  protected:
    uint32_t symbol_count_;
    uint32_t symbol_bit_count_;
  public:
    inline huf_adaptive_symbol_context() {
      symbol_count_ = 0;
      symbol_bit_count_ = 0;
    }
    inline void reset() {}

    inline uint32_t model_symbols() const {
      return symbol_count_;
    }
    inline uint32_t symbol_bit_count() const {
      return symbol_bit_count_;
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
  };
  
  /**
   *  Simple fixed Huffman coder with run length encoding of zeros.
   *  Up to 128 consecutive zero values are encoded in 8 bits, whereas
   *  the nonzero values require 2n+2 bits, n being the number of bits
   *  needed to represent the absolute value of the coded number.
   *  It is sensibly faster than arithmetic coding and works
   *  reasonably well for coding small sparse matrices (e.g., of quantized
   *  wavelet coefficients.). This is a variation of the codec used in
   *  the following paper:
   *
   *  Stefan Guthe, Michael Wand, Julius Gonser, and Wolfgang Straßer,
   *  Interactive Rendering of Large Volume Data Sets
   *  Proc. IEEE Visualization 2002.
   */
  class fixed_huffman_rle_codec {
  public:
    typedef huf_adaptive_symbol_context symbol_context_t;
    typedef null_context                bit_context_t;
    typedef null_context                int_context_t;
  protected:
    bool        is_encoding_;
    bool        is_decoding_;
    uint8_t     zero_count_;
    uint8_t    *buffer_;
    std::size_t buffer_capacity_;
    uint8_t    *byte_pointer_;
    std::size_t bit_index_;

  public: // constructor

    // Construct codec
    inline fixed_huffman_rle_codec() {
      is_encoding_ = 0;
      is_decoding_ = 0;
#if RLE_ZERO_RUN_ENABLED
      zero_count_ = 0;
#endif
      buffer_ = 0;
      buffer_capacity_ = 0;
      byte_pointer_ = 0;
      bit_index_ = 0;
    }

    // Destruct codec
    inline ~fixed_huffman_rle_codec() {
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
      buffer_ = b;
      buffer_capacity_ = n;
    }

    /// The current encoding/decoding buffer
    inline uint8_t *buffer() const {
      return buffer_;
    }

    /// The maximum capacity of current buffer
    inline std::size_t buffer_capacity() const {
      return buffer_capacity_;
    }

    /// The number of bits encoded or decoded so far
    inline std::size_t current_bit_count() const {
      return 8*(byte_pointer_-buffer_) + bit_index_;
    }

    /// The number of bytes encoded or decoded so far
    inline std::size_t current_byte_count() const {
      return (byte_pointer_-buffer_) + (bit_index_ ? 1 : 0);
    }
    
  protected: // helpers

    /// Encode a bit, storing into a buffer if it exists
    inline void output_bit(bool b) {
      assert(is_encoding());
      assert(bit_index_ < 8);
      if (buffer_) {
        assert(byte_pointer_ < buffer_ + buffer_capacity_);
        const uint8_t bit = uint8_t(1<<bit_index_);
        uint8_t byte = *byte_pointer_;
        byte &= ~bit;       // clear
        if (b) byte |= bit; // conditionally set
        *byte_pointer_ = byte;
      }
      ++bit_index_;
      if (bit_index_ == 8) {
        ++byte_pointer_;
        bit_index_ = 0;
      }
    }
    
#if RLE_ZERO_RUN_ENABLED
    /// Terminate current zero run
    inline void output_flush_zero_run() {
      assert(is_encoding());
      if (zero_count_) {
        assert(zero_count_> 0 && zero_count_ <=128);
        output_bit(0);
        int x = zero_count_ - 1;
        output_bit((x >> 0) & 1);
        output_bit((x >> 1) & 1);
        output_bit((x >> 2) & 1);
        output_bit((x >> 3) & 1);
        output_bit((x >> 4) & 1);
        output_bit((x >> 5) & 1);
        output_bit((x >> 6) & 1);
        output_bit((x >> 7) & 1);
        zero_count_ = 0;
      }
    }
#endif
    
    /// Output a zero
    inline void output_zero() {
      assert(is_encoding());
#if RLE_ZERO_RUN_ENABLED
      if (zero_count_) {
        ++zero_count_;
        if (zero_count_ == 128) {
          output_flush_zero_run();
        }
      } else {
        zero_count_ = 1;
      }
#else
      output_bit(0);
#endif
    }

    /// Output a positive number
    inline void output_positive(uint32_t x) {
      assert(is_encoding());
      assert(x > 0 && x < (1 << 30));

#if RLE_ZERO_RUN_ENABLED
      output_flush_zero_run();
#endif
      
      // Signal number
      output_bit(1);
      
      // Write number of digits as N 1 bits followed by 0
      for (int digits = 30; digits > 0; --digits) {
        if (x >= uint32_t(1 << digits)) {
          output_bit(1);
        }
      }
      output_bit(0);
      for (int digits = 29; digits >= 0; --digits) {
        if (x >= uint32_t(2 << digits)) {
          if (x & (1 << digits)) {
            output_bit(1); 
          } else {
            output_bit(0); 
          }
        }
      }
    }

    /// Output a signed number
    inline void output_number(int32_t x) {
      assert(is_encoding());
      if (x<0) {
        output_positive(-x*2);  
      } else if (x>0) {
        output_positive(x*2-1);
      } else {
        output_zero();
      }
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
      byte_pointer_ = buffer_;
      bit_index_ = 0;
#if RLE_ZERO_RUN_ENABLED
      zero_count_ = 0;
#endif
      assert(is_encoding());
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_symbol(symbol_context_t& ctx, int32_t x) {
      assert(is_encoding());
#if RLE_ZERO_RUN_ENABLED
      output_flush_zero_run();
#endif
      const uint32_t nbits = ctx.symbol_bit_count();
      for (uint32_t i=0; i<nbits; ++i) {
	output_bit((x>>i)&1);
      }
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_int(int_context_t& /*ctx*/, int32_t x) {
      assert(is_encoding());
      output_number(x);
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_int(int32_t x) {
      assert(is_encoding());
      output_number(x);
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_bit(bit_context_t& /*ctx*/, bool x) {
      assert(is_encoding());
#if RLE_ZERO_RUN_ENABLED
      output_flush_zero_run();
#endif
      output_bit(x);
    }

    /**
     *  Encode x in current buffer. It the buffer is NULL, the
     *  codec just updates the current bit count.
     */
    inline void encode_bit(bool x) {
      assert(is_encoding());
#if RLE_ZERO_RUN_ENABLED
      output_flush_zero_run();
#endif
      output_bit(x);
    }
    
    /**
     *  Stop the encoder, flushing output so as to end on a byte
     *  boundary
     */
    inline void stop_encoder() {
      assert(is_encoding());
#if RLE_ZERO_RUN_ENABLED
      output_flush_zero_run();
#endif
      while (bit_index_ !=0) output_bit(0); // End on byte boundary
      is_encoding_ = false;
      assert(!is_encoding());
    }

  public: // decoding helpers

    /// Return next bit in buffer and advance.
    inline bool input_bit() {
      assert(is_decoding());
      assert(buffer_);
      assert(byte_pointer_ < buffer_ + buffer_capacity_);
      assert(bit_index_ < 8);
      const bool b = bool((*byte_pointer_) & (uint8_t(1<<bit_index_)!=0));
      ++bit_index_;
      if (bit_index_ == 8) {
        ++byte_pointer_;
        bit_index_ = 0;
      }
      return b;
    }

    /// Return next number in buffer and advance
    inline int32_t input_number() {
      assert(is_decoding());
#if RLE_ZERO_RUN_ENABLED
      if (zero_count_) {
        --zero_count_;
        return 0;
      } else {
#endif
        bool b = input_bit();
        if (b==0) {
#if RLE_ZERO_RUN_ENABLED
          // zero run
          int x = 0;
          b = input_bit(); x |= (b?1:0)<<0;
          b = input_bit(); x |= (b?1:0)<<1;
          b = input_bit(); x |= (b?1:0)<<2;
          b = input_bit(); x |= (b?1:0)<<3;
          b = input_bit(); x |= (b?1:0)<<4;
          b = input_bit(); x |= (b?1:0)<<5;
          b = input_bit(); x |= (b?1:0)<<6;
          b = input_bit(); x |= (b?1:0)<<7;
          zero_count_ = x;
#endif
          return 0;
        } else {
          // nonzero number
          int digits = 0;
          while(input_bit()) {
            ++digits;
          }
          int32_t result = 1;
          for (int i = 0; i < digits; ++i) {
            result = (result<<1) + (input_bit()?1:0);
          }
          // recover sign
          if (result&1) {
            result += 1;
            result >>=1;
          } else {
            result >>=1;
            result = -result;
          }
          return result;
        }
#if RLE_ZERO_RUN_ENABLED
      }
#endif
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
#if RLE_ZERO_RUN_ENABLED
      zero_count_ = 0;
#endif
      byte_pointer_ = buffer_;
      bit_index_ = 0;
      assert(is_decoding());
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline uint32_t decode_symbol(symbol_context_t& ctx) {
      assert(is_decoding());
#if RLE_ZERO_RUN_ENABLED
      assert(zero_count_ == 0);
#endif
      const uint32_t nbits = ctx.symbol_bit_count();
      uint32_t x = 0;
      for (uint32_t i=0; i<nbits; ++i) {
	uint32_t bit_i = input_bit();
	x |= (bit_i<<i);
      }
      return x;
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline int32_t decode_int(int_context_t& /*ctx*/) {
      assert(is_decoding());
      return input_number();
    }

    /**
     *  Return next number in the encoded buffer and advance.
     */
    inline int32_t decode_int() {
      assert(is_decoding());
      return input_number();
    }
    
    /**
     *  Return next bit in the encoded buffer and advance.
     */
    inline bool decode_bit() {
      assert(is_decoding());
#if RLE_ZERO_RUN_ENABLED
      assert(zero_count_ == 0);
#endif
      return input_bit();
    }

    /**
     *  Return next bit in the encoded buffer and advance.
     */
    inline bool decode_bit(bit_context_t& /*ctx*/) {
      assert(is_decoding());
#if RLE_ZERO_RUN_ENABLED
      assert(zero_count_ == 0);
#endif
      return input_bit();
    }

    /**
     *  Stop decoding decoding on current buffer. 
     */
    inline void stop_decoder() {
      assert(is_decoding());
      is_decoding_ = false;
      assert(!is_decoding());
    }
      
  };
} // namespace sl

#endif
