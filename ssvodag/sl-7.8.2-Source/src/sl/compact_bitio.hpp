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
#ifndef COMPACT_BITIO_HPP
#define COMPACT_BITIO_HPP

#include <sl/cstdint.hpp>
#include <cassert>

namespace sl {

  /**
   *  A simple class for writing bits and integers to a memory buffer
   */
  class compact_bitio {
  protected:
    uint8_t*    buffer_;
    std::size_t buffer_capacity_;
    uint8_t*    byte_pointer_;
    uint8_t     bit_index_;
  public:

    compact_bitio() :
        buffer_(0), buffer_capacity_(0), byte_pointer_(0), bit_index_(0) {
    }

    ~compact_bitio() {
    }
    
    inline void start_encoder(void* buf, std::size_t buf_size) {
      buffer_ = static_cast<uint8_t*>(buf);
      buffer_capacity_ = buf_size;
      byte_pointer_ = buffer_;
      bit_index_ = 0;
    }
  
    inline void stop_encoder() {
      while (bit_index_) output_bit(0); // End on byte boundary
    }
    
    inline void start_decoder(const void* buf, std::size_t buf_size) {
      buffer_ = const_cast<uint8_t*>(static_cast<const uint8_t*>(buf));
      buffer_capacity_ = buf_size;
      byte_pointer_ = buffer_;
      bit_index_ = 0;
    }

    inline void stop_decoder() {
      // Nothing to do
    }

    inline std::size_t bit_count() const {
      return 8*(byte_pointer_-buffer_) + bit_index_;
    }
    
    inline std::size_t byte_count() const {
      return (byte_pointer_-buffer_) + (bit_index_ ? 1 : 0);
    }

  public: // Bit i/o
    
    inline void output_bit(bool b) {
      assert(bit_index_ < 8);
      if (buffer_) {
        assert(byte_pointer_ < buffer_ + buffer_capacity_);
        const uint8_t bit = uint8_t(1<<bit_index_);
        uint8_t byte = *byte_pointer_;
        byte &= ~bit;          // clear
        byte |= (b ? bit : 0); // conditionally set
        *byte_pointer_ = byte;
      }
      ++bit_index_;
      bit_index_ &= 7;
      if (bit_index_ == 0) ++byte_pointer_;
      //if (buffer_)  { std::cerr << (b ? '1' : '0'); }
    }

    inline bool input_bit() {
      assert(buffer_);
      assert(byte_pointer_ < buffer_ + buffer_capacity_);
      assert(bit_index_ < 8);
      const bool b = bool((*byte_pointer_) & (uint8_t(1<<bit_index_)!=0));
      ++bit_index_;
      bit_index_ &= 7;
      if (bit_index_ == 0) ++byte_pointer_;
      //if (buffer_) { std::cerr << (b ? '1' : '0'); }
      return b;
    }

  public: // Int i/o
    
    inline void output_byte_coded_uint(uint32_t x) {
      uint32_t bits = 8;
      bool b8 = (x<=255);
      output_bit(b8);
      if (!b8) {
        bits = 16;
        bool b16 = (x<=65535);
        output_bit(b16);
        if (!b16) {
          bits=32;
        }
      }
      for (uint32_t i=0; i<bits; ++i) {
        output_bit((x>>i)&1);
      }
      //if (buffer_) std::cerr << "OUT(" << x << ")" << std::endl;
    }
      
    inline uint32_t input_byte_coded_uint() {
      uint32_t bits=8;
      bool b8 = input_bit();
      if (!b8) {
        bits = 16;
        bool b16 = input_bit();
        if (!b16) {
          bits=32;
        }
      }
      uint32_t x = 0;
      for (uint32_t i=0; i<bits; ++i) {
        x |= (uint32_t(input_bit()?1:0)<<i);
      }
      //if (buffer_) std::cerr << "IN(" << x << ")" << std::endl;
      return x;      
    }

  
    inline void output_byte_coded_int(int32_t x) {
      uint32_t ux = (x<0) ? (((-x)<<1)-1) : x<<1;
      output_byte_coded_uint(ux);
    }
  
    inline int32_t input_byte_coded_int() {
      uint32_t ux = input_byte_coded_uint();
      return ((ux&1)==0) ? int32_t(ux>>1) : -int32_t((ux+1)>>1);
    }

  public: // elias int i/o

    inline void output_elias_positive_uint(uint32_t x) {
      assert(x > 0 && x < (1 << 30));
      
      // Write number of digits as N 1 bits followed by 0
      for (int digits = 30; digits > 0; --digits) {
        if (x >= uint32_t(1 << digits)) {
          output_bit(1);
        }
      }
      output_bit(0);
      
      // Write digits
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

    inline void output_elias_uint(uint32_t value) {
      output_elias_positive_uint(value+1);
    }

    inline void output_elias_int(int32_t value) {
      if (value<0) {
        output_elias_positive_uint(-2*value); // -1 -2 -3 ... -> 2 4 6 8
      } else {
        output_elias_positive_uint(2*value+1); // 0 1 2 ... -> 1 3 5 7
      }
    }

    inline uint32_t input_elias_positive_uint() {
      // nonzero number
      int32_t result = 1;

      int digits = 0;
      while(input_bit()) {
        ++digits;
      }
      for (int i = 0; i < digits; ++i) {
        result = (result<<1) + (input_bit()?1:0);
      }
      return result;
    }

    inline uint32_t input_elias_uint() {
      return input_elias_positive_uint() - 1;
    }

    inline int32_t input_elias_int() {
      uint32_t uvalue = input_elias_positive_uint();
      if (uvalue&1) {
        // 1 3 5 7 ... -> 0 1 2 3 ...
        return int32_t(uvalue-1)/2;
      } else {
        // 2 4 6 8 ... -> -1 -2 -3 -4 ...
        return -int32_t(uvalue)/2;
      }
    }
    
  public: // punctured elias int i/o

    inline void output_punctured_elias_positive_uint(uint32_t value) {
      // All bits are encoded from lsb to msb
      // 0 -> 0
      // 1 -> 10 if inner 1
      // 1 -> 11 if last 1
      assert(value > 0);
      bool end_of_int = false;
      do {
        bool bit = (value & 1);
        value >>= 1;
        if (bit) {
          output_bit(1);
          end_of_int = (value == 0);
          output_bit(end_of_int);
        } else {
          output_bit(0);
        }
      } while (!end_of_int);
    }

    inline void output_punctured_elias_uint(uint32_t value) {
      output_punctured_elias_positive_uint(value+1);
    }

    inline void output_punctured_elias_int(int32_t value) {
      if (value<0) {
        output_punctured_elias_positive_uint(-2*value); // -1 -2 -3 ... -> 2 4 6 8
      } else {
        output_punctured_elias_positive_uint(2*value+1); // 0 1 2 ... -> 1 3 5 7
      }
    }

    inline uint32_t input_punctured_elias_positive_uint() {
      // All bits are encoded from lsb to msb
      // 0 -> 0
      // 1 -> 10 if inner 1
      // 1 -> 11 if last 1
      uint32_t value = 0;
      uint32_t k = 0;
      bool end_of_int = false;
      do {
        bool bit = input_bit();
        if (bit) {
          value |= uint32_t(1<<k);
          end_of_int = input_bit();
        }
      } while (!end_of_int);

      return value;
    }

    inline uint32_t input_punctured_elias_uint() {
      return input_punctured_elias_positive_uint() - 1;
    }

    inline int32_t input_punctured_elias_int() {
      uint32_t uvalue = input_punctured_elias_positive_uint();
      if (uvalue&1) {
        // 1 3 5 7 ... -> 0 1 2 3 ...
        return int32_t(uvalue-1)/2;
      } else {
        // 2 4 6 8 ... -> -1 -2 -3 -4 ...
        return -int32_t(uvalue)/2;
      }
    }
    
  public: // rice int i/o

    inline uint32_t rice_uint_bitcount(uint32_t k, uint32_t value) const {
      uint32_t unarycode=value>>k;
      return unarycode + 1 + k;
    }
    
    inline void output_rice_uint(uint32_t k, uint32_t value) {
      uint32_t unarycode=value>>k;
      while (unarycode--) output_bit(0);
      output_bit(1);
      uint32_t bits=k;					
      while (bits--) { output_bit((value>>bits)&1); }
    }

    inline uint32_t input_rice_uint(uint32_t k) {
      uint32_t value;
      
      uint32_t bit;
      uint32_t unarycode=0;
      while (!(bit=input_bit())) ++unarycode;
      value=unarycode;
      uint32_t bits=k;					
      while (bits--) { value=((value<<1) | (input_bit()?1:0)); }
      return value;
    }

    inline uint32_t rice_k_from_sum_count(uint64_t sum, uint32_t count) {
      // k = min{j:2^j>mean=sum/count}
      uint32_t k=0;
      while ((count<<k)<sum) ++k;
      return k;
    }

    inline uint32_t rice_int_bitcount(uint32_t k, int32_t x) const  {
      uint32_t ux = (x<0) ? ((uint32_t(-x)<<1)-1) : x<<1;
      return rice_uint_bitcount(k, ux);
    }

    inline void output_rice_int(uint32_t k, int32_t x) {
      uint32_t ux = (x<0) ? ((uint32_t(-x)<<1)-1) : x<<1;
      output_rice_uint(k, ux);
    }
  
    inline int32_t input_rice_int(uint32_t k) {
      uint32_t ux = input_rice_uint(k);
      return ((ux&1)==0) ? int32_t(ux>>1) : -int32_t((ux+1)>>1);
    }
    
  };

}

#endif
