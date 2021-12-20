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
#ifndef SL_RANGE_CODEC_HPP
#define SL_RANGE_CODEC_HPP

#include <sl/cstdint.hpp>
#include <cassert>

namespace sl {

  
  /**
   * Range coder based upon
   *    G. N. N. Martin, "Range encoding: an algorithm for removing
   *    redundancy from a digitised message. Proc. Video and Data
   *    Recording Conference, 1979.
   *
   */
  class range_codec {
  protected:
    int low_;
    int range_;
    int outstanding_count_;
    int outstanding_byte_;
    uint8_t zero_state_[256];
    uint8_t  one_state_[256];
    uint8_t *bytestream_start_;
    uint8_t *bytestream_;
    uint8_t *bytestream_end_;
  protected:
    inline void rac_encoder_init(uint8_t *buf, std::size_t buf_size) {
      bytestream_start_= buf; 
      bytestream_= buf;
      bytestream_end_ = buf + buf_size;
      low_= 0;
      range_= 0xFF00;
      outstanding_count_= 0;
      outstanding_byte_= -1;
    }
    
    inline void rac_decoder_init(const uint8_t *buf, std::size_t buf_size) {
      rac_encoder_init(const_cast<uint8_t*>(buf), buf_size);
      low_  =(*bytestream_++)<<8;
      low_ += *bytestream_++;
    }
    
    inline int rac_encoder_terminate() {
      range_=0xFF;
      low_ +=0xFF;
      rac_encoder_renormalize();
      range_=0xFF;
      rac_encoder_renormalize();
      // Encode eof
      *bytestream_ = 0;
      ++bytestream_;
      *bytestream_ = 0;
      ++bytestream_;
      assert(low_   == 0);
      assert(range_ >= 0x100);
      return static_cast<int>(bytestream_ - bytestream_start_);
    }
    
    void rac_build_states(int factor, int max_p) {
      const int64_t one= int64_t(1)<<int64_t(32);
      int64_t p;
      int last_p8, p8;

      for (int i=0; i<256; ++i) {
        zero_state_[i] = 0;
        one_state_[i] = 0;
      }

      last_p8= 0;
      p= one/2;
      for(int i=0; i<128; i++){
        p8= (256*p + one/2) >> 32; //FIXME try without the one
        if(p8 <= last_p8) p8= last_p8+1;
        if(last_p8 && last_p8<256 && p8<=max_p)
          one_state_[last_p8]= p8;
        
        p+= ((one-p)*factor + one/2) >> 32;
        last_p8= p8;
      }

      for(int i=256-max_p; i<=max_p; i++){
        if(one_state_[i]) 
          continue;
        
        p= (i*one + 128) >> 8;
        p+= ((one-p)*factor + one/2) >> 32;
        p8= (256*p + one/2) >> 32; //FIXME try without the one
        if(p8 <= i) p8= i+1;
        if(p8 > max_p) p8= max_p;
        one_state_[i]=     p8;
      }
      
#if 0
      for (int i=max_p+1; i<256; ++i) {
	one_state_[i] = one_state_[max_p];
      }
      for (int i=0; i<256-max_p; ++i) {
	one_state_[i] = one_state_[256-max_p];
      }
#endif
      for(int i=1; i<256; i++) {
        zero_state_[i]= 256-one_state_[256-i];
      }

#if 0
      for(int i=0; i<256; i++) {
	std::cerr << "zero = " << (int)(zero_state_[i]) << " one = " << (int)(one_state_[i]) << std::endl;
      }
#endif
    }
    
    void rac_encoder_renormalize() {
      //FIXME optimize
      while(range_ < 0x100) {
        if (outstanding_byte_ < 0) {
          outstanding_byte_ = low_>>8;
        } else if (low_ <= 0xFF00) {
          *bytestream_++ = outstanding_byte_;
          for(;outstanding_count_; outstanding_count_--)
            *bytestream_++ = 0xFF;
          outstanding_byte_= low_>>8;
        } else if (low_ >= 0x10000) {
          *bytestream_++ = outstanding_byte_ + 1;
          for(;outstanding_count_; outstanding_count_--)
            *bytestream_++ = 0x00;
          outstanding_byte_= (low_>>8) & 0xFF;
        } else {
          outstanding_count_++;
        }
        
        low_ = (low_ & 0xFF)<<8;
        range_ <<= 8;
      }
    }
    
    inline void rac_put_bit(uint8_t * const state, int bit){
      int range1= (range_ * (*state)) >> 8;

      assert(*state);
      assert(range1 < range_);
      assert(range1 > 0);
      if(!bit){
        range_ -= range1;
        *state= zero_state_[*state];
      }else{
        low_ += range_ - range1;
        range_ = range1;
        *state= one_state_[*state];
      }
      rac_encoder_renormalize();
    }

    inline void rac_refill() {
      if (range_ < 0x100){
        range_ <<= 8;
        low_ <<= 8;
        if (bytestream_ < bytestream_end_) {
          low_+= bytestream_[0];
        }
        bytestream_++;
      }
    }

    inline int rac_get_bit(uint8_t * const state) {
      int range1= (range_ * (*state)) >> 8;
      range_ -= range1;

      if (low_ < range_){
        *state= zero_state_[*state];
        rac_refill();
        return 0;
      } else {
        low_ -= range_;
        *state= one_state_[*state];
        range_ = range1;
        rac_refill();
        return 1;
      }
    }

  public:

    inline range_codec() {
      rac_encoder_init(0,0); // clears all variables
      rac_build_states((int)(0.05*(int64_t(1)<<int64_t(32))), 128+64+32+16); // FIXME parameterize!
    }

    inline ~range_codec() {
    }

    // Set current encoding/deconding buffer
    inline void set_buffer(uint8_t* b, std::size_t n) {
      rac_encoder_init(b, n);
    }

    /// The current encoding/decoding buffer
    inline const uint8_t *buffer() const {
      return bytestream_start_;
    }
    
    /// The maximum capacity of current buffer
    inline std::size_t buffer_capacity() const {
      return bytestream_end_ - bytestream_start_;
    }

    /// The number of bytes encoded or decoded so far
    inline std::size_t current_byte_count() const {
      return bytestream_ - bytestream_start_;
    }

    /// The number of bits encoded or decoded so far
    inline std::size_t current_bit_count() const {
      return 8*current_byte_count();
    }
    
    /**
     *  Start encoding on current buffer.
     */
    inline void start_encoder() {
      rac_encoder_init(bytestream_start_, bytestream_end_-bytestream_start_);
    }

    /**
     * Encode a bit using the given context.
     */
    inline void encode_bit(uint8_t * const state, int bit) {
      rac_put_bit(state, bit);
    }

    /**
     *  Stop the encoder, flushing output so as to end on a byte
     *  boundary
     */
    inline void stop_encoder() {
      rac_encoder_terminate();
    }
    
    /**
     *  Start encoding on current buffer.
     */
    inline void start_decoder() {
      rac_decoder_init(bytestream_start_, bytestream_end_-bytestream_start_);
    }

    /**
     * Encode a bit using the given context.
     */
    inline bool decode_bit(uint8_t * const state) {
      return rac_get_bit(state)!=0;
    }

    /**
     *  Stop the decoder
     */
    inline void stop_decoder() {
      /** NOTHING TO DO **/
    }
    
  }; // class range_codec
  
} // namespace sl

#endif
