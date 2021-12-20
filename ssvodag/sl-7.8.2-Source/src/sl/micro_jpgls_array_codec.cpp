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
#include <sl/micro_jpgls_array_codec.hpp>
#include <sl/float_cast.hpp>
#include <sl/dense_array.hpp>
#include <sl/indexed_functions.hpp>
#include <sl/generative_types.hpp>
#include <sl/assert.hpp>
#include <sl/bitops.hpp>
#include <sl/numeric_traits.hpp>
#include <limits>
#include <algorithm>
#include <cassert>

namespace sl {
  
  // ---------------------------------------------------------------------------------
  // BIT I/O
  // ---------------------------------------------------------------------------------

  inline void micro_jpgls_array_codec::bitio_start_encoder(void* buf, std::size_t buf_size) {
    bitio_buffer_ = static_cast<uint8_t*>(buf);
    bitio_buffer_capacity_ = buf_size;
    bitio_byte_pointer_ = bitio_buffer_;
    bitio_bit_index_ = 0;
  }

  inline void micro_jpgls_array_codec::bitio_start_decoder(const void* buf, std::size_t buf_size) {
    bitio_buffer_ = const_cast<uint8_t*>(static_cast<const uint8_t*>(buf));
    bitio_buffer_capacity_ = buf_size;
    bitio_byte_pointer_ = bitio_buffer_;
    bitio_bit_index_ = 0;
  }
    
  inline void micro_jpgls_array_codec::bitio_output_bit(bool b) {
    assert(bitio_bit_index_ < 8);
    if (bitio_buffer_) {
      assert(bitio_byte_pointer_ < bitio_buffer_ + bitio_buffer_capacity_);
      const uint8_t bit = uint8_t(1<<bitio_bit_index_);
      uint8_t byte = *bitio_byte_pointer_;
      byte &= ~bit;          // clear
      byte |= (b ? bit : 0); // conditionally set
      *bitio_byte_pointer_ = byte;
    }
    ++bitio_bit_index_;
    bitio_bit_index_ &= 7;
    if (bitio_bit_index_ == 0) ++bitio_byte_pointer_;
    //if (bitio_buffer_)  { std::cerr << (b ? '1' : '0'); }
  }

  inline bool micro_jpgls_array_codec::bitio_input_bit() {
    assert(bitio_buffer_);
    assert(bitio_byte_pointer_ < bitio_buffer_ + bitio_buffer_capacity_);
    assert(bitio_bit_index_ < 8);
    const bool b = (bool)(((*bitio_byte_pointer_) & uint8_t(1<<bitio_bit_index_))!=0);
    ++bitio_bit_index_;
    bitio_bit_index_ &= 7;
    if (bitio_bit_index_ == 0) ++bitio_byte_pointer_;
    //if (bitio_buffer_) { std::cerr << (b ? '1' : '0'); }
    return b;
  }

  inline void micro_jpgls_array_codec::bitio_output_uint(uint32_t x) {
    uint32_t bits = 8;
    bool b8 = (x<=255);
    bitio_output_bit(b8);
    if (!b8) {
      bits = 16;
      bool b16 = (x<=65535);
      bitio_output_bit(b16);
      if (!b16) {
        bits=32;
      }
    }
    for (uint32_t i=0; i<bits; ++i) {
      bitio_output_bit((x>>i)&1);
    }
    //if (bitio_buffer_) std::cerr << "OUT(" << x << ")" << std::endl;
  }
  
  inline void micro_jpgls_array_codec::bitio_output_int(int32_t x) {
    uint32_t ux = (x<0) ? (-x*2-1) : x*2;
    bitio_output_uint(ux);
  }

  inline uint32_t micro_jpgls_array_codec::bitio_input_uint() {
    uint32_t bits=8;
    bool b8 = bitio_input_bit();
    if (!b8) {
      bits = 16;
      bool b16 = bitio_input_bit();
      if (!b16) {
        bits=32;
      }
    }
    uint32_t x = 0;
    for (uint32_t i=0; i<bits; ++i) {
      x |= (uint32_t(bitio_input_bit()?1:0)<<i);
    }
    //if (bitio_buffer_) std::cerr << "IN(" << x << ")" << std::endl;
    return x;      
  }

    
  inline int32_t micro_jpgls_array_codec::bitio_input_int() {
    uint32_t ux = bitio_input_uint();
    return (ux%2==0) ? (ux/2) : -int32_t((ux+1)/2);
  }
  
  inline void micro_jpgls_array_codec::bitio_stop_encoder() {
    while (bitio_bit_index_) bitio_output_bit(0); // End on byte boundary
  }

  inline void micro_jpgls_array_codec::bitio_stop_decoder() {
    // Nothing to do
  }

  inline std::size_t micro_jpgls_array_codec::bitio_bit_count() const {
    return 8*(bitio_byte_pointer_-bitio_buffer_) + bitio_bit_index_;
  }
  
  inline std::size_t micro_jpgls_array_codec::bitio_byte_count() const {
    return (bitio_byte_pointer_-bitio_buffer_) + (bitio_bit_index_ ? 1 : 0);
  }
  
  inline void micro_jpgls_array_codec::bitio_output_limited_golomb(uint32_t k,
                                                                   uint32_t glimit,
                                                                   uint32_t qbpp,
                                                                   uint32_t value) {
    assert(glimit > qbpp+1);
    uint32_t limit=glimit-qbpp-1;
    
    // A.5.3 Mapped-error encoding
    uint32_t unarycode=value>>k;				
    //if (bitio_buffer_) std::cerr << "UC=" << unarycode << ": ";
    if (unarycode < limit) {
      while (unarycode--) bitio_output_bit(0);
      bitio_output_bit(1);
      uint32_t bits=k;					
      while (bits--) { bitio_output_bit((value>>bits)&1); } 	
    } else {
      while (limit--) bitio_output_bit(0);
      bitio_output_bit(1);					
      uint32_t bits=qbpp;					
      while (bits--) { bitio_output_bit(((value-1)>>bits)&1); } 	
    }
    //if (bitio_buffer_) std::cerr << "OUT_G_" << k << "_" << glimit << "_" << qbpp << "(" << value << ")" << std::endl;
  }

  inline uint32_t micro_jpgls_array_codec::bitio_input_limited_golomb(uint32_t k,
                                                                      uint32_t glimit,
                                                                      uint32_t qbpp) {
    assert(glimit > qbpp+1);
    uint32_t limit=glimit-qbpp-1;
    
    uint32_t value;
    
    uint32_t bit;
    uint32_t unarycode=0;
    while (!(bit=bitio_input_bit())) ++unarycode;
    //if (bitio_buffer_) std::cerr << "DECODE: UC=" << unarycode << ": ";
    
    if (unarycode < limit) {
      value=unarycode;
      uint32_t bits=k;					
      while (bits--) { value=((value<<1) | (bitio_input_bit()?1:0)); }
    } else {
      value=0; // no contribution from unary code ... whole value is next
      uint32_t bits=qbpp;					
      while (bits--) { value=((value<<1) | (bitio_input_bit()?1:0)); }
      value+=1; // correct for limited case
    }
    //if (bitio_buffer_) std::cerr << "IN_G_" << k << "_" << glimit << "_" << qbpp << "(" << value << ")" << std::endl;
    return value;
  }
  
  // ---------------------------------------------------------------------------------
  // array access implementation
  // ---------------------------------------------------------------------------------
  
  template <class value_t>
  inline value_t micro_jpgls_array_codec::top_item(const dense_array<value_t,2,void>& array,
                                                   uint32_t i,
                                                   uint32_t j) {
    return (i==0) ? value_t(0) : array(i-1,j);
  }

  template <class value_t>
  inline value_t micro_jpgls_array_codec::left_item(const dense_array<value_t,2,void>& array,
                                                    uint32_t i,
                                                    uint32_t j) {
    return
      (j==0)
      ?
      ((i==0) ? value_t(0) :  array(i-1,j))
      :
      array(i,j-1);
  }

  template <class value_t>
  inline value_t micro_jpgls_array_codec::right_item(const dense_array<value_t,2,void>& array,
                                                            uint32_t i,
                                                            uint32_t j) {
    return array(i,std::min(std::size_t(j+1), array.extent()[1]-1));
  }

  template <class value_t>
  inline value_t micro_jpgls_array_codec::top_left_item(const dense_array<value_t,2,void>& array,
                                                               uint32_t i,
                                                               uint32_t j) {
    return
      (i==0)
      ?
      value_t(0) // FIXME
      :
      left_item(array, i-1, j);
  }

  template <class value_t>
  inline value_t micro_jpgls_array_codec::top_right_item(const dense_array<value_t,2,void>& array,
                                                                uint32_t i,
                                                                uint32_t j) {
    return
      (i==0)
      ?
      value_t(0) // FIXME
      :
      right_item(array, i-1, j);
  }

  // ---------------------------------------------------------------------------------
  // JLS implementation
  // ---------------------------------------------------------------------------------
  
  uint32_t micro_jpgls_array_codec::jls_j_[32] = {
    0, 0, 0, 0,
    1, 1, 1, 1,
    2, 2, 2, 2,
    3, 3, 3, 3,
    4, 4, 5, 5,
    6, 6, 7, 7,
    8, 9,10,11,
   12,13,14,15 };
  
  uint32_t micro_jpgls_array_codec::jls_jn_[32] = {
    (1<< 0), (1<< 0), (1<< 0), (1<< 0),
    (1<< 1), (1<< 1), (1<< 1), (1<< 1),
    (1<< 2), (1<< 2), (1<< 2), (1<< 2),
    (1<< 3), (1<< 3), (1<< 3), (1<< 3),
    (1<< 4), (1<< 4), (1<< 5), (1<< 5),
    (1<< 6), (1<< 6), (1<< 7), (1<< 7),
    (1<< 8), (1<< 9), (1<<10), (1<<11),
    (1<<12), (1<<13), (1<<14), (1<<15) };
     
  void micro_jpgls_array_codec::jls_init(int32_t in_min, int32_t in_max, int32_t eps) { 
    // A.2.1.1
    jls_minval_ = in_min;
    jls_maxval_ = in_max-in_min;
    jls_near_   = eps;
    jls_near_times2_plus1_ = 2*jls_near_+1;
    jls_range_  = jls_near_ ? int32_t(1+(jls_maxval_+2*jls_near_)/(jls_near_times2_plus1_)) : int32_t(jls_maxval_+1);
    jls_qbpp_   = 0; while ((1<<jls_qbpp_)<=jls_range_) ++jls_qbpp_;
    jls_bpp_    = 2; while ((1<<jls_bpp_)<=jls_maxval_+1) ++jls_bpp_;
    jls_limit_  = 2*(jls_bpp_+std::max(int32_t(8),jls_bpp_)); // FIXME ?
#if 0
    std::cerr << "jls_minval : " << jls_minval_ << std::endl;
    std::cerr << "jls_maxval : " << jls_maxval_ << std::endl;
    std::cerr << "jls_near   : " << jls_near_ << std::endl;
    std::cerr << "jls_range  : " << jls_range_ << std::endl;
    std::cerr << "jls_qbpp   : " << jls_qbpp_ << std::endl;
    std::cerr << "jls_bpp    : " << jls_bpp_ << std::endl;
    std::cerr << "jls_limit  : " << jls_limit_ << std::endl;
#endif 
    // A.2.1.2
    int32_t a_init = std::max(int32_t(2), int32_t((jls_range_+(1<<5))/(1<<6)));
    for (uint32_t i=0; i<jls_context_count; ++i) {
      jls_n_[i] = 1;
      jls_b_[i] = 0;
      jls_c_[i] = 0;
      jls_a_[i] = a_init;
    }
    
    // C
    jls_reset_ = 64;
    jls_maxc_ = std::max(2, int32_t((1<<jls_bpp_)/2));
    jls_minc_ = -jls_maxc_ -1;
  }

  template <class int_t>
  void micro_jpgls_array_codec::jls_init(const dense_array<int_t,2,void>& in, int32_t eps) {
    // Compute range
    int32_t in_min = in(0,0);
    int32_t in_max = in_min;
    for (std::size_t i=0; i<in.extent()[0]; ++i) {
      for (std::size_t j=0; j<in.extent()[1]; ++j) {
        int32_t in_ij = int32_t(in(i,j));
        in_min = std::min(in_min, in_ij);
        in_max = std::max(in_max, in_ij);
      }
    }
    jls_init(in_min, in_max, eps);
  }

  inline int32_t micro_jpgls_array_codec::jls_quantize(int32_t x, int32_t eps, int32_t eps_times2_plus1) {
    return 
      eps ?
      ((x>=0) ? (x+eps)/(eps_times2_plus1) : -(eps-x)/(eps_times2_plus1))
      :
      (x);
  }

  inline int32_t micro_jpgls_array_codec::jls_dequantize(int32_t x, int32_t eps_times2_plus1) {
    return x*eps_times2_plus1;
  }
  
  inline int32_t micro_jpgls_array_codec::jls_quantize(int32_t x) const {
    return jls_quantize(x, jls_near_, jls_near_times2_plus1_);
  }
    
  inline int32_t micro_jpgls_array_codec::jls_dequantize(int32_t x) const {
    return jls_dequantize(x, jls_near_times2_plus1_);
  }
  
  inline int32_t micro_jpgls_array_codec::jls_predictor(int32_t a,
                                                        int32_t b,
                                                        int32_t c) {
    int32_t max_b_a = std::max(b, a);
    int32_t min_b_a = std::min(b, a);
    
    if (c >= max_b_a) {
      return min_b_a;
    } else if (c <= min_b_a) {
      return max_b_a;
    } else {        
      return a+b-c;
    }
  }
    
  inline int32_t micro_jpgls_array_codec::jls_context(int32_t a,
                                                      int32_t b,
                                                      int32_t c,
                                                      int32_t d,
                                                      int32_t eps) {
    int32_t g1 = d-b;
    int32_t g2 = b-c;
    int32_t g3 = c-a;
    
    int32_t c1 = (g1<-eps)?-1:(g1>eps)?1:0;
    int32_t c2 = (g2<-eps)?-1:(g2>eps)?1:0;
    int32_t c3 = (g3<-eps)?-1:(g3>eps)?1:0;
    
    int32_t sign = 1;
    if (c1<0) {
      sign = -1;
      c1 = -c1;
      c2 = -c2;
      c3 = -c3;
    } else if (c1==0) {
      if (c2<0) {
        sign = -1;
        c2 = -c2;
        c3 = -c3;
      } else if (c2==0) {
        if (c3<0) {
          sign = -1;
          c3 = -c3;
        }
      }
    }

    // T = 1 -> ((2*T+1)^3+1)/2  = 14 contexts
    // Q1: 0 Q2: 0     Q3  0..1   => 1x1x2: [0..1]
    // Q1: 0 Q2: 1     Q3 -1/0/1  => 1x1x3: [2..4]
    // Q1: 1 Q2:-1/0/1 Q3 -1/0/1  => 1x3x3: [5..13]     
    int32_t context_idx;
    if (c1==0) {
      if (c2==0) {
        context_idx = c3; // [0..1]
      } else {
        context_idx = 2+(c3+1); // [2..4]
      }
    } else {
      // Q1: 1 Q2:-1/0/1 Q3 -1/0/1  => 1x3x3: [5..13]     
      context_idx = 5+(c2+1)+(c3+1)*3;
    }
    assert(context_idx < jls_context_count);
    return sign * context_idx;
  }
  
  inline void micro_jpgls_array_codec::jls_update(uint32_t q, int32_t errval) {
    jls_a_[q] += ((errval >= 0) ? (errval) : (-errval));
    jls_b_[q] += jls_dequantize(errval);
    if (jls_n_[q] == jls_reset_) {
      jls_a_[q] >>= 1;
      jls_b_[q] >>= 1;
      jls_n_[q] >>= 1;
    }
    jls_n_[q] += 1;
    
    if (jls_b_[q] <= -jls_n_[q]) {
      jls_b_[q] += jls_n_[q];
      if (jls_b_[q] <= -jls_n_[q]) {
        jls_b_[q] = -jls_n_[q]+1;
      }
      if (jls_c_[q] > jls_minc_) {
        jls_c_[q] -= 1;
      }
    } else if (jls_b_[q]>0) {
      jls_b_[q] -= jls_n_[q];
      if (jls_b_[q] > 0) {
        jls_b_[q] = 0;
      }
      if (jls_c_[q] < jls_maxc_) {
        jls_c_[q] += 1;
      }
    }
  }

  template <class int_t>
  void micro_jpgls_array_codec::jls_encode(const dense_array<int_t,2,void>& in,
                                           int32_t eps,
                                           void *buf,
                                           std::size_t buf_size,
                                           std::size_t *actual_size,
                                           double *actual_error,
                                           error_kind_t error_kind) {
    // Encode array shape
    const uint32_t h = in.extent()[0];
    const uint32_t w = in.extent()[1];
    
    // Encode array data
    double e2 = 0;
    double eamax = 0;
    
    bitio_start_encoder(buf, buf_size);
    {
      bool is_square = h==w;
      bitio_output_bit(is_square);
      bitio_output_uint(h);
      if (!is_square) bitio_output_uint(w);
      dense_array<int_t,2,void> array_prime(h,w);
      
      jls_init(in, eps);
      bitio_output_int(jls_minval_);
      bitio_output_int(jls_maxval_);
      bitio_output_uint(jls_near_);

      const int32_t nearval = jls_near_;
      const int32_t nearval_times2_plus1 = nearval*2+1;
      const int32_t offset = jls_minval_;

      uint32_t rle_index = 0;
      for (uint32_t i=0; i<h; ++i) {
        int32_t ra = int32_t(left_item(array_prime, i,0))-offset;
        int32_t rb = int32_t(top_item(array_prime, i,0))-offset;
        int32_t rc = int32_t(top_left_item(array_prime, i, 0))-offset;
        for (uint32_t j=0; i<h && j<w; ++j) {
          int32_t  rd = int32_t(top_right_item(array_prime, i, j))-offset;
          int32_t  context_ij = jls_context(ra, rb, rc, rd, nearval);

          // Run mode processing
          if ((i>0) && (j>0) && (context_ij == 0)) { // FIXME
            // Detect run length
            uint32_t rle_count = 0;
            uint32_t run_i_end = i;
            uint32_t run_j_end = j;
            int32_t x_ij = int32_t(in(run_i_end,run_j_end))-offset;
            while ((run_i_end<h && run_j_end<w) &&
                   (abs(x_ij-ra) <= nearval)) {
              ++rle_count;
              ++run_j_end; if (run_j_end==w) { run_j_end=0; ++run_i_end; }
              if (run_i_end<h) x_ij = int32_t(in(run_i_end,run_j_end))-offset;
            }
            // Encode run length
            uint32_t encoded_rle_count = rle_count;
            while (encoded_rle_count >= jls_jn_[rle_index]) {
              bitio_output_bit(1);
              encoded_rle_count -= jls_jn_[rle_index];
              if (rle_index<31) ++rle_index;
            }
            bitio_output_bit(0);
            for (uint32_t k=0; k<jls_j_[rle_index]; ++k) {
              bitio_output_bit((encoded_rle_count>>k) & 1);
            }
            if (rle_index>0) --rle_index;
            
            // Reconstruct
            //if (bitio_buffer_) std::cerr << "EPS=" << nearval << " => ENCODED RL(" << i << " " << j << "): " << rle_count << std::endl;
            for (uint32_t k=0; k<rle_count; ++k) {
              int32_t reconstructed_ij = ra;
              array_prime(i,j) = int_t(reconstructed_ij+offset);
              //if (bitio_buffer_) std::cerr << "  ENC(" << i << " " << j << "): RL=" << ra << std::endl;
              // Update error
              double e = double(array_prime(i,j))-double(in(i,j));
              e2 += e*e;
              eamax = std::max(abs(e), eamax);
              ++j; if (j==w) { j=0; ++i; }
            }
            
            // Update variables
            if (rle_count && (i<h)) {
              ra = int32_t(left_item(array_prime, i,j))-offset;
              rb = int32_t(top_item(array_prime, i,j))-offset;
              rc = int32_t(top_left_item(array_prime, i, j))-offset;
              rd = int32_t(top_right_item(array_prime, i, j))-offset;
              context_ij     = jls_context(ra, rb, rc, rd, nearval);
            }
          }

          // Regular mode processing or run interruption processing
          if (i<h && j<w) {
            int32_t  sign_ij        = (context_ij<0)?-1:1;
            int32_t  context_idx_ij = sign_ij*context_ij;
              
            int32_t  predicted_ij = jls_predictor(ra, rb, rc);
            predicted_ij += sign_ij * jls_c_[context_ij];
            predicted_ij = median(predicted_ij, int32_t(0), int32_t(jls_maxval_));
            
            // A5.1 Golomb coding variable computation
            int32_t golomb_k = 0;
            int32_t nq = jls_n_[context_idx_ij];
            int32_t aq = jls_a_[context_idx_ij];
            int32_t bq = jls_b_[context_idx_ij];
            while ((nq<<golomb_k)<aq) ++golomb_k;
            
            int32_t  delta_ij = sign_ij*((int32_t(in(i,j))-offset)-predicted_ij);
            delta_ij = jls_quantize(delta_ij, nearval, nearval_times2_plus1);
            
            int32_t reconstructed_ij = predicted_ij+sign_ij*jls_dequantize(delta_ij, nearval_times2_plus1);
            reconstructed_ij = median(reconstructed_ij, int32_t(0), int32_t(jls_maxval_));
            array_prime(i,j) = int_t(reconstructed_ij+offset);
              
            // A4.5 Modulo reduction of prediction error
            if (delta_ij<0) delta_ij += jls_range_;
            if (delta_ij>=((jls_range_+1)/2)) delta_ij -= jls_range_;
            
            // A5.2 Error mapping
            int32_t m_delta_ij;
            if ((nearval==0) && (golomb_k==0) && (2*bq<=-nq)) {
              m_delta_ij = (delta_ij>=0) ? (2*delta_ij+1) : -2*(delta_ij+1);
            } else {
              m_delta_ij = (delta_ij>=0) ? (2*delta_ij) : (-2*delta_ij-1);
            }
            bitio_output_limited_golomb(golomb_k, jls_limit_, jls_qbpp_, m_delta_ij);
            
            //if (bitio_buffer_) std::cerr << "ENC(" << i << " " << j << "): ctx = " << context_ij << " delta_ij = " << delta_ij << " -> m_delta_ij = " << m_delta_ij << std::endl;
            
            // A6 Update variables
            jls_update(context_idx_ij, delta_ij);
            
            // Update error
            double e = double(array_prime(i,j))-double(in(i,j));
            e2 += e*e;
            eamax = std::max(abs(e), eamax);
            
            // Advance
            rc = rb;
            rb = rd;
            ra = int32_t(array_prime(i,j))-offset;
          }
        }
      }
    }
    bitio_stop_encoder();
    
    // Return size and error
    //if (bitio_buffer_) std::cerr << "ENCODED => " << bitio_byte_count() << " bytes" << std::endl;
    if (actual_size)  *actual_size = bitio_byte_count();
    if (actual_error) *actual_error = ((error_kind == Error_kind_amax) ? eamax : std::sqrt(e2/double(h*w)));
  }
  
      
  template <class int_t>
  double micro_jpgls_array_codec::jls_encoding_error(const dense_array<int_t,2,void>& in,
                                                     int32_t eps,
                                                     error_kind_t error_kind) {
    std::size_t actual_size;
    double      actual_error;
    jls_encode(in, eps, 0, 0, &actual_size, &actual_error, error_kind);
    return actual_error;
  }
  
  template <class int_t>
  void micro_jpgls_array_codec::jls_decode(dense_array<int_t,2,void>& array_prime,
                                           const void *buf,
                                           std::size_t buf_size) {
    bitio_start_decoder(buf, buf_size);
    {
      const bool is_square = bitio_input_bit();
      const uint32_t h = bitio_input_uint();
      const uint32_t w = (is_square) ? h : bitio_input_uint();
      array_prime.resize(subscript_t(h,w));
        
      const int32_t  jls_minval = bitio_input_int();
      const int32_t  jls_maxval = bitio_input_int();
      const int32_t  jls_near   = bitio_input_uint();
      jls_init(jls_minval, jls_maxval+jls_minval, jls_near);

      const int32_t nearval = jls_near_;
      const int32_t nearval_times2_plus1 = nearval*2+1;
      const int32_t offset = jls_minval_;

      uint32_t rle_index = 0; // Index for run mode order
      for (uint32_t i=0; i<h; ++i) {
        int32_t ra = int32_t(left_item(array_prime, i,0))-offset;
        int32_t rb = int32_t(top_item(array_prime, i,0))-offset;
        int32_t rc = int32_t(top_left_item(array_prime, i, 0))-offset;
        for (uint32_t j=0; i<h && j<w; ++j) {
          int32_t  rd = int32_t(top_right_item(array_prime, i, j))-offset;
          int32_t  context_ij     = jls_context(ra, rb, rc, rd, nearval);
          if ((i>0) && (j>0) && (context_ij == 0)) { // FIXME
            // Run mode processing
              
            // Decode run length
            uint32_t rle_count = 0;
            while (bitio_input_bit()) {
              rle_count += jls_jn_[rle_index];
              if (rle_index<31) ++rle_index;
            }
            for (uint32_t k=0; k<jls_j_[rle_index]; ++k) {
              rle_count += (uint32_t(bitio_input_bit() ? 1 : 0) << k);
            }
            if (rle_index>0) --rle_index;
            
            // Reconstruct
            //if (bitio_buffer_) std::cerr << "EPS=" << nearval << " => DECODED RL(" << i << " " << j << "): " << rle_count << std::endl;
            for (uint32_t k=0; k<rle_count; ++k) {
              int32_t reconstructed_ij = ra;
              array_prime(i,j) = int_t(reconstructed_ij+offset);
              //if (bitio_buffer_) std::cerr << "  DEC(" << i << " " << j << "): RL=" << ra << std::endl;
              ++j; if (j==w) { j=0; ++i; }
            }
            
            // Update variables
            if (rle_count && (i<h)) {
              ra = int32_t(left_item(array_prime, i,j))-offset;
              rb = int32_t(top_item(array_prime, i,j))-offset;
              rc = int32_t(top_left_item(array_prime, i, j))-offset;
              rd = int32_t(top_right_item(array_prime, i, j))-offset;
              context_ij = jls_context(ra, rb, rc, rd, nearval);
            }
          }
          
          // Regular mode processing or run interruption processing
          if (i<h && j<w) {
            int32_t  sign_ij        = (context_ij<0)?-1:1;
            int32_t  context_idx_ij = sign_ij*context_ij;
              
            int32_t  predicted_ij = jls_predictor(ra, rb, rc);
            predicted_ij += sign_ij * jls_c_[context_ij];
            predicted_ij = median(predicted_ij, int32_t(0), int32_t(jls_maxval_));

            // A5.1 Golomb coding variable computation
            int32_t golomb_k = 0;
            int32_t nq = jls_n_[context_idx_ij];
            int32_t aq = jls_a_[context_idx_ij];
            int32_t bq = jls_b_[context_idx_ij];
            while ((nq<<golomb_k)<aq) ++golomb_k;

            int32_t m_delta_ij = bitio_input_limited_golomb(golomb_k, jls_limit_, jls_qbpp_);
            
            // A5.2 Error unmapping
            int32_t delta_ij;
            if ((nearval==0) && (golomb_k==0) && (2*bq<=-nq)) {
              delta_ij = ((m_delta_ij%2)!=0) ? (m_delta_ij-1)/2 : -m_delta_ij/2-1;
            } else {
              delta_ij = ((m_delta_ij%2)==0) ? (m_delta_ij/2) : -(m_delta_ij+1)/2;
            }
            //if (bitio_buffer_) std::cerr << "DEC(" << i << " " << j << "): ctx = " << context_ij << " m_delta_ij = " << m_delta_ij << " -> delta_ij = " << delta_ij << std::endl;
            
            int32_t reconstructed_ij = predicted_ij+sign_ij*jls_dequantize(delta_ij, nearval_times2_plus1);
            if (reconstructed_ij<-nearval) {
              reconstructed_ij += jls_dequantize(jls_range_, nearval_times2_plus1);
            } else if (reconstructed_ij > jls_maxval_+nearval) {
              reconstructed_ij -= jls_dequantize(jls_range_, nearval_times2_plus1);
            }
            reconstructed_ij = median(reconstructed_ij, int32_t(0), int32_t(jls_maxval_));
            array_prime(i,j) = int_t(reconstructed_ij+offset);
            
            // A6 Update variables
            jls_update(context_idx_ij, delta_ij);
            
            // Advance
            rc = rb;
            rb = rd;
            ra = int32_t(array_prime(i,j))-offset;
          }
        }
      }
    }
    bitio_stop_decoder();
    //if (bitio_buffer_) std::cerr << "DECODED => " << bitio_byte_count() << " bytes" << std::endl;
  }

  //----------------------------------------------------------------------------
  // General codec entry points
  //-----------------------------------------------------------------------------

  template <class int_t>
  uint32_t micro_jpgls_array_codec::threshold_from_target_error(const dense_array<int_t,2,void>& array,
                                                                double target_error,
                                                                error_kind_t error_kind) {
    // Check lossless compression
    if (target_error < 1e-6) return 0;

    // Search quantization space
    uint32_t steps = 0;
    
    // Finest possible value
    uint32_t threshold_fine  = uint32_t(target_error); // FIXME 
    double   error_fine = jls_encoding_error(array, threshold_fine, error_kind);
    ++steps;
    //std::cerr << steps << ": " << "T=" << threshold_fine << " => E= " << error_fine << std::endl;
    if (error_fine >=target_error) return threshold_fine;
      
    // Coarsest possible value
    uint32_t threshold_coarse = uint32_t(amax(array));
    if (error_kind == Error_kind_amax) threshold_coarse=std::max(threshold_coarse, 1+threshold_fine);
    
    //std::cerr << steps << ": " << "T=" << threshold_coarse << std::endl;
    if (threshold_coarse <= threshold_fine) return threshold_coarse; // null matrix
    double   error_coarse = jls_encoding_error(array, threshold_coarse, error_kind);
    ++steps;
    //std::cerr << steps << ": " << "T=" << threshold_coarse << " => E= " << error_coarse << std::endl;
    if (error_coarse<=target_error) return threshold_coarse;
    
    // Bracket threshold range with binary search
    assert(error_fine <= target_error);
    assert(target_error <= error_coarse);
    while (threshold_coarse>threshold_fine+1) {
      uint32_t threshold_mid = (threshold_fine+threshold_coarse)/2;
      double   error_mid = jls_encoding_error(array, threshold_mid, error_kind);
      ++steps;
      //std::cerr << steps << ": " << "T=" << threshold_mid << " => E= " << error_mid << std::endl;
      if (error_mid<=target_error) {
        threshold_fine = threshold_mid;
        error_fine = error_mid;
      } else {
        threshold_coarse = threshold_mid;
        error_coarse = error_mid;
      }
    }
    // Conservative
    if (error_coarse <= 1.1*target_error) {
      return threshold_coarse;
    } else {
      return threshold_fine;
    }
  }

  template <class int_t>
  void micro_jpgls_array_codec::int_compress(const dense_array<int_t,2,void>& array,
                                             std::size_t target_size,
                                             double target_error,
                                             void* buf,
                                             std::size_t buf_size,
                                             std::size_t* actual_size, 
                                             double *actual_error,
                                             error_kind_t error_kind) {
    // Guess threshold from target error
    uint32_t threshold_fine = threshold_from_target_error(array, target_error, error_kind);
    
    // Compress and check size
    std::vector<uint8_t> compression_buffer(4*8*array.count());
    jls_encode(array, threshold_fine, &compression_buffer[0], buf_size, actual_size, actual_error, error_kind);
    if (*actual_size <= target_size) {
      // Converged on error
    } else {
      // Buffer overflow!!
      std::cerr << "BUFFER OVERFLOW!!" << std::endl;
      
      uint32_t threshold_coarse = std::max(threshold_fine, uint32_t(30));
      if (threshold_coarse>threshold_fine) {
        std::size_t size_fine = *actual_size;
        std::size_t size_coarse;
        jls_encode(array, threshold_coarse,  &compression_buffer[0], buf_size, &size_coarse, NULL, error_kind);
        if (size_coarse>target_size) {
          // Buffer overflow!!
          SL_FAIL("Codec buffer overflow");
        } else {
          // Binary search: 
          while (threshold_coarse>threshold_fine+1) {
            uint32_t threshold_mid = (threshold_fine+threshold_coarse)/2;
            std::size_t size_mid;
            jls_encode(array, threshold_fine, &compression_buffer[0], buf_size, &size_mid, NULL, error_kind);
            if (size_mid<=target_size) {
              threshold_coarse = threshold_mid;
              size_coarse = size_mid;
            } else {
              threshold_fine = threshold_mid;
              size_fine = size_mid;
            }
          }
          // Choose largest size passing threshold
          if (size_fine<=target_size) {
            jls_encode(array, threshold_fine, &compression_buffer[0], buf_size, actual_size, actual_error, error_kind);
          } else {
            jls_encode(array, threshold_coarse, &compression_buffer[0], buf_size, actual_size, actual_error, error_kind);
          }
        }
      }
    }
    
    // Here, we have the result in the codec buffer, copy to output buffer
    if (*actual_size > buf_size) SL_FAIL("Codec buffer overflow");
    
    uint8_t* obuf = static_cast<uint8_t*>(buf);
    for (std::size_t i=0; i<*actual_size; ++i) {
      obuf[i] = compression_buffer[i];
    }
  }
  
  //----------------------------------------------------------------------------
  // Specific coded entrypoints, specialized per type
  //-----------------------------------------------------------------------------
    
  /**
   *  Decompress buffer and store result into array.
   */
  template <class int_t>
  void micro_jpgls_array_codec::int_decompress(dense_array<int_t,2,void>& array,
                                               const void* buf,
                                               std::size_t buf_size) {
    assert(buf);
    assert(buf_size >= 2);
    
    jls_decode(array, buf, buf_size);
  }
  
  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::int8_compress(const int8_array2_t& array,
                                              std::size_t target_size,
                                              double target_error,
                                              void* buf,
                                              std::size_t buf_size,
                                              std::size_t* actual_size, 
                                              double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }
    
  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::int8_compress_amax(const int8_array2_t& array,
                                                   std::size_t target_size,
                                                   double target_error,
                                                   void* buf,
                                                   std::size_t buf_size,
                                                   std::size_t* actual_size, 
                                                   double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::int8_decompress(int8_array2_t& array,
                                                const void* buf,
                                                std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::int16_compress(const int16_array2_t& array,
                                               std::size_t target_size,
                                               double target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t* actual_size, 
                                               double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::int16_compress_amax(const int16_array2_t& array,
                                                    std::size_t target_size,
                                                    double target_error,
                                                    void* buf,
                                                    std::size_t buf_size,
                                                    std::size_t* actual_size, 
                                                    double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::int16_decompress(int16_array2_t& array,
                                                 const void* buf,
                                                 std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::int32_compress(const int32_array2_t& array,
                                               std::size_t target_size,
                                               double target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t* actual_size, 
                                               double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::int32_compress_amax(const int32_array2_t& array,
                                                    std::size_t target_size,
                                                    double target_error,
                                                    void* buf,
                                                    std::size_t buf_size,
                                                    std::size_t* actual_size, 
                                                    double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::int32_decompress(int32_array2_t& array,
                                                 const void* buf,
                                                 std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::uint8_compress(const uint8_array2_t& array,
                                               std::size_t target_size,
                                               double target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t* actual_size, 
                                               double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::uint8_compress_amax(const uint8_array2_t& array,
                                                    std::size_t target_size,
                                                    double target_error,
                                                    void* buf,
                                                    std::size_t buf_size,
                                                    std::size_t* actual_size, 
                                                    double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::uint8_decompress(uint8_array2_t& array,
                                                 const void* buf,
                                                 std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::uint16_compress(const uint16_array2_t& array,
                                                std::size_t target_size,
                                                double target_error,
                                                void* buf,
                                                std::size_t buf_size,
                                                std::size_t* actual_size, 
                                                double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::uint16_compress_amax(const uint16_array2_t& array,
                                                     std::size_t target_size,
                                                     double target_error,
                                                     void* buf,
                                                     std::size_t buf_size,
                                                     std::size_t* actual_size, 
                                                     double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::uint16_decompress(uint16_array2_t& array,
                                                  const void* buf,
                                                  std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::uint32_compress(const uint32_array2_t& array,
                                                std::size_t target_size,
                                                double target_error,
                                                void* buf,
                                                std::size_t buf_size,
                                                std::size_t* actual_size, 
                                                double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void micro_jpgls_array_codec::uint32_compress_amax(const uint32_array2_t& array,
                                                     std::size_t target_size,
                                                     double target_error,
                                                     void* buf,
                                                     std::size_t buf_size,
                                                     std::size_t* actual_size, 
                                                     double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::uint32_decompress(uint32_array2_t& array,
                                                  const void* buf,
                                                  std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::float_compress(const float_array2_t& array,
                                               std::size_t  target_size,
                                               double       target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t* actual_size,
                                               double      *actual_error,
                                               error_kind_t error_kind) {
    assert(buf);
    assert(buf_size >= 8);
      
    enum { Int32_quantize_bits = 20 };
    const float Int32_quantize_scale =(float)(1<<Int32_quantize_bits);
      
    // Quantize to int32    
    const uint8_t h = array.extent()[0];
    const uint8_t w = array.extent()[1];
    float mx = float(amax(array));
    float qthr = std::max((mx/Int32_quantize_scale),1e-6f);
    if (target_error) {
      // This is here just to speed-up bisection
      qthr = std::max(qthr, float(target_error/1024.0f));
    }
    float inv_qthr = 1.0f/qthr;
    int32_array2_t iarray(h,w);
    for (uint32_t i=0; i<h; ++i) {
      for (uint32_t j=0; j<w; ++j) {
        const float x = array(i,j)*inv_qthr;
#if 0
        iarray(i,j) = int32_t(x<0.0f ? x-0.5f : x+0.5f); // round to nearest integer
#else
        iarray(i,j) = int32_t(fast_round_to_integer(x));
#endif
      }
    }
    //std::cerr << "ENCODED: QTHR: " << qthr << " 1/QTHR=" << inv_qthr << std::endl;
    // Encode
    uint8_t* obuf = static_cast<uint8_t*>(buf);
    obuf += be_encode<float>(qthr, obuf); buf_size-=sizeof(float);
    int_compress(iarray,
                 std::min(target_size,buf_size),
                 target_error*inv_qthr,
                 obuf,
                 buf_size,
                 actual_size,
                 actual_error,
                 error_kind);
    *actual_size += sizeof(float);
    *actual_error = *actual_error * qthr;
  }
      
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::float_compress(const float_array2_t& array,
                                               std::size_t  target_size,
                                               double       target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t* actual_size,
                                               double      *actual_error) {
    float_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::float_compress_amax(const float_array2_t& array,
                                                    std::size_t  target_size,
                                                    double       target_error,
                                                    void* buf,
                                                    std::size_t buf_size,
                                                    std::size_t* actual_size,
                                                    double      *actual_error) {
    float_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
    
  /**
   *  Decompress buffer and store result into array.
   */
  void micro_jpgls_array_codec::float_decompress(float_array2_t& array,
                                                 const void* buf,
                                                 std::size_t buf_size) {
    assert(buf);
    assert(buf_size >= 2);
      
    // Decode
    const uint8_t* ibuf = static_cast<const uint8_t*>(buf);
    float qthr = be_decode<float>(ibuf); ibuf += sizeof(float);
    int32_array2_t iarray;
    int_decompress(iarray, ibuf, buf_size-sizeof(float));
      
    //std::cerr << "DECODED: QTHR: " << qthr << std::endl;
      
    // Dequantize from int32    
    const uint8_t h = iarray.extent()[0];
    const uint8_t w = iarray.extent()[1];
    array.resize(subscript_t(h,w));
    for (uint32_t i=0; i<h; ++i) {
      for (uint32_t j=0; j<w; ++j) {
        array(i,j) = float(iarray(i,j))*qthr;
      }
    }
  }

}
