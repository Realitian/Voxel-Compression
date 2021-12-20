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
#ifndef SL_EMBEDDED_ZEROTREE_CODEC_HPP
#define SL_EMBEDDED_ZEROTREE_CODEC_HPP

#include <sl/dense_array.hpp>
#include <sl/indexed_functions.hpp>
#include <sl/bitops.hpp>
#include <sl/arithmetic_codec.hpp>
#include <sl/fixed_ac_int_codec.hpp>
#include <sl/fixed_rc_int_codec.hpp>
#include <sl/encdec.hpp>
#include <queue>
#include <vector>
#include <cassert>

namespace sl {

  /**
   *  Embededded zerotree codec for signed integer matrices.
   *
   *  Compression technique based on "Embedded Image Coding Using Zerotrees of Wavelet Coefficients"
   *  by Jerome M. Shapiro, IEEE Transactions on Signal Processing, Vol.41, No.12,
   *  December 1993: 3445-3462.
   *
   *  Zerotree detection based on "A Fast Technique for Identifying Zerotrees in the EZW Algorithm"
   *  by Jerome M. Shapiro, IEEE International Conference on Acoustics, Speech, and Signal Processing, 1996.
   *  ICASSP-96. Conference Proceedings, Volume 3, 7-10 May 1996:1455 - 1458
   *
   *  Compression to targer RMS based on "Limiting the Distortion of a Wavelet
   *  Image Codec", Joonas Lehtinen, Acta Cybernetica, Vol 14, 1999: 341-356.
   */
  template <class INT_T>
  class embedded_zerotree_codec {
  public:
    typedef INT_T                              int_t;
    typedef SL_SUMTYPENAME(int_t)              energy_t;
    typedef dense_array<int_t, 2, void>        int_matrix_t;
    typedef typename int_matrix_t::subscript_t subscript_t;
  protected:
    
    // Code alphabet
    typedef enum {
      EZW_ZEROTREE_ROOT  = 0, /* binary 00 */
      EZW_POSITIVE       = 1, /* binary 01 */
      EZW_ISOLATED_ZERO  = 2, /* binary 10 */
      EZW_NEGATIVE       = 3  /* binary 11 */
    } ezw_code_t;

    typedef std::vector<int_t>         int_vector_t;
    typedef std::vector<subscript_t>   subscript_vector_t;
    
  protected:
    std::size_t           bitio_current_buffer_size_;

    typedef fixed_rc_int_codec                 entropy_codec_t;
    typedef entropy_codec_t::symbol_context_t  symbol_context_t;
    typedef entropy_codec_t::bit_context_t     bit_context_t;
    typedef entropy_codec_t::int_context_t     int_context_t;

    symbol_context_t      dominant_data_model_[4];
    bit_context_t         subordinate_bit_model_;
    entropy_codec_t       entropy_codec_;
    std::vector<uint8_t>  entropy_codec_buffer_;

    energy_t              current_encoding_square_error_;
    std::size_t           encoded_symbols_;
    std::size_t           decoded_symbols_;
    
    static inline int_t iabs(int_t x) { return x>0?x:-x; }
    static inline energy_t isqr(energy_t x) { return x*x; }
    
  public:

    embedded_zerotree_codec() : encoded_symbols_(0), decoded_symbols_(0) {
      bitio_current_buffer_size_ = 0;
      dominant_data_model_[0].set_alphabet(4);
      dominant_data_model_[1].set_alphabet(4);
      dominant_data_model_[2].set_alphabet(4);
      dominant_data_model_[3].set_alphabet(4);
    }

    ~embedded_zerotree_codec() {
    }

  protected: // Scanning queue helpers
    
    typedef std::vector< std::pair<subscript_t, bool> > scanning_queue_t;
    scanning_queue_t      scanning_queue_;
    std::size_t           scanning_queue_capacity_;
    std::size_t           scanning_queue_size_;
    std::size_t           scanning_queue_head_;
    std::size_t           scanning_queue_tail_;
    
    inline bool scanning_queue_empty() const {
      return scanning_queue_size_ == 0;
    }

    inline const std::pair<subscript_t, bool>& scanning_queue_front() const {
      return scanning_queue_[scanning_queue_head_];
    }

    inline void scanning_queue_pop() {
      assert(scanning_queue_size_ > 0);
      ++scanning_queue_head_; if (scanning_queue_head_ == scanning_queue_capacity_) scanning_queue_head_ = 0;
      --scanning_queue_size_;
    }

    inline void scanning_queue_push(const std::pair<subscript_t, bool>& x) {
      assert(scanning_queue_size_ < scanning_queue_capacity_);

      scanning_queue_[scanning_queue_tail_] = x;
      ++scanning_queue_tail_; if (scanning_queue_tail_ == scanning_queue_capacity_) scanning_queue_tail_ = 0;
      ++scanning_queue_size_;
    }

    inline void scanning_queue_init(std::size_t h, std::size_t w) {
      scanning_queue_capacity_ = 1+2*h*w;
      scanning_queue_size_     = 0;
      scanning_queue_head_     = 0;
      scanning_queue_tail_     = 0;
      if (scanning_queue_.size() < scanning_queue_capacity_) scanning_queue_.resize(scanning_queue_capacity_);
      scanning_queue_push(std::make_pair(subscript_t(0,0), false));
    }

  protected: // Compressed header

    class compressed_data_header {
    protected:
      bool        encode_extent_;
      std::size_t h_;
      std::size_t w_;
      uint8_t     log2_threshold_;
      std::size_t symbol_count_;

    protected:
      
      static std::size_t encode(uint8_t* buf, std::size_t x) {
        std::size_t header_bytes = 0;
        do {
          uint8_t y = int(x & 0x7F); x >>= 7;
          if (x) y |= 0x80;
          buf[header_bytes] = y;
          ++header_bytes;
        } while (x);
        return header_bytes;
      }

      static std::size_t decode(const uint8_t* buf, std::size_t &x) {
        std::size_t header_bytes = 0;
        std::size_t shift = 0; x = 0;
        uint8_t y;
        do {
          y = buf[header_bytes]; ++header_bytes;
          x |= std::size_t(y & 0x7F) << shift;
          shift += 7;
        } while (y & 0x80);
        return header_bytes;
      }
      
    public:

      inline bool encode_extent() const { return encode_extent_; }
      inline std::size_t h() const { return h_; }
      inline std::size_t w() const { return w_; }
      inline uint8_t log2_threshold() const { return log2_threshold_; }
      inline std::size_t symbol_count() const { return symbol_count_; }
      
      compressed_data_header(std::size_t h,
                             std::size_t w,
                             uint8_t log2_threshold,
                             std::size_t symbol_count,
                             bool encode_extent) :
          encode_extent_(encode_extent),
          h_(h),
          w_(w),
          log2_threshold_(log2_threshold),
          symbol_count_(symbol_count) {
        assert(log2_threshold <= 32);
      }
      
      compressed_data_header(const void* header_buffer,
                             std::size_t h = std::size_t(-1),
                             std::size_t w = std::size_t(-1)) {
        const bool needs_extent = (h == std::size_t(-1)) || (w == std::size_t(-1));

        const uint8_t* buf = (const uint8_t*)header_buffer;
        uint8_t hdr = be_decode<uint8_t>(buf); buf += 1;
        const bool square_matrix = (hdr & 0x80)!=0;
        encode_extent_           = (hdr & 0x40)!=0;
        log2_threshold_          = (hdr & (~(0x80|0x40)));
        // Sanity check
        if (encode_extent_ != needs_extent) SL_FAIL("Decode/encode mismatch");
        buf += decode(buf, symbol_count_);
        if (needs_extent) {
          buf += decode(buf, h_); w_ = h_;
          if (!square_matrix) buf += decode(buf, w_);
        } else {
          h_ = h;
          w_ = w;
        }
      }

      std::size_t encode(void* header_buffer) const {
        uint8_t* buf = (uint8_t*)header_buffer;
        const bool square_matrix = h() == w();
        uint8_t hdr = uint8_t(log2_threshold());
        if (square_matrix)   hdr |= 0x80;
        if (encode_extent()) hdr |= 0x40; // just for checking
        buf += be_encode(hdr, buf); 
        buf += encode(buf, symbol_count());
        if (encode_extent()) {
          buf += encode(buf, h());
          if (!square_matrix) buf += encode(buf, w());
        }
        return buf-(uint8_t*)header_buffer;
      }

      inline std::size_t compressed_size() const {
        uint8_t fake_buffer[128];
        return encode(fake_buffer);
      }
    };
      
  protected: // backend encoder helpers
    
    inline void bitio_reset() {
      dominant_data_model_[0].reset();
      dominant_data_model_[1].reset();
      dominant_data_model_[2].reset();
      dominant_data_model_[3].reset();
      subordinate_bit_model_.reset();
    }

    inline void bitio_start_encoding(std::size_t max_code_bytes) {
      bitio_current_buffer_size_ = max_code_bytes;
      std::size_t min_ac_buffer_size = std::max(bitops<uint32_t>::next_power2(uint32_t(max_code_bytes+16)), uint32_t(1024)); // To avoid overflows
      if (entropy_codec_.buffer_capacity() < min_ac_buffer_size) {
	entropy_codec_buffer_.resize(min_ac_buffer_size);
        entropy_codec_.set_buffer(&(entropy_codec_buffer_[0]), min_ac_buffer_size);
      }
      entropy_codec_.start_encoder();
      bitio_reset();
      encoded_symbols_ = 0;
    }

    // Stop encoding anc copy data into buffer
    inline std::size_t bitio_stop_encoding(const compressed_data_header& hdr,
                                           uint8_t* buf) {
      // copy back into current_buffer
      entropy_codec_.stop_encoder();
      std::size_t code_byte_count = entropy_codec_.current_byte_count();      
      if (hdr.symbol_count() == 0) code_byte_count = 0; 
      const std::size_t header_byte_count = hdr.compressed_size();
      const std::size_t total_byte_count = header_byte_count + code_byte_count;
      assert(total_byte_count<bitio_current_buffer_size_);
      // Write header
      buf += hdr.encode(buf);
      // Write compressed symbols
      for (std::size_t i=0; i<code_byte_count; ++i) {
        buf[i] = entropy_codec_.buffer()[i];
      }
      return total_byte_count;
    }
      
    std::size_t bitio_encode_byte_count() const {
      return entropy_codec_.current_byte_count();
    }

    inline void bitio_encode(ezw_code_t code,
                             bool parent_is_significant,
                             bool previous_is_significant) {
      ++encoded_symbols_;
      //entropy_codec_.encode_symbol(dominant_data_model_[0], code);
      entropy_codec_.encode_symbol(dominant_data_model_[(int(parent_is_significant&1)<<1)|
							(int(previous_is_significant&1))],
				   code);
      // FIXME { char str[4] = { 't', 'p', 'z', 'n' }; std::cerr << str[code]; }
    }

    inline void bitio_encode(bool bit) {
      ++encoded_symbols_;
      entropy_codec_.encode_bit(subordinate_bit_model_, bit);
      // FIXME std::cerr << code ? '1' : '0';
    }

  protected: // decoding
    
    inline compressed_data_header bitio_start_decoding(uint8_t* code_buffer,
                                                       std::size_t code_bytes,
                                                       std::size_t h,
                                                       std::size_t w) {
      compressed_data_header hdr(code_buffer, h, w);
      encoded_symbols_ = hdr.symbol_count();
      decoded_symbols_ = 0;
      if (encoded_symbols_) {
        const std::size_t hdr_size = hdr.compressed_size();
        
        bitio_current_buffer_size_ = code_bytes-hdr_size;
	std::size_t min_ac_buffer_size = std::max(bitops<uint32_t>::next_power2(uint32_t(code_bytes+16)), uint32_t(1024)); // To avoid overflows
        if (entropy_codec_.buffer_capacity() < min_ac_buffer_size) {
	  entropy_codec_buffer_.resize(min_ac_buffer_size);
	  entropy_codec_.set_buffer(&(entropy_codec_buffer_[0]), min_ac_buffer_size);
        }
        for (std::size_t i=0; i<bitio_current_buffer_size_; ++i) {
          entropy_codec_buffer_[i] = code_buffer[i+hdr_size];
        }
        entropy_codec_.start_decoder();
        bitio_reset();
      }
      return hdr;
    }

    inline void bitio_stop_decoding() {
      if (encoded_symbols_) entropy_codec_.stop_decoder();
    }

    inline ezw_code_t bitio_decode(bool parent_is_significant,
                                   bool previous_is_significant) {
      ++decoded_symbols_;
      //ezw_code_t code = ezw_code_t(entropy_codec_.decode_symbol(dominant_data_model_[0]));
      ezw_code_t code = ezw_code_t(entropy_codec_.decode_symbol(dominant_data_model_[(int(parent_is_significant&1)<<1)|
										     (int(previous_is_significant&1))]));
      // FIXME { char str[4] = { 't', 'p', 'z', 'n' }; std::cerr << str[code]; }
      return code;
    }

    inline bool bitio_decode() {
      ++decoded_symbols_;
      bool bit = entropy_codec_.decode_bit(subordinate_bit_model_);
      // FIXME std::cerr << code ? '1' : '0';
      return bit;
    }
    
  public:

    /**
     *  Compress array and encode it in given buffer. Compression terminates when
     *  the size equals or exceeds target size or when the root mean square error
     *  falls below the target value. The input array is assumed to encode a
     *  hierarchical transformation using a pyramidal structure with root at (0,0)
     *  and children of (i,j) at (2i..2i+1)-(2j-2j+1). If encode_extent is true,
     *  the encoded data is totally self contained, otherwise the user should provide
     *  the matrix extent when decoding.
     */
    void compress(const int_matrix_t& array,
                  std::size_t target_size,
                  double target_rms,
                  void* buf,
                  std::size_t buf_size,
                  std::size_t *actual_size,
                  double *actual_rms,
                  bool encode_extent = true) {
      assert(buf);
      assert(buf_size);
      assert(buf_size >= target_size);
      assert(actual_size);
      assert(actual_rms);

      const std::size_t max_extent_bytes = encode_extent ? 2 + 2 : 0; 
      const std::size_t max_header_bytes = max_extent_bytes + 1 + 3; // extent()[0] | extent()[1] | log2(threshold) | encoded_symbols;
      const std::size_t max_code_bytes   = buf_size - max_header_bytes; 
      uint8_t* data_buffer = (uint8_t*)buf;

      // Initialize error estimate
      current_encoding_square_error_ = energy(array);
      
      // Find starting threshold
      int_t array_amax = amax(array);
      uint8_t log2_threshold = 1;
      int_t threshold = 1 << log2_threshold;
      while (array_amax >= 2*threshold) {
        ++log2_threshold;
        threshold = 1 << log2_threshold;
      }

      // Compute zerotree map
      int_matrix_t zerotree_map(array.extent());
      zerotree_map_in(zerotree_map, array);
      
      // Encode
      bitio_start_encoding(max_code_bytes);
      {
        const std::size_t target_code_size = target_size > max_header_bytes ? target_size-max_header_bytes : 0;
        const energy_t target_square_error = energy_t(target_rms*target_rms*double(array.count()));
        
        int_vector_t subordinate_list; subordinate_list.reserve(array.count());
        while (threshold!=0) {
          ezw_encode_pass(array,threshold,
                          zerotree_map,
                          target_code_size,target_square_error,subordinate_list);
          threshold >>= 1;
        }
      }
      std::size_t total_bytes = bitio_stop_encoding(compressed_data_header(array.extent()[0],
                                                                          array.extent()[1],
                                                                          log2_threshold,
                                                                          encoded_symbols_,
                                                                          encode_extent),
                                                    data_buffer);
      
      // Update size/rms
      *actual_size = total_bytes;
      *actual_rms  = std::sqrt(double(current_encoding_square_error_)/double(array.count()));
    }

    /**
     *  Compress array and encode it in given buffer. Compression terminates when
     *  the size equals or exceeds target size. The input array is assumed to encode a
     *  hierarchical transformation using a pyramidal structure with root at (0,0)
     *  and children of (i,j) at (2i..2i+1)-(2j-2j+1). If encode_extent is true,
     *  the encoded data is totally self contained, otherwise the user should provide
     *  the matrix extent when decoding.
     */
    void compress_to_target_size(const int_matrix_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t *actual_size,
                                 double *actual_rms,
                                 bool encode_extent = true) {
      assert(buf);
      assert(buf_size);
      assert(buf_size >= target_size);
      assert(actual_size);
      assert(actual_rms);
      compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms, encode_extent); 
    }

    /**
     *  Compress array and encode it in given buffer. Compression terminates when
     *  the root mean square error falls below the target value. The input array is
     *  assumed to encode a hierarchical transformation using a pyramidal structure with
     *  root at (0,0) and children of (i,j) at (2i..2i+1)-(2j-2j+1). If encode_extent is true,
     *  the encoded data is totally self contained, otherwise the user should provide
     *  the matrix extent when decoding.
     */
    void compress_to_target_rms(const int_matrix_t& array,
                                double target_rms,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t *actual_size,
                                double *actual_rms,
                                bool encode_extent = true) {
      assert(buf);
      assert(buf_size);
      assert(actual_size);
      assert(actual_rms);
      compress(array, buf_size > 32 ? buf_size -16 : buf_size, target_rms, buf, buf_size, actual_size, actual_rms, encode_extent); 
    }

    /**
     *  Decompress buffer and store result into array. If extents are different than std::size_t(-1),
     *  it is assumed that the encoder did not store the extend in the buffer, and the provided
     *  extent is used to resize the matrix. Otherwise, the right extent is read from the buffer.
     */
    void decompress(int_matrix_t& array,
                    const void* buf,
                    std::size_t buf_size,
                    std::size_t extent_0 = std::size_t(-1),
                    std::size_t extent_1 = std::size_t(-1)) {
      assert(buf);
      assert(buf_size);

      uint8_t* code_buffer = (uint8_t*)buf;
      compressed_data_header hdr = bitio_start_decoding(code_buffer, buf_size, extent_0, extent_1);
      {
        array.resize(subscript_t(hdr.h(), hdr.w()));
        array.clear();

        subscript_vector_t subordinate_list; subordinate_list.reserve(array.count());
        int_t threshold = 1 << hdr.log2_threshold();
        while (threshold!=0) {
          ezw_decode_pass(array,threshold,subordinate_list);
          threshold >>= 1;
        }
      }
      bitio_stop_decoding();
    }

  protected: // EZW coding

    /**
     *  Compute map for zerotree detection based on
     *  "A Fast Technique for Identifying Zerotrees in the EZW Algorithm"
     *  by Jerome M. Shapiro, IEEE International Conference on Acoustics, Speech, and Signal Processing, 1996.
     *  ICASSP-96. Conference Proceedings, Volume 3, 7-10 May 1996:1455 - 1458
     */
    void zerotree_map_in(int_matrix_t& zerotree_map,
                         const int_matrix_t& m) const {
      const std::size_t h = m.extent()[0];
      const std::size_t w = m.extent()[1];
      
      // Initialize the zerotree map by replacing each coefficient with the
      // largest power of two smaller than the corresponding magnitude.
      for (std::size_t i=0; i<h; ++i) {
        for (std::size_t j=0; j<w; ++j) {
          const int_t abs_m_ij = iabs(m(i,j));
          int_t zmap_ij = 0;
          if (abs_m_ij>0) {
            zmap_ij = bitops<int_t>::prev_power2(abs_m_ij);
          }
          zerotree_map(i,j) = zmap_ij;
        }
      }
      
      // For each coarser scale, replace each coefficient with
      // the bitwise OR of itself and its children
      zerotree_map_pull(zerotree_map, 0, 0);
    }
    
    /**
     * Compute all values of zerotree map rooted at j, j.
     * For each coarser scale, replace each coefficient with
     * the bitwise OR of itself and its children
     */
    void zerotree_map_pull(int_matrix_t& zerotree_map,
                           std::size_t i,
                           std::size_t j) const {
      assert(zerotree_map.good_subscript(subscript_t(i,j)));

      // Filter children, and OR with this
      const std::size_t min_child_i = i << 1;
      const std::size_t max_child_i = min_child_i+1;
      const std::size_t min_child_j = j << 1;
      const std::size_t max_child_j = min_child_j+1;
      for (std::size_t child_j=min_child_j; child_j<=max_child_j; ++child_j) {
        for (std::size_t child_i=min_child_i; child_i<=max_child_i; ++child_i) {
          if ((!(child_i==i && child_j==j)) &&
                zerotree_map.good_subscript(subscript_t(child_i,child_j))) {
            zerotree_map_pull(zerotree_map, child_i, child_j);
            zerotree_map(i,j) |= zerotree_map(child_i, child_j);
          }
        }
      }
    }
    
    bool stop_encoding(std::size_t target_size,
                       energy_t target_square_error) const {
      return
        (bitio_encode_byte_count()>=target_size) ||
        (current_encoding_square_error_ <= target_square_error);
    }
 
    static bool is_zerotree(const int_matrix_t& m,
                            std::size_t i,
                            std::size_t j,
                            int_t threshold) {
      assert(m.good_subscript(subscript_t(i,j)));
      const int_t abs_m_ij = iabs(m(i,j));
      
      bool result =
        (abs_m_ij < threshold) ||  // is insignificant at this level or
        (abs_m_ij >= 2*threshold); // was significant at previous levels 
      
      const std::size_t min_child_i = i << 1;
      const std::size_t max_child_i = min_child_i+1;
      const std::size_t min_child_j = j << 1;
      const std::size_t max_child_j = min_child_j+1;
      for (std::size_t child_j=min_child_j; child_j<=max_child_j && result; ++child_j) {
        for (std::size_t child_i=min_child_i; child_i<=max_child_i && result; ++child_i) {
          if ((!(child_i==i && child_j==j)) && m.good_subscript(subscript_t(child_i,child_j))) {
            result = is_zerotree(m,child_i,child_j,threshold);
          }
        }
      }
      return result;
    }
    
    ezw_code_t ezw_code(const int_matrix_t& m,
                        const int_matrix_t& zerotree_map,
                        std::size_t i,
                        std::size_t j,
                        int_t threshold) const {
      ezw_code_t result = EZW_ISOLATED_ZERO;
      const int_t m_ij = m(i,j);
      const int_t abs_m_ij = iabs(m_ij);
      
      if (abs_m_ij >= threshold && abs_m_ij < 2*threshold) {
        result = (m_ij>=0) ? EZW_POSITIVE : EZW_NEGATIVE;
#if 0
      } else if (is_zerotree(m,i,j,threshold)) {
         result = EZW_ZEROTREE_ROOT;
#else
      } else if ((zerotree_map(i,j) & threshold) == 0) {
        result = EZW_ZEROTREE_ROOT;
#endif
      }
      return result;
    }
        
    void ezw_encode_pass(const int_matrix_t& m,
                         int_t threshold,
                         const int_matrix_t& zerotree_map,
                         std::size_t target_size,
                         energy_t target_square_error,
                         int_vector_t& subordinate_list) {
      const int_t estimate  = threshold + threshold/2;

      // Dominant pass
      scanning_queue_init(m.extent()[0], m.extent()[1]);
      bool previous_is_significant = false;
      while (!scanning_queue_empty() && !stop_encoding(target_size, target_square_error)) {
        const subscript_t idx                   = scanning_queue_front().first;
        bool              parent_is_significant = scanning_queue_front().second;
        scanning_queue_pop();
        
        const ezw_code_t code = ezw_code(m,zerotree_map,idx[0],idx[1],threshold);
        bitio_encode(code, parent_is_significant, previous_is_significant);

        previous_is_significant = (code == EZW_POSITIVE || code == EZW_NEGATIVE);
        
        if (code == EZW_POSITIVE || code == EZW_NEGATIVE) {
          const int_t abs_m_idx = iabs(m(idx));
          subordinate_list.push_back(abs_m_idx);
          current_encoding_square_error_ += -isqr(abs_m_idx) + isqr(energy_t(abs_m_idx) - energy_t(estimate));
        }
        
        // Process quadtree children
        if (code!=EZW_ZEROTREE_ROOT) {
          parent_is_significant = code != EZW_ISOLATED_ZERO;
          const std::size_t min_child_i = idx[0] << 1;
          const std::size_t max_child_i = min_child_i+1;
          const std::size_t min_child_j = idx[1] << 1;
          const std::size_t max_child_j = min_child_j+1;
          for (std::size_t child_j=min_child_j; child_j<=max_child_j; ++child_j) {
            for (std::size_t child_i=min_child_i; child_i<=max_child_i; ++child_i) {
              const subscript_t child_idx = subscript_t(child_i, child_j);
              if ((child_idx != idx) && m.good_subscript(child_idx)) {
                scanning_queue_push(std::make_pair(child_idx, parent_is_significant));
              }
            }
          }
        }
      }
      // FIXME std::cerr << "\n";

      // Subordinate pass
      threshold >>= 1;
      if (threshold>0) {
        const int_t delta_one  = -threshold + threshold + threshold/2;
        const int_t delta_zero = -threshold + threshold/2;
        const std::size_t subordinate_list_size =subordinate_list.size(); 
        for (std::size_t i=0; i<subordinate_list_size && !stop_encoding(target_size, target_square_error); ++i) {
          const bool bit = ((subordinate_list[i] & threshold)!=0);
          bitio_encode(bit);
          const int_t target_value = subordinate_list[i];
          const int_t old_approx = (target_value & (~((threshold<<1)-1))) + threshold; 
          const int_t new_approx = bit ? (old_approx+delta_one) : (old_approx+delta_zero);
          current_encoding_square_error_ += -isqr(energy_t(old_approx)-energy_t(target_value)) + isqr(energy_t(new_approx)-energy_t(target_value));
        }
        // FIXME std::cerr << "\n";
      }
    }
    
  protected: // EZW decoding

    bool stop_decoding() const {
      return
        (decoded_symbols_ == encoded_symbols_);        
    }
    
    void ezw_decode_pass(int_matrix_t& m,
                         int_t threshold,
                         subscript_vector_t& subordinate_list) {
      const int_t estimate  = threshold + threshold/2;

      scanning_queue_init(m.extent()[0], m.extent()[1]);
      bool previous_is_significant = false;
      while (!scanning_queue_empty() && !stop_decoding()) {
        const subscript_t idx                   = scanning_queue_front().first;
        bool              parent_is_significant = scanning_queue_front().second;
        scanning_queue_pop();
        
        const ezw_code_t code = bitio_decode(parent_is_significant, previous_is_significant);

        if (code == EZW_POSITIVE) { 
          subordinate_list.push_back(idx);
          m(idx) =  estimate;
        } else if (code == EZW_NEGATIVE) {
          subordinate_list.push_back(idx);
          m(idx) = -estimate;
        }

        previous_is_significant = (code == EZW_POSITIVE || code == EZW_NEGATIVE);

        // Process quadtree children
        if (code!=EZW_ZEROTREE_ROOT) {
          parent_is_significant = code != EZW_ISOLATED_ZERO;
          const std::size_t min_child_i = idx[0] << 1;
          const std::size_t max_child_i = min_child_i+1;
          const std::size_t min_child_j = idx[1] << 1;
          const std::size_t max_child_j = min_child_j+1;
          for (std::size_t child_j=min_child_j; child_j<=max_child_j; ++child_j) {
            for (std::size_t child_i=min_child_i; child_i<=max_child_i; ++child_i) {
              const subscript_t child_idx = subscript_t(child_i, child_j);
              if ((child_idx != idx) && m.good_subscript(child_idx)) {
                scanning_queue_push(std::make_pair(child_idx, parent_is_significant));
              }
            }
          }
        }
      }
      // FIXME std::cerr << "\n";

      // Subordinate pass
      threshold >>= 1;
      if (threshold>0) {
        const int_t delta_one  = -threshold + threshold + threshold/2;
        const int_t delta_zero = -threshold + threshold/2;
        const std::size_t subordinate_list_size = subordinate_list.size();
        for (std::size_t i=0; i<subordinate_list_size && !stop_decoding(); ++i) {
          const bool bit = bitio_decode();
          const subscript_t& idx = subordinate_list[i];
          if (bit) {
            if (m(idx)>=0) {
              m(idx) += delta_one;
            } else {
              m(idx) -= delta_one;
            }
          } else {
            if (m(idx)>=0) {
              m(idx) += delta_zero;
            } else {
              m(idx) -= delta_zero;
            }
          }
        }          
        // FIXME std::cerr << "\n";
      }
    }
  };
}

#endif
