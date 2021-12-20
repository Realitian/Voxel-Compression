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
#ifndef SL_MICRO_JPGLS_ARRAY_CODEC_HPP
#define SL_MICRO_JPGLS_ARRAY_CODEC_HPP

#include <sl/array_codec.hpp>

namespace sl {

  /**
   *  Simple/fast compression of arbitrary scalar arrays using a
   *  straightforward delta codec. Useful for lossless or
   *  near-losses compression with reasonably fast decompression.
   *
   *  The encoding method is derived from the 'regular mode' of
   *  JPEG-LS, described in:
   *
   *  FCD 14495, Lossless and near-lossless coding of continuous
   *  tone still images (JPEG-LS)(ISO/IEC JTC1/SC29 WG1 (JPEG/JBIG))
   *
   *  JPEG-LS FCD 14495 http://www.jpeg.org/public/fcd14495p.pdf
   */
  class micro_jpgls_array_codec: public array_codec {
  public:
    typedef micro_jpgls_array_codec       this_t;
    typedef array_codec                    super_t;
    typedef super_t::int8_array2_t         int8_array2_t;
    typedef super_t::uint8_array2_t        uint8_array2_t;
    typedef super_t::int16_array2_t        int16_array2_t;
    typedef super_t::uint16_array2_t       uint16_array2_t;
    typedef super_t::int32_array2_t        int32_array2_t;
    typedef super_t::uint32_array2_t       uint32_array2_t;
    typedef super_t::float_array2_t        float_array2_t;
    typedef float_array2_t::subscript_t    subscript_t;

    typedef enum { Error_kind_rms, Error_kind_amax } error_kind_t;

  protected: // BIT I/O
    
    uint8_t    *bitio_buffer_;
    std::size_t bitio_buffer_capacity_;
    uint8_t    *bitio_byte_pointer_;
    uint8_t     bitio_bit_index_;

    void bitio_start_encoder(void* buf, std::size_t buf_size);

    void bitio_start_decoder(const void* buf, std::size_t buf_size);
    
    void bitio_output_bit(bool b);

    bool bitio_input_bit();

    void bitio_output_uint(uint32_t x);

    void bitio_output_int(int32_t x);

    uint32_t bitio_input_uint();
    
    int32_t bitio_input_int();
    
    void bitio_stop_encoder();

    void bitio_stop_decoder();

    std::size_t bitio_bit_count() const;
    
    std::size_t bitio_byte_count() const;
    
    void bitio_output_limited_golomb(uint32_t k,
                                            uint32_t glimit,
                                            uint32_t qbpp,
                                            uint32_t value);

    uint32_t bitio_input_limited_golomb(uint32_t k,
                                               uint32_t glimit,
                                               uint32_t qbpp);

  protected: // Array access helpers
        
    template <class value_t>
    static value_t top_item(const dense_array<value_t,2,void>& array,
                                   uint32_t i,
                                   uint32_t j);

    template <class value_t>
    static value_t left_item(const dense_array<value_t,2,void>& array,
                                 uint32_t i,
                                 uint32_t j);

    template <class value_t>
    static value_t right_item(const dense_array<value_t,2,void>& array,
                                     uint32_t i,
                                     uint32_t j);
    template <class value_t>
    static value_t top_left_item(const dense_array<value_t,2,void>& array,
                                      uint32_t i,
                                      uint32_t j);

    template <class value_t>
    static value_t top_right_item(const dense_array<value_t,2,void>& array,
                                         uint32_t i,
                                         uint32_t j);

  protected: // JPEG-LS entry points

    enum {
      jls_context_top     = 1, // [-1,0,1]
      jls_context_regions = 2*jls_context_top+1,
      jls_context_count   = (jls_context_regions*jls_context_regions*jls_context_regions+1)/2,
    }; 

    int32_t jls_minval_; // Minimum value in input array
    int32_t jls_maxval_; // Maximum possible array value + minval
    int32_t jls_near_;   // Difference bound for near lossless encoding
    int32_t jls_near_times2_plus1_;
    
    int32_t jls_range_;  // Range of prediction error representation
    int32_t jls_qbpp_;   // Number of bits needed to represent a mapped error value
    int32_t jls_bpp_;    // Number of bits needed to represent jls_maxval_, with a minimum of 2
    int32_t jls_limit_;  // The value of glimit for a sample encoded in regular mode

    int32_t jls_n_[jls_context_count];      // Frequency of occurrence of each context
    int32_t jls_a_[jls_context_count];      // Accumulated prediction error magnitude 
    int32_t jls_b_[jls_context_count];      // Bias counter
    int32_t jls_c_[jls_context_count];      // Prediction correction value

    static uint32_t jls_j_[32];    // Run length codes
    static uint32_t jls_jn_[32];   // 1<<jls_j
    
    int32_t jls_reset_;     // Histogram update cycle length
    int32_t jls_minc_;      // Max negative value for bias update
    int32_t jls_maxc_;      // Max positive value for bias update
    
    void jls_init(int32_t in_min, int32_t in_max, int32_t eps);

    template <class int_t>
    void jls_init(const dense_array<int_t,2,void>& in, int32_t eps);
    
    static int32_t jls_quantize(int32_t x, int32_t eps, int32_t eps_times2_plus1);        
    static int32_t jls_dequantize(int32_t x, int32_t eps_times2_plus1);
    
    int32_t jls_quantize(int32_t x) const;

    int32_t jls_dequantize(int32_t x) const;

    static int32_t jls_predictor(int32_t a,
                                        int32_t b,
                                        int32_t c);
    
    static int32_t jls_context(int32_t a,
                                      int32_t b,
                                      int32_t c,
                                      int32_t d,
                                      int32_t eps);
    
    void jls_update(uint32_t q, int32_t errval);

    template <class int_t>
    void jls_encode(const dense_array<int_t,2,void>& in,
                    int32_t eps,
                    void *buf,
                    std::size_t buf_size,
                    std::size_t *actual_size,
                    double *actual_error,
                    error_kind_t error_kind);    
      
    template <class int_t>
    double jls_encoding_error(const dense_array<int_t,2,void>& in,
                              int32_t eps,
                              error_kind_t error_kind);
    
    template <class int_t>
    void jls_decode(dense_array<int_t,2,void>& array_prime,
                    const void *buf,
                    std::size_t buf_size);
    
  public: // Construction and destruction

    /// Construct codec
    micro_jpgls_array_codec() {
    }

    /// Destruct codec
    ~micro_jpgls_array_codec() {
    }

    virtual std::string description() const {
      return std::string("u-jpgls");
    }

  protected:
                              
    template <class int_t>
    uint32_t threshold_from_target_error(const dense_array<int_t,2,void>& array,
                                         double target_error,
                                         error_kind_t error_kind);

    template <class int_t>
    void int_compress(const dense_array<int_t,2,void>& array,
                      std::size_t target_size,
                      double target_error,
                      void* buf,
                      std::size_t buf_size,
                      std::size_t* actual_size, 
                      double *actual_error,
                      error_kind_t error_kind);  
    
    /**
     *  Decompress buffer and store result into array.
     */
    template <class int_t>
    void int_decompress(dense_array<int_t,2,void>& array,
                        const void* buf,
                        std::size_t buf_size);

  public:

    /**
     * Compress array to buffer
     */
    void int8_compress(const int8_array2_t& array,
                       std::size_t target_size,
                       double target_error,
                       void* buf,
                       std::size_t buf_size,
                       std::size_t* actual_size, 
                       double *actual_error);
    
    /**
     * Compress array to buffer
     */
    void int8_compress_amax(const int8_array2_t& array,
                            std::size_t target_size,
                            double target_error,
                            void* buf,
                            std::size_t buf_size,
                            std::size_t* actual_size, 
                            double *actual_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    void int8_decompress(int8_array2_t& array,
                         const void* buf,
                         std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    void int16_compress(const int16_array2_t& array,
                        std::size_t target_size,
                        double target_error,
                        void* buf,
                        std::size_t buf_size,
                        std::size_t* actual_size, 
                        double *actual_error);

    /**
     * Compress array to buffer
     */
    void int16_compress_amax(const int16_array2_t& array,
                             std::size_t target_size,
                             double target_error,
                             void* buf,
                             std::size_t buf_size,
                             std::size_t* actual_size, 
                             double *actual_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    void int16_decompress(int16_array2_t& array,
                          const void* buf,
                          std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    void int32_compress(const int32_array2_t& array,
                        std::size_t target_size,
                        double target_error,
                        void* buf,
                        std::size_t buf_size,
                        std::size_t* actual_size, 
                        double *actual_error);

    /**
     * Compress array to buffer
     */
    void int32_compress_amax(const int32_array2_t& array,
                             std::size_t target_size,
                             double target_error,
                             void* buf,
                             std::size_t buf_size,
                             std::size_t* actual_size, 
                             double *actual_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    void int32_decompress(int32_array2_t& array,
                          const void* buf,
                          std::size_t buf_size);

  public:

    /**
     * Compress array to buffer
     */
    void uint8_compress(const uint8_array2_t& array,
                       std::size_t target_size,
                       double target_error,
                       void* buf,
                       std::size_t buf_size,
                       std::size_t* actual_size, 
                       double *actual_error);

    /**
     * Compress array to buffer
     */
    void uint8_compress_amax(const uint8_array2_t& array,
                             std::size_t target_size,
                             double target_error,
                             void* buf,
                             std::size_t buf_size,
                             std::size_t* actual_size, 
                             double *actual_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    void uint8_decompress(uint8_array2_t& array,
                         const void* buf,
                         std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    void uint16_compress(const uint16_array2_t& array,
                        std::size_t target_size,
                        double target_error,
                        void* buf,
                        std::size_t buf_size,
                        std::size_t* actual_size, 
                        double *actual_error);

    /**
     * Compress array to buffer
     */
    void uint16_compress_amax(const uint16_array2_t& array,
                              std::size_t target_size,
                              double target_error,
                              void* buf,
                              std::size_t buf_size,
                              std::size_t* actual_size, 
                              double *actual_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    void uint16_decompress(uint16_array2_t& array,
                          const void* buf,
                          std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    void uint32_compress(const uint32_array2_t& array,
                         std::size_t target_size,
                         double target_error,
                         void* buf,
                         std::size_t buf_size,
                         std::size_t* actual_size, 
                         double *actual_error);

    /**
     * Compress array to buffer
     */
    void uint32_compress_amax(const uint32_array2_t& array,
                              std::size_t target_size,
                              double target_error,
                              void* buf,
                              std::size_t buf_size,
                              std::size_t* actual_size, 
                              double *actual_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    void uint32_decompress(uint32_array2_t& array,
                          const void* buf,
                          std::size_t buf_size);

  public: // Float

    /**
     *  Decompress buffer and store result into array.
     */
    void float_compress(const float_array2_t& array,
                        std::size_t  target_size,
                        double       target_error,
                        void* buf,
                        std::size_t buf_size,
                        std::size_t* actual_size,
                        double      *actual_error,
                        error_kind_t error_kind);
      
    /**
     *  Decompress buffer and store result into array.
     */
    void float_compress(const float_array2_t& array,
                        std::size_t  target_size,
                        double       target_error,
                        void* buf,
                        std::size_t buf_size,
                        std::size_t* actual_size,
                        double      *actual_error);

    /**
     *  Decompress buffer and store result into array.
     */
    void float_compress_amax(const float_array2_t& array,
                             std::size_t  target_size,
                             double       target_error,
                             void* buf,
                             std::size_t buf_size,
                             std::size_t* actual_size,
                             double      *actual_error);
    
    /**
     *  Decompress buffer and store result into array.
     */
    void float_decompress(float_array2_t& array,
                          const void* buf,
                          std::size_t buf_size);

  };
  
} // namespace sl

#endif

