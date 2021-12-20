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
#ifndef SL_HIDWT_ARRAY_CODEC_HPP
#define SL_HIDWT_ARRAY_CODEC_HPP

#include <sl/dense_array.hpp>
#include <sl/indexed_functions.hpp>
#include <sl/generative_types.hpp>
#include <sl/assert.hpp>
#include <sl/fixed_ac_int_codec.hpp>
#include <sl/fixed_rc_int_codec.hpp>
#include <sl/fixed_huffman_rle_codec.hpp>
#include <sl/array_codec.hpp>
#include <sl/bitops.hpp>
#include <sl/encdec.hpp>
#include <limits>
#include <algorithm>

namespace sl {
  
  /**
   *  Simple/fast compression of arbitrary scalar arrays using a
   *  hierarchical discrete wavelet transform that maps integers
   *  to integers.
   */
  class hidwt_array_codec: public array_codec {
  public:
    typedef hidwt_array_codec this_t;
    typedef array_codec                      super_t;
    typedef super_t::int8_array2_t           int8_array2_t;
    typedef super_t::uint8_array2_t          uint8_array2_t;
    typedef super_t::int16_array2_t          int16_array2_t;
    typedef super_t::uint16_array2_t         uint16_array2_t;
    typedef super_t::int32_array2_t          int32_array2_t;
    typedef super_t::uint32_array2_t         uint32_array2_t;
    typedef super_t::float_array2_t          float_array2_t;
    typedef float_array2_t::subscript_t      subscript_t;

    typedef fixed_rc_int_codec               entropy_codec_t;
    typedef entropy_codec_t::bit_context_t   bit_context_t;
    typedef entropy_codec_t::int_context_t   int_context_t;
    
  public: // Construction and destruction

    /// Construct codec
    inline hidwt_array_codec() {
    }

    /// Destruct codec
    inline ~hidwt_array_codec() {
    }

    virtual std::string description() const { return std::string("HIDWT/C22"); }

  protected:

    entropy_codec_t      entropy_codec_;
    std::vector<uint8_t> entropy_codec_buffer_;
        
  public: // Integer types - helpers

    typedef enum { Error_kind_rms, Error_kind_amax } error_kind_t;

  protected: // Integer types - helpers
    
    /**
     * Compress array to current codec buffer
     */
    template <class int_t>
    void int_quantized_compress(const dense_array<int_t,2,void>& array,
                                uint32_t threshold,
                                std::size_t *actual_size, 
                                double *actual_error,
                                error_kind_t error_kind = Error_kind_rms);

    template <class int_t>
    void int_compress(const dense_array<int_t,2,void>& array,
                      std::size_t target_size,
                      double target_rms_error,
                      void* buf,
                      std::size_t buf_size,
                      std::size_t *actual_size, 
                      double *actual_error,
                      error_kind_t error_kind = Error_kind_rms);
    
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
    virtual void int8_compress(const int8_array2_t& array,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t *actual_size, 
                               double *actual_rms);

    /**
     * Compress array to buffer
     */
    virtual void int8_compress_amax(const int8_array2_t& array,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t *actual_size, 
                                    double *actual_amax_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    virtual void int8_decompress(int8_array2_t& array,
                                 const void* buf,
                                 std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    virtual void int16_compress(const int16_array2_t& array,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t *actual_size, 
                                double *actual_rms);
    
    /**
     * Compress array to buffer
     */
    virtual void int16_compress_amax(const int16_array2_t& array,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t *actual_size, 
                                     double *actual_amax_error);

    /**
     *  Decompress buffer and store result into array.
     */
    virtual void int16_decompress(int16_array2_t& array,
                                  const void* buf,
                                  std::size_t buf_size);


    /**
     * Compress array to buffer
     */
    virtual void int32_compress(const int32_array2_t& array,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t *actual_size, 
                                double *actual_rms);

    /**
     * Compress array to buffer
     */
    virtual void int32_compress_amax(const int32_array2_t& array,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t *actual_size, 
                                     double *actual_amax_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    virtual void int32_decompress(int32_array2_t& array,
                                  const void* buf,
                                  std::size_t buf_size);

  public:

    /**
     * Compress array to buffer
     */
    virtual void uint8_compress(const uint8_array2_t& array,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t *actual_size, 
                                double *actual_rms);
    /**
     * Compress array to buffer
     */
    virtual void uint8_compress_amax(const uint8_array2_t& array,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t *actual_size, 
                                     double *actual_amax_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    virtual void uint8_decompress(uint8_array2_t& array,
                                  const void* buf,
                                  std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    virtual void uint16_compress(const uint16_array2_t& array,
                                 std::size_t target_size,
                                 double target_rms_error,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t *actual_size, 
                                 double *actual_rms);
    /**
     * Compress array to buffer
     */
    virtual void uint16_compress_amax(const uint16_array2_t& array,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t *actual_size, 
                                     double *actual_amax_error);
        
    /**
     *  Decompress buffer and store result into array.
     */
    virtual void uint16_decompress(uint16_array2_t& array,
                                   const void* buf,
                                   std::size_t buf_size);

    /**
     * Compress array to buffer
     */
    virtual void uint32_compress(const uint32_array2_t& array,
                                 std::size_t target_size,
                                 double target_rms_error,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t *actual_size, 
                                 double *actual_rms);
    /**
     * Compress array to buffer
     */
    virtual void uint32_compress_amax(const uint32_array2_t& array,
                                      std::size_t target_size,
                                      double target_amax_error,
                                      void* buf,
                                      std::size_t buf_size,
                                      std::size_t *actual_size, 
                                      double *actual_amax_error);

    /**
     *  Decompress buffer and store result into array.
     */
    virtual void uint32_decompress(uint32_array2_t& array,
                                   const void* buf,
                                   std::size_t buf_size);

  public:
    
   /**
     *  Decompress buffer and store result into array.
     */
   void float_compress(const float_array2_t& wm,
		       std::size_t  target_size,
		       double       target_rms,
		       void* buf,
		       std::size_t buf_size,
		       std::size_t *actual_size,
		       double      *actual_rms,
		       error_kind_t error_kind);

    /**
     *  Decompress buffer and store result into array.
     */
    virtual void float_compress(const float_array2_t& wm,
                                std::size_t  target_size,
                                double       target_rms,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t *actual_size,
                                double      *actual_rms);

    /**
     *  Decompress buffer and store result into array.
     */
    virtual void float_compress_amax(const float_array2_t& wm,
                                     std::size_t  target_size,
                                     double       target_rms,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t *actual_size,
                                     double      *actual_rms);
    
    /**
     *  Decompress buffer and store result into array.
     */
    virtual void float_decompress(float_array2_t& array,
                                  const void* buf,
                                  std::size_t buf_size);
    
  };

  
}

#endif
