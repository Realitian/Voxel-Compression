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
#ifndef SL_QUANTIZED_ARRAY_CODEC_HPP
#define SL_QUANTIZED_ARRAY_CODEC_HPP

#include <sl/array_codec.hpp>
#include <sl/fixed_rc_int_codec.hpp>
#include <sl/compact_bitio.hpp>

namespace sl {

  class quantized_array_codec: public array_codec {
  public:
    typedef quantized_array_codec          this_t;
    typedef array_codec                    super_t;
    typedef super_t::int8_array2_t         int8_array2_t;
    typedef super_t::uint8_array2_t        uint8_array2_t;
    typedef super_t::int16_array2_t        int16_array2_t;
    typedef super_t::uint16_array2_t       uint16_array2_t;
    typedef super_t::int32_array2_t        int32_array2_t;
    typedef super_t::uint32_array2_t       uint32_array2_t;
    typedef super_t::float_array2_t        float_array2_t;
    typedef float_array2_t::subscript_t    subscript_t;
    
    typedef fixed_rc_int_codec             bitio_t;

  protected:

    mutable bitio_t bitio_;

    // This flag is here for backward compatibility reasons!
    bool            is_compressing_header_;
    bool            is_delta_encoding_;

  public: // Encoding/Decoding
    
    inline quantized_array_codec(bool b = false) : is_compressing_header_(b), is_delta_encoding_(false) {}
    inline virtual ~quantized_array_codec() {}

    inline void set_is_compressing_header(bool b) { is_compressing_header_ = b; }
    inline bool is_compressing_header() const { return is_compressing_header_; }

    inline void set_is_delta_encoding(bool b) { is_delta_encoding_ = b; }
    inline bool is_delta_encoding() const { return is_delta_encoding_; }

    virtual std::string description() const { return "QUANT"; }

    virtual void float_compress(const float_array2_t& array,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error);

    virtual void float_compress_amax(const float_array2_t& iarray,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t* actual_size,
                                     double* actual_amax_error);

    virtual void float_decompress(float_array2_t& array,
                                  const void* buf,
                                  std::size_t buf_size);

    
    virtual void int8_compress(const int8_array2_t& iarray,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t* actual_size,
                               double* actual_rms_error);

    virtual void int8_compress_amax(const int8_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error);

    virtual void int8_decompress(int8_array2_t& iarray,
                                 const void* buf,
                                 std::size_t buf_size);

    virtual void uint8_compress(const uint8_array2_t& iarray,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error);

    virtual void uint8_compress_amax(const uint8_array2_t& iarray,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t* actual_size,
                                     double* actual_amax_error);

    virtual void uint8_decompress(uint8_array2_t& iarray,
                                  const void* buf,
                                  std::size_t buf_size);
    
    virtual void int16_compress(const int16_array2_t& iarray,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t* actual_size,
                               double* actual_rms_error);

    virtual void int16_compress_amax(const int16_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error);

    virtual void int16_decompress(int16_array2_t& iarray,
                                  const void* buf,
                                  std::size_t buf_size);

    virtual void uint16_compress(const uint16_array2_t& iarray,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error);
    
    virtual void uint16_compress_amax(const uint16_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error);

    virtual void uint16_decompress(uint16_array2_t& iarray,
                                   const void* buf,
                                   std::size_t buf_size);
        
    virtual void int32_compress(const int32_array2_t& iarray,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t* actual_size,
                               double* actual_rms_error);

    virtual void int32_compress_amax(const int32_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error);

    virtual void int32_decompress(int32_array2_t& iarray,
                                  const void* buf,
                                  std::size_t buf_size);

    virtual void uint32_compress(const uint32_array2_t& iarray,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error);
    
    virtual void uint32_compress_amax(const uint32_array2_t& iarray,
                                      std::size_t target_size,
                                      double target_amax_error,
                                      void* buf,
                                      std::size_t buf_size,
                                      std::size_t* actual_size,
                                      double* actual_amax_error);

    virtual void uint32_decompress(uint32_array2_t& iarray,
                                   const void* buf,
                                   std::size_t buf_size);
  };
  
}

#endif
