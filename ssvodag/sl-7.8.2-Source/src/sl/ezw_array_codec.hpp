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
#ifndef SL_EZW_ARRAY_CODEC_HPP
#define SL_EZW_ARRAY_CODEC_HPP

#include <sl/wavelet_array_codec.hpp>
#include <sl/embedded_zerotree_codec.hpp>

namespace sl {

  /**
   *  Lossy compression of arbitrary scalar arrays using a wavelet
   *  transformed followed by arithmetic encoding of zerotrees of
   *  wavelet coefficients.
   *
   *  Based on "Embedded Image Coding Using Zerotrees of Wavelet Coefficients"
   *  by Jerome M. Shapiro, IEEE Transactions on Signal Processing, Vol.41, No.12,
   *  December 1993: 3445-3462.
   *
   *  Compression to targer RMS based on "Limiting the Distortion of a Wavelet
   *  Image Codec", Joonas Lehtinen, Acta Cybernetica, Vol 14, 1999: 341-356.
   */
  class ezw_array_codec: public wavelet_array_codec {
  public:
    typedef ezw_array_codec                this_t;
    typedef wavelet_array_codec            super_t;
    typedef super_t::int8_array2_t         int8_array2_t;
    typedef super_t::uint8_array2_t        uint8_array2_t;
    typedef super_t::int16_array2_t        int16_array2_t;
    typedef super_t::uint16_array2_t       uint16_array2_t;
    typedef super_t::int32_array2_t        int32_array2_t;
    typedef super_t::uint32_array2_t       uint32_array2_t;
    typedef super_t::float_array2_t        float_array2_t;
    typedef super_t::subscript_t           subscript_t;
    typedef super_t::wavelet_transform_t   wavelet_transform_t;
    typedef super_t::transform_kind_t      transform_kind_t;

    typedef embedded_zerotree_codec<int32_t> int32_embedded_zerotree_codec_t;
    
  protected: // Compressed header

    class compressed_data_header {
    protected:
      std::size_t      h_;
      std::size_t      w_;
      transform_kind_t transform_kind_;
      float            signal_mean_;
      float            xform_amax_;
    public:
      compressed_data_header(std::size_t h,
                             std::size_t w,
                             const transform_kind_t& transform_kind,
                             float signal_mean,
                             float xform_amax)
          :
          h_(h),
          w_(w),
          transform_kind_(transform_kind),
          signal_mean_(signal_mean),
          xform_amax_(xform_amax) {
      }

      compressed_data_header(const void* header_buffer) {
        const uint8_t* buf = (const uint8_t*)header_buffer;
        uint8_t hdr = be_decode<uint8_t>(buf); buf += 1;
        const bool square_matrix  = (hdr & 0x80)!=0;
        const bool h_16bit        = (hdr & 0x40)!=0;
        const bool w_16bit        = (hdr & 0x20)!=0;
        const bool zero_mean      = (hdr & 0x10)!=0;
        transform_kind_ = transform_kind_t(hdr & (~(0x80|0x40|0x20|0x10)));
        if (h_16bit) { h_ = be_decode<uint16_t>(buf); buf += 2; } else { h_ = be_decode<uint8_t>(buf); buf += 1; }
        if (square_matrix) {
          w_ = h_;
        } else {
          if(w_16bit) { w_ = be_decode<uint16_t>(buf); buf += 2; } else { w_ = be_decode<uint8_t>(buf); buf += 1; }
        }
        if (zero_mean) {
          signal_mean_ = 0.0f;
        } else {
          signal_mean_ = be_decode<float>(buf); buf+=4; // FIXME do not encode/decode a zero
        }
        xform_amax_     = be_decode<float>(buf); buf+=4;
        assert(std::size_t(buf-(const uint8_t*)header_buffer) == compressed_size());
      }

      void encode(void* header_buffer) {
        uint8_t* buf = (uint8_t*)header_buffer;
        const bool square_matrix = h() == w();
        const bool h_16bit = h() > 0xff;
        const bool w_16bit = w() > 0xff;
        const bool zero_mean = (signal_mean_ == 0.0f);
        uint8_t hdr = uint8_t(transform_kind());
        if (square_matrix) hdr |= 0x80;
        if (h_16bit) hdr |= 0x40;
        if (w_16bit) hdr |= 0x20;
        if (zero_mean) hdr |= 0x10;

        buf += be_encode(hdr, buf);
        if (h_16bit) { buf += be_encode(uint16_t(h()), buf); } else { buf += be_encode(uint8_t(h()), buf); }
        if (!square_matrix) {
          if (w_16bit) { buf += be_encode(uint16_t(w()), buf); } else { buf += be_encode(uint8_t(w()), buf); }
        }
        if (!zero_mean) buf += be_encode(signal_mean(), buf); 
        buf += be_encode(xform_amax(), buf);
        assert(std::size_t(buf-(uint8_t*)header_buffer) == compressed_size());
      }

      inline std::size_t h() const { return h_; }
      inline std::size_t w() const { return w_; }
      inline transform_kind_t transform_kind() const { return transform_kind_; }
      inline float signal_mean() const { return signal_mean_; }
      inline float xform_amax() const { return xform_amax_; }

      inline std::size_t compressed_size() const {
        const bool square_matrix = h() == w();
        const bool h_16bit = h() > 0xff;
        const bool w_16bit = w() > 0xff;
        const bool zero_mean = (signal_mean_ == 0.0f);
        return
          1 +
          (h_16bit ? 2 : 1) +
          (square_matrix ? 0 : (w_16bit ? 2 : 1)) +
          (zero_mean ? 0 : 4) +
          4;
      }
    };
    
  protected: // Instance variables

    int32_embedded_zerotree_codec_t    embedded_zerotree_codec_;
    
  protected: // Helpers

    static inline float int32_quantize_threshold(float amax) {
      enum { Int32_quantize_bits = 23 };
      const float Int32_quantize_scale     =(float)(1<<Int32_quantize_bits);
      return amax ? amax/Int32_quantize_scale : 1.0f;
    }
      
  public: // Construction and destruction

    /// Construct codec
    ezw_array_codec() {
    }

    /// Destruct codec
    ~ezw_array_codec() {
    }

    virtual std::string description() const { return std::string("EZW") + "/" + this->transform_name_[transform_kind_]; }
    
  public:

    /**
     *  Compress array and encode it in given buffer using given transform kind.
     *  Compression terminates when the size equals or exceeds target size or when the root mean
     *  square error falls below the target value. The array is transformed by the current
     *  transform kind and then compressed by arithmetic encoding zerotrees of
     *  wavelet coefficients.
     */
    void float_wavelet_compress(transform_kind_t tkind,
                                const float_array2_t& array,
                                std::size_t target_size,
                                double target_rms,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t *actual_size,
                                double *actual_rms) {
      assert(tkind < TRANSFORM_KIND_COUNT);
      assert(buf);
      assert(buf_size >= target_size);

      const uint16_t h = array.extent()[0];
      const uint16_t w = array.extent()[1];
      const uint16_t h_wavelet = wavelet_size(h,transform_table_[tkind]);
      const uint16_t w_wavelet = wavelet_size(w,transform_table_[tkind]);
      
      // Transform to zero-centered range
      
      float signal_mean = float(mean(array));
      float signal_amax = float(amax(array));
      if (abs(signal_mean)<0.1*target_rms ||
          abs(signal_mean)<0.1*signal_amax) {
        signal_mean = 0.0; // FIXME force zero mean to avoid encoding it...
      }
      
      float_array2_t xform(h_wavelet,w_wavelet);
      for (std::size_t i=0; i<h; ++i) {
        for (std::size_t j=0; j<w; ++j) {
          xform(i,j) = array(i,j) - signal_mean;
        }
      }

      // Extend to power of 2 dimension
      symmetric_extend_top_left(xform, h, w);
      
      // Wavelet transform
      transform_table_[tkind]->forward2d(xform,h_wavelet,w_wavelet);

      // Clear coefficients referring purely to virtual non pow2 part
      // FIXME clear_extra_wavelet_coefficients(xform, h, w, transform_table_[tkind]);
      
      // Convert to integer and compute distortion due to quantization
      double quantization_rms = 0.0f;
      const float xform_amax = amax(xform);
      const float qthr = int32_quantize_threshold(xform_amax);
      
      int32_array2_t ixform(h_wavelet,w_wavelet);
      for (std::size_t i=0; i<h_wavelet; ++i) {
        for (std::size_t j=0; j<w_wavelet; ++j) {
          ixform(i,j) = int32_quantize(xform(i,j), qthr);
          quantization_rms += sqr((double)xform(i,j)-(double)int32_dequantize(ixform(i,j), qthr));
        }
      }
      quantization_rms = std::sqrt(quantization_rms/double(h_wavelet*w_wavelet));

      // Encode header
      compressed_data_header hdr(h, w, tkind, signal_mean, xform_amax);

      uint8_t* header_buffer = (uint8_t*)buf;
      hdr.encode(header_buffer);
      const std::size_t header_bytes   = hdr.compressed_size();
      
      // Encode wavelet-transformed data
      const std::size_t max_code_bytes = buf_size - header_bytes; 
      const std::size_t target_code_size = target_size > header_bytes ? target_size-header_bytes : 0;
      uint8_t* code_buffer = header_buffer + header_bytes;

      std::size_t actual_code_size;
      double      actual_code_rms;
      embedded_zerotree_codec_.compress(ixform,
                                        target_code_size,
                                        double(int32_quantize(std::max((float)(target_rms-quantization_rms), 0.0f), qthr)),
                                        code_buffer,
                                        max_code_bytes,
                                        &actual_code_size,
                                        &actual_code_rms,
                                        false); // NO EXTENT
      *actual_size = actual_code_size + header_bytes;
      *actual_rms  = quantization_rms + int32_dequantize(int32_t(actual_code_rms), qthr);
    }

    
    /**
     *  Decompress buffer and store result into array.
     */
    void float_decompress(float_array2_t& array,
                          const void* buf,
                          std::size_t buf_size) {
      assert(buf);
      assert(buf_size >= 2);

      // Decode header
      const uint8_t* header_buffer = (const uint8_t*)buf;
      compressed_data_header hdr(header_buffer);
      const std::size_t header_bytes   = hdr.compressed_size();
      
      // Decode wavelet-transformed data
      const std::size_t code_bytes     = buf_size - header_bytes;
      const uint8_t* code_buffer       = header_buffer + header_bytes;

      const std::size_t h_wavelet = wavelet_size(hdr.h(), transform_table_[hdr.transform_kind()]);
      const std::size_t w_wavelet = wavelet_size(hdr.w(), transform_table_[hdr.transform_kind()]);

      int32_array2_t ixform;
      embedded_zerotree_codec_.decompress(ixform, code_buffer, code_bytes, h_wavelet, w_wavelet);

      float_array2_t xform(h_wavelet,w_wavelet);
      float qthr = int32_quantize_threshold(hdr.xform_amax());
      for (std::size_t i=0; i<h_wavelet; ++i) {
        for (std::size_t j=0; j<w_wavelet; ++j) {
          xform(i,j) = int32_dequantize(ixform(i,j), qthr);
        }
      }
      
      // Invert wavelet transform
      transform_table_[hdr.transform_kind()]->backward2d(xform,h_wavelet,w_wavelet);
      
      // Back-transform from zero-centered range
      array.resize(subscript_t(hdr.h(),hdr.w()));
      for (std::size_t i=0; i<hdr.h(); ++i) {
        for (std::size_t j=0; j<hdr.w(); ++j) {
          array(i,j) = xform(i,j) + hdr.signal_mean();
        }
      }
    }

  };
}


#endif
