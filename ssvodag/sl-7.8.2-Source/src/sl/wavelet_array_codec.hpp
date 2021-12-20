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
#ifndef SL_WAVELET_ARRAY_CODEC_HPP
#define SL_WAVELET_ARRAY_CODEC_HPP

#include <sl/array_codec.hpp>
#include <sl/wavelet_transform.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/bitops.hpp>
#include <sl/encdec.hpp>
#include <cassert>
#include <limits>
#include <algorithm>

namespace sl {

  /**
   * Base class for lossy and losselss array codecs
   * based on normalized wavelet transforms
   */
  class wavelet_array_codec: public array_codec {
  public:
    typedef wavelet_array_codec this_t;
    typedef array_codec                    super_t;
    typedef super_t::int8_array2_t         int8_array2_t;
    typedef super_t::uint8_array2_t        uint8_array2_t;
    typedef super_t::int16_array2_t        int16_array2_t;
    typedef super_t::uint16_array2_t       uint16_array2_t;
    typedef super_t::int32_array2_t        int32_array2_t;
    typedef super_t::uint32_array2_t       uint32_array2_t;
    typedef super_t::float_array2_t        float_array2_t;
    typedef float_array2_t::subscript_t    subscript_t;

    typedef wavelet_transform<float>       wavelet_transform_t;
    
    typedef enum {
      TRANSFORM_KIND_IDENTITY    = 0,
      TRANSFORM_KIND_HAAR        = 1,
      TRANSFORM_KIND_DAUBECHIES4 = 2,
      TRANSFORM_KIND_SYMMLET8    = 3,
      TRANSFORM_KIND_AUTO        = 4,
      TRANSFORM_KIND_COUNT = TRANSFORM_KIND_AUTO
    } transform_kind_t;
    
  protected: // Instance variables

    std::string             transform_name_[TRANSFORM_KIND_COUNT];
    wavelet_transform_t*    transform_table_[TRANSFORM_KIND_COUNT];
    transform_kind_t        transform_kind_;

  protected: // Helpers

    static inline uint16_t pow2size(uint16_t n) {
      uint16_t result = 1;
      while (result<n) result <<=1;
      return result;
    }

    static inline uint16_t wavelet_size(uint16_t n, const wavelet_transform_t* wxform) {
      uint16_t result = 1;
      while (result<n) result <<=1;
      if (result>n) {
        const uint16_t k_analysis  = uint16_t(std::max(-wxform->analysis_low_first_index(),
                                                          -wxform->analysis_high_first_index()));
        const uint16_t k_synthesis = uint16_t(std::max(-wxform->synthesis_low_first_index(),
                                                       -wxform->synthesis_high_first_index()));
        
        const uint16_t k_wavelet   = std::max(k_analysis, k_synthesis);

        result = std::min(uint16_t((result>>1)+(n>>1)+k_wavelet+1), result);
      }
      return result;
    }

    template <class T>
    static void symmetric_extend_top_left(dense_array<T,2,void>& array,
                                          std::size_t h,
                                          std::size_t w) {
      std::size_t h2 = array.extent()[0];
      std::size_t w2 = array.extent()[1];
      
      // Extend by mirroring non pow2 upper left
      for (std::size_t i=h; i<h2; ++i) {
        for (std::size_t j=0; j<w; ++j) {
          array(i,j) = array(h-2-(i-h),j);
        }
      }
      for (std::size_t j=w; j<w2; ++j) {
        for (std::size_t i=0; i<h2; ++i) {
          array(i,j) = array(i, w-2-(j-w));
        }
      }
    }

    template <class T>
    static void clear_extra_wavelet_coefficients(dense_array<T,2,void>& array,
                                                 std::size_t h,
                                                 std::size_t w,
                                                 const wavelet_transform_t* wxform) {
      std::size_t h2 = array.extent()[0];
      std::size_t w2 = array.extent()[1];
      std::size_t k_analysis  = std::size_t(std::max(-wxform->analysis_low_first_index(),
                                                     -wxform->analysis_high_first_index()));
      std::size_t k_synthesis = std::size_t(std::max(-wxform->synthesis_low_first_index(),
                                                     -wxform->synthesis_high_first_index()));
      
      std::size_t k_wavelet   = std::max(k_analysis, k_synthesis);

      for (std::size_t h_wavelet=pow2size(h); h_wavelet>2; h_wavelet>>=1) {
        std::size_t h_detail = h_wavelet>>1;
        for (std::size_t i=h_detail+h/2+1+k_wavelet; i<std::min(h_wavelet,h2); ++i) {
          for (std::size_t j=0; j<w2; ++j) {
            array(i,j) = scalar_math<T>::zero();
          }
        }
      }
      
      for (std::size_t w_wavelet=pow2size(w); w_wavelet>2; w_wavelet>>=1) {
        std::size_t w_detail = w_wavelet>>1;
        for (std::size_t j=w_detail+w/2+1+k_wavelet; j<std::min(w_wavelet,w2); ++j) {
          for (std::size_t i=0; i<h2; ++i) {
            array(i,j) = scalar_math<T>::zero();
          }
        }
      }
    }
      
  public: // Construction and destruction

    /// Construct codec
    wavelet_array_codec() {
      transform_table_[TRANSFORM_KIND_IDENTITY]    = new identity_wavelet_transform<float>(); 
      transform_table_[TRANSFORM_KIND_HAAR]        = new normalized_haar_wavelet_transform<float>(); 
      transform_table_[TRANSFORM_KIND_DAUBECHIES4] = new normalized_daubechies4_wavelet_transform<float>();
      transform_table_[TRANSFORM_KIND_SYMMLET8]    = new normalized_symmlet8_wavelet_transform<float>();
      transform_name_[TRANSFORM_KIND_IDENTITY]     = std::string("Identity");
      transform_name_[TRANSFORM_KIND_HAAR]         = std::string("Haar");
      transform_name_[TRANSFORM_KIND_DAUBECHIES4]  = std::string("Daub4");
      transform_name_[TRANSFORM_KIND_SYMMLET8]     = std::string("Symmlet8");
      transform_kind_ = TRANSFORM_KIND_DAUBECHIES4;
    }

    /// Destruct codec
    ~wavelet_array_codec() {
      for (std::size_t i=0; i<TRANSFORM_KIND_COUNT; ++i) delete transform_table_[i];
    }

  public:

    /// Set current transform kind for compression to k
    inline void set_transform_kind(transform_kind_t k) {
      assert(k<=TRANSFORM_KIND_COUNT);
      transform_kind_ = k;
    }

    /**
     *  The current transform kind used for compression. If auto, different kinds are
     * tested and the 'best' one is chosen.
     */
    transform_kind_t transform_kind() const {
      return transform_kind_;
    }

    inline transform_kind_t first_transform_kind() const {
      return TRANSFORM_KIND_IDENTITY;
    }

    inline transform_kind_t last_transform_kind() const {
      return TRANSFORM_KIND_SYMMLET8;
    }
    
  protected:

    /**
     *  Compress array and encode it in given buffer using current transform kind.
     *  Compression terminates when the size equals or exceeds target size or when the root mean
     *  square error falls below the target value. The array is transformed by the current
     *  transform kind and then compressed by the method defined in subclass. If the
     *  current transform kind is AUTO, various possibilities
     *  are tested and the best one  is chosen.
     */
    virtual void float_compress(const float_array2_t& array,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      assert(buf);
      assert(buf_size >= target_size);

      transform_kind_t tkind = transform_kind();
      if (tkind != TRANSFORM_KIND_AUTO) {
        // Compress using current transform kind
        float_wavelet_compress(tkind,
                               array,
                               target_size, target_rms_error,
                               buf, buf_size,
                               actual_size,
                               actual_rms_error);
      } else {
        // Try Identity first
        float_wavelet_compress(TRANSFORM_KIND_IDENTITY,
                               array,
                               target_size,
                               target_rms_error,
                               buf, buf_size,
                               actual_size,
                               actual_rms_error);
        // Try others, select them if improved
        std::vector<uint8_t> tmp_buffer(buf_size);
        std::size_t          tmp_size = 0;
        double               tmp_rms_error  = 0;
        for (std::size_t i = 1; i<TRANSFORM_KIND_COUNT; ++i) {
          transform_kind_t tkind = transform_kind_t(i);
          float_wavelet_compress(tkind,
                                 array,
                                 target_size, target_rms_error,
                                 &(tmp_buffer[0]), buf_size,
                                 &tmp_size,
                                 &tmp_rms_error);
          bool improved = false;

          if (*actual_size <= target_size+4 && *actual_rms_error <= 1.1*target_rms_error) {
            // *BEST* passes criteria
            if (tmp_size <= target_size+4 && tmp_rms_error <= 1.1*target_rms_error) {
              // *CURRENT* and *BEST* pass all criteria
              if (tmp_size < *actual_size) { // FIXME - use a better selection scheme
                improved = true;
              }
            }
          } else if (tmp_size <= target_size+4 && tmp_rms_error <= 1.1*target_rms_error) {
            // *CURRENT* passes all criteria, *BEST* misses them
            improved = true;
          } else {
            // *NONE* passes criteria
            if (tmp_size <= target_size+4 && *actual_size <= target_size+4) {
              if (tmp_rms_error < *actual_rms_error) {
                // RMS_ERROR improvement
                improved = true;
              }
            } else if (tmp_rms_error <= 1.1*target_rms_error &&  *actual_rms_error <= 1.1*target_rms_error) {
              if (tmp_size < *actual_size) {
                // Size improvement
                improved = true;
              }
            }
          }

          if (improved) {
            std::copy(tmp_buffer.begin(), tmp_buffer.begin()+tmp_size, (uint8_t*)buf);
            *actual_size = tmp_size;
            *actual_rms_error  = tmp_rms_error;
          }
        }
      }
    }

    virtual void float_wavelet_compress(transform_kind_t tkind,
                                        const float_array2_t& array,
                                        std::size_t target_size,
                                        double target_rms_error,
                                        void* buf,
                                        std::size_t buf_size,
                                        std::size_t* actual_size,
                                        double* actual_rms_error) = 0;

  };
  
}

#endif
