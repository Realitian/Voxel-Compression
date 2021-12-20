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
#ifndef SL_ARRAY_CODEC_HPP
#define SL_ARRAY_CODEC_HPP

#include <sl/assert.hpp>
#include <sl/float_cast.hpp>
#include <sl/dense_array.hpp>
#include <sl/indexed_functions.hpp>
#include <sl/generative_types.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/bitops.hpp>
#include <sl/encdec.hpp>
#include <limits>
#include <algorithm>

namespace sl {

  /**
   *  Base class for all lossy and lossless array codecs
   */
  class array_codec {
  public:
    typedef array_codec this_t;
    
    typedef dense_array<int8_t,2,void>    int8_array2_t;
    typedef dense_array<uint8_t,2,void>   uint8_array2_t;
    typedef dense_array<int16_t,2,void>   int16_array2_t;
    typedef dense_array<uint16_t,2,void>  uint16_array2_t;
    typedef dense_array<int32_t,2,void>   int32_array2_t;
    typedef dense_array<uint32_t,2,void>  uint32_array2_t;
    typedef dense_array<float,2,void>     float_array2_t;
    typedef float_array2_t::subscript_t   subscript_t;

  public:

    inline array_codec() {}
    inline virtual ~array_codec() {}

    virtual std::string description() const { return "base"; }
    
  protected: // Integer to float conversion and vice-versa

    static inline int32_t int32_quantize(float x, float threshold) {
      if (x>threshold) {
        return int32_t(x/threshold);
      } else if (x<-threshold) {
        return -int32_t(-x/threshold);
      } else {
        // zero
        return int32_t(0);
      }
    }
    
    static inline float int32_dequantize(int32_t ix, float threshold) {
      return float(ix * threshold);
    }

    static inline float int32_quantize_square_error(float x, float threshold) {
      float x_prime = int32_dequantize(int32_quantize(x, threshold), threshold);
      return (x-x_prime)*(x-x_prime);
    }

    template <class INT_T>
    static void int_array_from_float_array_in(dense_array<INT_T,2,void>& iarray,
                                              const dense_array<float,2,void>& farray) {
      typedef INT_T int_t;
      const int_t imin = std::numeric_limits<int_t>::min();
      const int_t imax = std::numeric_limits<int_t>::max();
      
      std::size_t h = farray.extent()[0];
      std::size_t w = farray.extent()[1];
      iarray.resize(subscript_t(h,w));
      for (std::size_t i=0; i<h; ++i) {
        for (std::size_t j=0; j<w; ++j) {
          float fa_ij = farray(i,j);
          if (fa_ij<0) {
            iarray(i,j) = int_t(median(fa_ij-0.5f,float(imin),float(imax)));
          } else {
            iarray(i,j) = int_t(median(fa_ij+0.5f,float(imin),float(imax)));
          }
        }
      }
    }
    
    template <class INT_T>
    static void float_array_from_int_array_in(dense_array<float,2,void>& farray,
                                              const dense_array<INT_T,2,void>& iarray) {
      std::size_t h = iarray.extent()[0];
      std::size_t w = iarray.extent()[1];
      farray.resize(subscript_t(h,w));
      for (std::size_t i=0; i<h; ++i) {
        for (std::size_t j=0; j<w; ++j) {
          farray(i,j) = float(iarray(i,j));
        }
      }
    }

  protected:

    virtual void float_compress(const float_array2_t& array,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) = 0;

    virtual void float_decompress(float_array2_t& array,
                                  const void* buf,
                                  std::size_t buf_size) = 0;


    template <class T>
    void base_compress_to_target_amax_error(const dense_array<T, 2, void>& array,
                                            std::size_t /*target_size*/,
                                            double target_amax_error,
                                            void* buf,
                                            std::size_t buf_size,
                                            std::size_t* actual_size,
                                            double* actual_amax_error) {
      typedef dense_array<T, 2, void> array2_t;
      array2_t array_prime;

      double tol_0          = 0.0f;
      double rms_tol_0      = 0.0f;
      std::size_t sz_tol_0 = 0;
      double err_tol_0     = 0;
      double tol_1          = *actual_amax_error;
      double rms_tol_1      = 0.0f;
      std::size_t sz_tol_1 = 0;
      compress_to_target_rms_error(array,rms_tol_1,buf,buf_size,&sz_tol_1,&rms_tol_1);
      decompress(array_prime, buf, buf_size);
      double err_tol_1 = amax_diff(array,array_prime);
      std::size_t steps = 0;
      while (err_tol_1 < target_amax_error && steps < 6) { // FIXME
        ++steps;
        tol_0 = tol_1; rms_tol_0 = rms_tol_1; err_tol_0 = err_tol_1; sz_tol_0 = sz_tol_1;
        rms_tol_1 *= 2.0f;
        compress_to_target_rms_error(array,tol_1,buf,buf_size,&sz_tol_1,&rms_tol_1);
        decompress(array_prime, buf, buf_size);
        err_tol_1 = amax_diff(array,array_prime);
      }
      if (err_tol_1 <= target_amax_error) {
        *actual_amax_error = err_tol_1;
        *actual_size = sz_tol_1;
      } else {
        // err_tol_0 < target_amax_error < err_tol_1
        double  tol_mid = tol_1;
        double  rms_tol_mid  = rms_tol_1;
        double err_tol_mid = err_tol_1;
        std::size_t sz_tol_mid = 0;
        std::size_t steps = 0;
        while ((err_tol_1-err_tol_0)> 0.1*target_amax_error &&
               steps < 6) { // FIXME
          ++steps;
          tol_mid = 0.5 * (tol_1 + tol_0);
          compress_to_target_rms_error(array,tol_mid,buf,buf_size,&sz_tol_mid,&rms_tol_mid);
          decompress(array_prime, buf, buf_size);
          err_tol_mid = amax_diff(array,array_prime);
          if (err_tol_mid < target_amax_error) {
            tol_0 = tol_mid; rms_tol_0 = rms_tol_mid; err_tol_0 = err_tol_mid; sz_tol_0 = sz_tol_mid;
          } else {
	    tol_1 = tol_mid; rms_tol_1 = rms_tol_mid; err_tol_1 = err_tol_mid; sz_tol_1 = sz_tol_mid;
          }
        }
        if (err_tol_mid < target_amax_error) {
          compress_to_target_rms_error(array,tol_1,buf,buf_size,&sz_tol_1,&rms_tol_1);
        }
        *actual_amax_error = err_tol_1;
        *actual_size = sz_tol_1;
      }
      SL_USEVAR(rms_tol_0);
      SL_USEVAR(sz_tol_0);
    }

  protected:

    virtual void float_compress_amax(const float_array2_t& iarray,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t* actual_size,
                                     double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }
    
    virtual void int8_compress(const int8_array2_t& iarray,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t* actual_size,
                               double* actual_rms_error) {
      float_array2_t farray;
      float_array_from_int_array_in(farray, iarray);
      float_compress(farray, target_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    virtual void int8_compress_amax(const int8_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    virtual void int8_decompress(int8_array2_t& iarray,
                                 const void* buf,
                                 std::size_t buf_size) {
      float_array2_t farray;
      float_decompress(farray, buf, buf_size);
      int_array_from_float_array_in(iarray, farray);
    }

    virtual void uint8_compress(const uint8_array2_t& iarray,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      float_array2_t farray;
      float_array_from_int_array_in(farray, iarray);
      float_compress(farray, target_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    virtual void uint8_compress_amax(const uint8_array2_t& iarray,
                                     std::size_t target_size,
                                     double target_amax_error,
                                     void* buf,
                                     std::size_t buf_size,
                                     std::size_t* actual_size,
                                     double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    virtual void uint8_decompress(uint8_array2_t& iarray,
                                  const void* buf,
                                  std::size_t buf_size) {
      float_array2_t farray;
      float_decompress(farray, buf, buf_size);
      int_array_from_float_array_in(iarray, farray);
    }

  protected:
    
    virtual void int16_compress(const int16_array2_t& iarray,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t* actual_size,
                               double* actual_rms_error) {
      float_array2_t farray;
      float_array_from_int_array_in(farray, iarray);
      float_compress(farray, target_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    virtual void int16_compress_amax(const int16_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    virtual void int16_decompress(int16_array2_t& iarray,
                                  const void* buf,
                                  std::size_t buf_size) {
      float_array2_t farray;
      float_decompress(farray, buf, buf_size);
      int_array_from_float_array_in(iarray, farray);
    }

    virtual void uint16_compress(const uint16_array2_t& iarray,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      float_array2_t farray;
      float_array_from_int_array_in(farray, iarray);
      float_compress(farray, target_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }
    
    virtual void uint16_compress_amax(const uint16_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    virtual void uint16_decompress(uint16_array2_t& iarray,
                                   const void* buf,
                                   std::size_t buf_size) {
      float_array2_t farray;
      float_decompress(farray, buf, buf_size);
      int_array_from_float_array_in(iarray, farray);
    }
    
  protected:
    
    virtual void int32_compress(const int32_array2_t& iarray,
                               std::size_t target_size,
                               double target_rms_error,
                               void* buf,
                               std::size_t buf_size,
                               std::size_t* actual_size,
                               double* actual_rms_error) {
      float_array2_t farray;
      float_array_from_int_array_in(farray, iarray);
      float_compress(farray, target_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    virtual void int32_compress_amax(const int32_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    virtual void int32_decompress(int32_array2_t& iarray,
                                  const void* buf,
                                  std::size_t buf_size) {
      float_array2_t farray;
      float_decompress(farray, buf, buf_size);
      int_array_from_float_array_in(iarray, farray);
    }

    virtual void uint32_compress(const uint32_array2_t& iarray,
                                std::size_t target_size,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      float_array2_t farray;
      float_array_from_int_array_in(farray, iarray);
      float_compress(farray, target_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }
    
    virtual void uint32_compress_amax(const uint32_array2_t& iarray,
                                    std::size_t target_size,
                                    double target_amax_error,
                                    void* buf,
                                    std::size_t buf_size,
                                    std::size_t* actual_size,
                                    double* actual_amax_error) {
      base_compress_to_target_amax_error(iarray, target_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    virtual void uint32_decompress(uint32_array2_t& iarray,
                                   const void* buf,
                                   std::size_t buf_size) {
      float_array2_t farray;
      float_decompress(farray, buf, buf_size);
      int_array_from_float_array_in(iarray, farray);
    }

  public:

    void compress_to_target_size(const int8_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      int8_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const int8_array2_t& array,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      int8_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }
    
    void compress_to_target_amax_error(const int8_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      int8_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    void decompress(int8_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      int8_decompress(array, buf, buf_size);
    }
    
  public:
    
    void compress_to_target_size(const uint8_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      uint8_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const uint8_array2_t& array,
                                      double target_rms_error,
                                      void* buf,
                                      std::size_t buf_size,
                                      std::size_t* actual_size,
                                      double* actual_rms_error) {
      uint8_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_amax_error(const uint8_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      uint8_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    void decompress(uint8_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      uint8_decompress(array, buf, buf_size);
    }
   
  public:

    void compress_to_target_size(const int16_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      int16_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const int16_array2_t& array,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      int16_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }
    
    void compress_to_target_amax_error(const int16_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      int16_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    void decompress(int16_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      int16_decompress(array, buf, buf_size);
    }

  public:
    
    void compress_to_target_size(const uint16_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      uint16_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const uint16_array2_t& array,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      uint16_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_amax_error(const uint16_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      uint16_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    void decompress(uint16_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      uint16_decompress(array, buf, buf_size);
    }
    
  public:

    void compress_to_target_size(const int32_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      int32_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const int32_array2_t& array,
                                double target_rms_error,
                                void* buf,
                                std::size_t buf_size,
                                std::size_t* actual_size,
                                double* actual_rms_error) {
      int32_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_amax_error(const int32_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      int32_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    void decompress(int32_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      int32_decompress(array, buf, buf_size);
    }

  public:
    
    void compress_to_target_size(const uint32_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      uint32_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const uint32_array2_t& array,
                                      double target_rms_error,
                                      void* buf,
                                      std::size_t buf_size,
                                      std::size_t* actual_size,
                                      double* actual_rms_error) {
      uint32_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_amax_error(const uint32_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      uint32_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }

    void decompress(uint32_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      uint32_decompress(array, buf, buf_size);
    }
    
  public:

    void compress_to_target_size(const float_array2_t& array,
                                 std::size_t target_size,
                                 void* buf,
                                 std::size_t buf_size,
                                 std::size_t* actual_size,
                                 double* actual_rms_error) {
      float_compress(array, target_size, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_rms_error(const float_array2_t& array,
                                      double target_rms_error,
                                      void* buf,
                                      std::size_t buf_size,
                                      std::size_t* actual_size,
                                      double* actual_rms_error) {
      float_compress(array, buf_size, target_rms_error, buf, buf_size, actual_size, actual_rms_error);
    }

    void compress_to_target_amax_error(const float_array2_t& array,
                                       double target_amax_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t* actual_size,
                                       double* actual_amax_error) {
      float_compress_amax(array, buf_size, target_amax_error, buf, buf_size, actual_size, actual_amax_error);
    }
    
    void decompress(float_array2_t& array,
                    const void* buf,
                    std::size_t buf_size) {
      float_decompress(array, buf, buf_size);
    }

  public:

    template <class T>
    void compress(const dense_array<T, 2, void>& array,
                  void* buf,
                  std::size_t buf_size,
                  std::size_t* actual_size,
                  double* actual_rms_error) {
      compress_to_target_rms_error(array, 0.0, buf, buf_size, actual_size, actual_rms_error);
    }
    
  };

} // namespace sl

#endif
