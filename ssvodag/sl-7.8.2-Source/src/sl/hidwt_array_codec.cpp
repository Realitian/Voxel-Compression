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
#include <sl/hidwt_array_codec.hpp>
#include <sl/fixed_size_square_matrix.hpp>

#define FIXME_NO_CLAMPING 0 
#define FIXME_NO_RESCALE  0 

namespace sl {

  // -----------------------------------------------------------------
  // Misc helpers

  typedef hidwt_array_codec::error_kind_t hidwt_error_kind_t;
  
  static inline int32_t jls_quantize(int32_t x, uint32_t eps, uint32_t eps_times2_plus1) {
    return 
      eps ?
      ((x>=0) ? (x+int32_t(eps))/int32_t(eps_times2_plus1) : -(int32_t(eps)-x)/int32_t(eps_times2_plus1))
      :
      (x);
  }

  static inline int32_t jls_dequantize(int32_t x, uint32_t eps_times2_plus1) {
    return x*int32_t(eps_times2_plus1);
  }

#if 1
  static inline int32_t pow2quantize(int32_t x, uint32_t log2_threshold) {
    return x>=0 ? (x>>int32_t(log2_threshold)) : -((-x)>>int32_t(log2_threshold));
  }

  static inline int32_t pow2dequantize(int32_t x, uint32_t log2_threshold) {
    const int32_t half_interval = log2_threshold ? ((int32_t(1)<<int32_t(log2_threshold-1))-1) : 0;
    return
      x
      ?
      ((x>0) ? ((x<<int32_t(log2_threshold)) + half_interval) : -(((-x)<<int32_t(log2_threshold)) + half_interval))
      :
      0;
  }
#else
  static inline int32_t pow2quantize(int32_t x, uint32_t log2_threshold) {
    const uint32_t eps = (1<<log2_threshold)-1;
    const uint32_t eps_times2_plus1 = eps*2+1;
    return jls_quantize(x, eps, eps_times2_plus1);
  }
  
  static inline int32_t pow2dequantize(int32_t x, uint32_t log2_threshold) {
    const uint32_t eps = (1<<log2_threshold)-1;
    const uint32_t eps_times2_plus1 = eps*2+1;
    return jls_dequantize(x, eps_times2_plus1);
  }
#endif
  
#if 1
    // Precise but slower
    static inline int32_t rounded_half(int32_t x) {
      return x>=0 ? ((x+1)>>1) : -((-x+1)>>1);
    }

    static inline int32_t rounded_quarter(int32_t x) {
      return x>=0 ? ((x+2)>>2) : -((-x+2)>>2);
    }
#else
    // Biased, but faster
    static inline int32_t rounded_half(int32_t x) {
      return (x+1)>>1;
    }

    static inline int32_t rounded_quarter(int32_t x) {
      return (x+2)>>2;
    }
#endif

  template <std::size_t sz>
  static inline int32_t iamax(const fixed_size_square_matrix<sz, int32_t>& M) {
    int32_t result = 0;
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	result = std::max(result, abs(M(i,j)));
      }
    }
    return result;
  } 

  template <std::size_t sz>
  static inline double iamax(const fixed_size_square_matrix<sz, int32_t>& M,
			    const fixed_size_square_matrix<sz, int32_t>& M_prime) {
    double result = 0.0;
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	result = std::max(result, sl::abs(double(M(i,j))-double(M_prime(i,j))));
      }
    }
    return result;
  }

  template <std::size_t sz>
  static inline double irms(const fixed_size_square_matrix<sz, int32_t>& M,
			    const fixed_size_square_matrix<sz, int32_t>& M_prime) {
    double result = 0.0;
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	double e_ij = M_prime(i,j)-M(i,j);
	result += e_ij*e_ij;
      }
    }
    return std::sqrt(result/double(sz*sz));
  }

  template <std::size_t sz>
  static inline double ierror(const fixed_size_square_matrix<sz, int32_t>& M,
                              const fixed_size_square_matrix<sz, int32_t>& M_prime,
                              hidwt_error_kind_t error_kind) {
    return (error_kind == hidwt_array_codec::Error_kind_rms) ? irms(M,M_prime) : iamax(M,M_prime);
  } 
  
  template <std::size_t sz>
  static inline void clamp(fixed_size_square_matrix<sz, int32_t>& M,
			   int32_t M_min,
			   int32_t M_max) {
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	M(i,j) = median(M_min,M(i,j),M_max);
      }
    }
  }

  template <std::size_t sz>
  static inline void pow2quantize_in(fixed_size_square_matrix<sz, int32_t>& QM,
				     const fixed_size_square_matrix<sz, int32_t>& M,
				     uint32_t log2_threshold) {
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	QM(i,j) = pow2quantize(M(i,j), log2_threshold);
      }
    }
  }

  template <std::size_t sz>
  static inline void pow2dequantize_in(fixed_size_square_matrix<sz, int32_t>& M,
				    const fixed_size_square_matrix<sz, int32_t>& QM,
				    uint32_t log2_threshold) {
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	M(i,j) = pow2dequantize(QM(i,j), log2_threshold);
      }
    }
  }

  template <std::size_t sz>
  void rescaled_pow2quantize_in(fixed_size_square_matrix<sz, int32_t>& QM,
				const fixed_size_square_matrix<sz, int32_t>& M,
				uint32_t log2_threshold) {
#if FIXME_NO_RESCALE
    pow2quantize_in(QM, M, log2_threshold);
#else

    int32_t LL_scale = 2; // 2^2 = 4
    int32_t LH_scale = 1; // 2^1 = 2
    int32_t HL_scale = 1; // 2^1 = 2
    int32_t HH_scale = 0; // 2^0 = 1
    for (uint32_t level_sz = sz; level_sz>=8; level_sz>>=1) {
      uint32_t level_half_sz = level_sz>>1;

      uint32_t LL_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - LL_scale);
      uint32_t LH_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - LH_scale);
      uint32_t HL_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - HL_scale);
      uint32_t HH_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - HH_scale);

      for (uint32_t i=0; i<level_half_sz; ++i) {
	const uint32_t i_l = i, i_h = i+level_half_sz;
	for (uint32_t j=0; j<level_half_sz; ++j) {
	  const uint32_t j_l = j, j_h = j+level_half_sz;
	  QM(i_l, j_l) = pow2quantize(M(i_l,j_l), LL_log2_threshold);
	  QM(i_l, j_h) = pow2quantize(M(i_l,j_h), LH_log2_threshold);
	  QM(i_h, j_l) = pow2quantize(M(i_h,j_l), HL_log2_threshold);
	  QM(i_h, j_h) = pow2quantize(M(i_h,j_h), HH_log2_threshold);
	}
      }
      const int32_t old_LL_scale = LL_scale;
      LL_scale = old_LL_scale + 1;
      LH_scale = old_LL_scale + 0;
      HL_scale = old_LL_scale + 0;
      HH_scale = old_LL_scale - 1;
    }
#endif
  }

  template <std::size_t sz>
  void rescaled_pow2dequantize_in(fixed_size_square_matrix<sz, int32_t>& M,
				  const fixed_size_square_matrix<sz, int32_t>& QM,
				  uint32_t log2_threshold) {

#if FIXME_NO_RESCALE
    pow2dequantize_in(M, QM, log2_threshold);
#else
    int32_t LL_scale = 2; // 2^2 = 4
    int32_t LH_scale = 1; // 2^1 = 2
    int32_t HL_scale = 1; // 2^1 = 2
    int32_t HH_scale = 0; // 2^0 = 1
    for (uint32_t level_sz = sz; level_sz>=8; level_sz>>=1) {
      uint32_t level_half_sz = level_sz>>1;

      uint32_t LL_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - LL_scale);
      uint32_t LH_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - LH_scale);
      uint32_t HL_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - HL_scale);
      uint32_t HH_log2_threshold = std::max(int32_t(0), int32_t(log2_threshold) - HH_scale);

      for (uint32_t i=0; i<level_half_sz; ++i) {
	const uint32_t i_l = i, i_h = i+level_half_sz;
	for (uint32_t j=0; j<level_half_sz; ++j) {
	  const uint32_t j_l = j, j_h = j+level_half_sz;
	  M(i_l, j_l) = pow2dequantize(QM(i_l,j_l), LL_log2_threshold);
	  M(i_l, j_h) = pow2dequantize(QM(i_l,j_h), LH_log2_threshold);
	  M(i_h, j_l) = pow2dequantize(QM(i_h,j_l), HL_log2_threshold);
	  M(i_h, j_h) = pow2dequantize(QM(i_h,j_h), HH_log2_threshold);
	}
      }
      const int32_t old_LL_scale = LL_scale;
      LL_scale = old_LL_scale + 1;
      LH_scale = old_LL_scale + 0;
      HL_scale = old_LL_scale + 0;
      HH_scale = old_LL_scale - 1;
    }
#endif
  }

  // -----------------------------------------------------------------

  template <std::size_t sz>
  static inline void jls_quantize_in(fixed_size_square_matrix<sz, int32_t>& QM,
				     const fixed_size_square_matrix<sz, int32_t>& M,
				     uint32_t threshold) {
    uint32_t threshold_times2_plus1 = threshold*2+1;
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	QM(i,j) = jls_quantize(M(i,j), threshold, threshold_times2_plus1);
      }
    }
  }

  template <std::size_t sz>
  static inline void jls_dequantize_in(fixed_size_square_matrix<sz, int32_t>& M,
				    const fixed_size_square_matrix<sz, int32_t>& QM,
				    uint32_t threshold) {
    uint32_t threshold_times2_plus1 = threshold*2+1;
    for (uint32_t i=0; i<sz; ++i) {
      for (uint32_t j=0; j<sz; ++j) {
	M(i,j) = jls_dequantize(QM(i,j), threshold_times2_plus1);
      }
    }
  }

  template <std::size_t sz>
  void rescaled_jls_quantize_in(fixed_size_square_matrix<sz, int32_t>& QM,
				const fixed_size_square_matrix<sz, int32_t>& M,
				uint32_t threshold) {
#if FIXME_NO_RESCALE
    jls_quantize_in(QM, M, threshold);
#else

    uint32_t LL_scale = 4; // 2^2 = 4
    uint32_t LH_scale = 2; // 2^1 = 2
    uint32_t HL_scale = 2; // 2^1 = 2
    uint32_t HH_scale = 1; // 2^0 = 1
    for (uint32_t level_sz = sz; level_sz>=8; level_sz>>=1) {
      uint32_t level_half_sz = level_sz>>1;

      uint32_t LL_threshold = threshold/LL_scale; uint32_t LL_threshold_times2_plus1 = LL_threshold*2+1;
      uint32_t LH_threshold = threshold/LH_scale; uint32_t LH_threshold_times2_plus1 = LH_threshold*2+1;
      uint32_t HL_threshold = threshold/HL_scale; uint32_t HL_threshold_times2_plus1 = HL_threshold*2+1;
      uint32_t HH_threshold = threshold/HH_scale; uint32_t HH_threshold_times2_plus1 = HH_threshold*2+1;

      for (uint32_t i=0; i<level_half_sz; ++i) {
	const uint32_t i_l = i, i_h = i+level_half_sz;
	for (uint32_t j=0; j<level_half_sz; ++j) {
	  const uint32_t j_l = j, j_h = j+level_half_sz;
	  QM(i_l, j_l) = jls_quantize(M(i_l,j_l), LL_threshold, LL_threshold_times2_plus1);
	  QM(i_l, j_h) = jls_quantize(M(i_l,j_h), LH_threshold, LH_threshold_times2_plus1);
	  QM(i_h, j_l) = jls_quantize(M(i_h,j_l), HL_threshold, HL_threshold_times2_plus1);
	  QM(i_h, j_h) = jls_quantize(M(i_h,j_h), HH_threshold, HH_threshold_times2_plus1);
	}
      }
      const int32_t old_LL_scale = LL_scale;
      LL_scale = old_LL_scale/2;
      LH_scale = old_LL_scale;
      HL_scale = old_LL_scale;
      HH_scale = old_LL_scale*2;
    }
#endif
  }

  template <std::size_t sz>
  void rescaled_jls_dequantize_in(fixed_size_square_matrix<sz, int32_t>& M,
				  const fixed_size_square_matrix<sz, int32_t>& QM,
				  uint32_t threshold) {

#if FIXME_NO_RESCALE
    jls_dequantize_in(M, QM, threshold);
#else
    uint32_t LL_scale = 4; // 2^2 = 4
    uint32_t LH_scale = 2; // 2^1 = 2
    uint32_t HL_scale = 2; // 2^1 = 2
    uint32_t HH_scale = 1; // 2^0 = 1
    for (uint32_t level_sz = sz; level_sz>=8; level_sz>>=1) {
      uint32_t level_half_sz = level_sz>>1;

      uint32_t LL_threshold = threshold/LL_scale; uint32_t LL_threshold_times2_plus1 = LL_threshold*2+1;
      uint32_t LH_threshold = threshold/LH_scale; uint32_t LH_threshold_times2_plus1 = LH_threshold*2+1;
      uint32_t HL_threshold = threshold/HL_scale; uint32_t HL_threshold_times2_plus1 = HL_threshold*2+1;
      uint32_t HH_threshold = threshold/HH_scale; uint32_t HH_threshold_times2_plus1 = HH_threshold*2+1;

      for (uint32_t i=0; i<level_half_sz; ++i) {
	const uint32_t i_l = i, i_h = i+level_half_sz;
	for (uint32_t j=0; j<level_half_sz; ++j) {
	  const uint32_t j_l = j, j_h = j+level_half_sz;
	  M(i_l, j_l) = jls_dequantize(QM(i_l,j_l), LL_threshold_times2_plus1);
	  M(i_l, j_h) = jls_dequantize(QM(i_l,j_h), LH_threshold_times2_plus1);
	  M(i_h, j_l) = jls_dequantize(QM(i_h,j_l), HL_threshold_times2_plus1);
	  M(i_h, j_h) = jls_dequantize(QM(i_h,j_h), HH_threshold_times2_plus1);
	}
      }
      const int32_t old_LL_scale = LL_scale;
      LL_scale = old_LL_scale/2;
      LH_scale = old_LL_scale;
      HL_scale = old_LL_scale;
      HH_scale = old_LL_scale*2;
    }
#endif
  }

  // -----------------------------------------------------------------

  template <class xform_t, std::size_t sz>
  static inline double pow2quantized_xform_error(const xform_t& xform,
                                                 const fixed_size_square_matrix<sz, int32_t>& M,
                                                 uint32_t log2_threshold,
                                                 hidwt_error_kind_t error_kind) {
    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;

    int32_matrix_t Mw= M, QMw, M_prime; //QMw= tags::not_initialized(), M_prime= tags::not_initialized();
    xform.forward(Mw);
    rescaled_pow2quantize_in(QMw, Mw, log2_threshold);
    rescaled_pow2dequantize_in(M_prime, QMw, log2_threshold); 
    xform.backward(M_prime); 
#if !FIXME_NO_CLAMPING
    int32_t M_amax = iamax(M);
    clamp(M_prime, -M_amax, M_amax);
#endif
   return ierror(M, M_prime, error_kind);
  }

  template <class xform_t, std::size_t sz>
  static inline double jls_quantized_xform_error(const xform_t& xform,
                                                 const fixed_size_square_matrix<sz, int32_t>& M,
                                                 uint32_t threshold,
                                                 hidwt_error_kind_t error_kind) {
    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;

    int32_matrix_t Mw= M, QMw,M_prime; //QMw= tags::not_initialized(), M_prime= tags::not_initialized();
    xform.forward(Mw);
    rescaled_jls_quantize_in(QMw, Mw, threshold);
    rescaled_jls_dequantize_in(M_prime, QMw, threshold); 
    xform.backward(M_prime); 
#if !FIXME_NO_CLAMPING
    int32_t M_amax = iamax(M);
    clamp(M_prime, -M_amax, M_amax);
#endif
   return ierror(M, M_prime, error_kind);
  }

  // -----------------------------------------------------------------
  
  template <uint32_t log2_sz>
  class hidwt_matrix_transform {
  public:
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };

    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    typedef fixed_size_square_matrix<half_sz, int32_t> int32_half_size_matrix_t;
    
    static inline void decompose(const int32_matrix_t& M,
                                 int32_half_size_matrix_t& LL,
                                 int32_half_size_matrix_t& LH,
                                 int32_half_size_matrix_t& HL,
                                 int32_half_size_matrix_t& HH) {
      for (uint32_t i=0; i<half_sz; ++i) {
        for (uint32_t j=0; j<half_sz; ++j) {
          LL(i, j) = M(i, j);
          LH(i, j) = M(i, j+half_sz);
          HL(i, j) = M(i+half_sz, j);
          HH(i, j) = M(i+half_sz, j+half_sz);
        }
      }
    }
    
    static inline void recompose(int32_matrix_t& M,
                                 const int32_half_size_matrix_t& LL,
                                 const int32_half_size_matrix_t& LH,
                                 const int32_half_size_matrix_t& HL,
                                 const int32_half_size_matrix_t& HH) {
      for (uint32_t i=0; i<half_sz; ++i) {
        for (uint32_t j=0; j<half_sz; ++j) {
          M(i, j) = LL(i, j);
          M(i, j+half_sz) = LH(i, j);
          M(i+half_sz, j) = HL(i, j);
          M(i+half_sz, j+half_sz) = HH(i, j);
        }
      }
    }

  };

  // -----------------------------------------------------------------
    /**
   *  Transforms M using one step of the integer Haar wavelet transformation.
   */
  template <uint32_t log2_sz>
  class hidwt_haar_matrix_transform: public hidwt_matrix_transform<log2_sz> {
  public:
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };

    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    
    static inline void forward(int32_matrix_t& M) {
      int32_t l[half_sz];
      int32_t h[half_sz];

      // Levels
      for (uint32_t level_sz = sz; level_sz>=8; level_sz>>=1) {
	const uint32_t level_half_sz = level_sz>>1;

	// Rows
	for (uint32_t k=0; k<level_sz; ++k) {
	  // H
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    h[i] = M(k, 2*i+1) - M(k, 2*i);
	  }
	  // L
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    l[i] = M(k, 2*i) + rounded_half(h[i]);
	  }
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(k, i)               = l[i];
	    M(k, i+level_half_sz) = h[i];
	  }
	}

	// Columns
	for (uint32_t k=0; k<level_sz; ++k) {
	  // H
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    h[i] = M(2*i+1, k) - M(2*i, k);
	  }
	  // L
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    l[i] = M(2*i, k) + rounded_half(h[i]);
	  }
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(i, k)         = l[i];
	    M(i+level_half_sz, k) = h[i];
	  }
	}
      }
    }

    static inline void backward(int32_matrix_t& M) {
      int32_t l[half_sz];
      int32_t h[half_sz];
      
      for (uint32_t level_sz = 8; level_sz<=sz; level_sz<<=1) {
	const uint32_t level_half_sz = level_sz>>1;
	
	// Columns
	for (uint32_t k=0; k<level_sz; ++k) {
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    l[i] = M(i, k);
	    h[i] = M(i+level_half_sz, k);
	  }

	  // Even
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(2*i, k) = l[i] - rounded_half(h[i]);
	  }
	  
	  // Odd
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(2*i+1, k) = h[i] + M(2*i, k);
	  }
	}
	
	// Rows
	for (uint32_t k=0; k<level_sz; ++k) {
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    l[i] = M(k, i);
	    h[i] = M(k, i+level_half_sz);
	  }
	  
	  // Even
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(k, 2*i) = l[i] - rounded_half(h[i]);
	  }
	  
	  // Odd
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(k, 2*i+1) = h[i] + M(k, 2*i);
	  }
	}
      }
    }
  };

  /**
   *  Transforms M using one step of the (2,2)  integer wavelet transformation of
   *  Calderbank, Daubechies, Sweldens, Yeo: Lossless image compression using
   *  wavelet transforms that map integers to integers.
   *  This wavelet is a simple linear interpolation transform that predicts
   *  that an odd element will be on a line between its two even neighbors.
   *  The difference between this prediction and the actual value becomes
   *  the detail coefficient, while the scaling function is simply the
   *  average of the original even and odd elements.
   */
  template <uint32_t log2_sz>
  class hidwt_c22_matrix_transform: public hidwt_matrix_transform<log2_sz> {
  public:
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };

    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    
    static inline void forward(int32_matrix_t& M) {
      
      int32_t l[half_sz];
      int32_t h[half_sz];
      
      // Levels
      for (uint32_t level_sz = sz; level_sz>=8; level_sz>>=1) {
	const uint32_t level_half_sz = level_sz>>1;

	// Rows
	for (uint32_t k=0; k<level_sz; ++k) {
	  // H
	  for (uint32_t i=0; i<level_half_sz-1; ++i) {
	    h[i] = M(k, 2*i+1) - rounded_half(M(k, 2*i)+M(k, 2*i+2));
	  }
	  h[level_half_sz-1] = M(k, level_sz-1) - rounded_half(3*M(k, level_sz-2)-M(k, level_sz-4)); // Linear extrapolation
	  
	  // L
	  l[0] = M(k, 0) + rounded_quarter(2*h[0]-h[1]+h[0]);
	  M(k, 0)     = l[0];
	  for (uint32_t i=1; i<level_half_sz; ++i) {
	    l[i] = M(k, 2*i) + rounded_quarter(h[i-1]+h[i]);
	    M(k, i) = l[i];
	  }
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(k, i+level_half_sz) = h[i];
	  }
	}
	
	// Columns
	for (uint32_t k=0; k<level_sz; ++k) {
	  // H
	  for (uint32_t i=0; i<level_half_sz-1; ++i) {
	    h[i] = M(2*i+1, k) - rounded_half(M(2*i, k)+M(2*i+2, k));
	  }
	  h[level_half_sz-1] = M(level_sz-1, k) - rounded_half(3*M(level_sz-2, k)-M(level_sz-4, k)); // Linear extrapolation
	  
	  // L
	  l[0] = M(0, k) + rounded_quarter(2*h[0]-h[1]+h[0]);
	  M(0, k)         = l[0];
	  for (uint32_t i=1; i<level_half_sz; ++i) {
	    l[i] = M(2*i, k) + rounded_quarter(h[i-1]+h[i]);
	    M(i, k)         = l[i];
	  }
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    M(i+level_half_sz, k) = h[i];
	  }
	}
      }
    }

    static inline void backward(int32_matrix_t& M) {
      int32_t l[half_sz];
      int32_t h[half_sz];
      
      for (uint32_t level_sz = 8; level_sz<=sz; level_sz<<=1) {
	const uint32_t level_half_sz = level_sz>>1;
	// Columns
	for (uint32_t k=0; k<level_sz; ++k) {
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    l[i] = M(i, k);
	    h[i] = M(i+level_half_sz, k);
	  }
	  
	  // Even
	  M(0, k) = l[0] - rounded_quarter(2*h[0]-h[1]+h[0]);
	  for (uint32_t i=1; i<level_half_sz; ++i) {
	    M(2*i, k) = l[i] - rounded_quarter(h[i-1]+h[i]);
	  }

	  // Odd
	  for (uint32_t i=0; i<level_half_sz-1; ++i) {
	    M(2*i+1, k) = h[i] + rounded_half(M(2*i, k)+M(2*i+2, k));
	  }
	  M(level_sz-1, k) = h[level_half_sz-1] + rounded_half(3*M(level_sz-2, k)-M(level_sz-4, k)); // Linear extrapolation
	}
	
	// ROWS
	for (uint32_t k=0; k<level_sz; ++k) {
	  for (uint32_t i=0; i<level_half_sz; ++i) {
	    l[i] = M(k, i);
	    h[i] = M(k, i+level_half_sz);
	  }
	  
	  // Even
	  M(k, 0) = l[0] - rounded_quarter(2*h[0]-h[1]+h[0]);
	  for (uint32_t i=1; i<level_half_sz; ++i) {
	    M(k, 2*i) = l[i] - rounded_quarter(h[i-1]+h[i]);
	  }
	  
	  // Odd
	  for (uint32_t i=0; i<level_half_sz-1; ++i) {
	    M(k, 2*i+1) = h[i] + rounded_half(M(k, 2*i)+M(k, 2*i+2));
	  }
	  M(k, level_sz-1) = h[level_half_sz-1] + rounded_half(3*M(k, level_sz-2)-M(k, level_sz-4)); // Linear extrapolation
	}
      }
    }
  }; 

  // -----------------------------------------------------------------
  
  /**
   * Compressor: quadtree encoding of a quantized matrix
   */
  template <class entropy_codec, uint32_t log2_sz>
  class hidwt_quantized_matrix_codec {
  public:
    typedef entropy_codec entropy_codec_t;
    typedef typename entropy_codec_t::bit_context_t    bit_context_t;
    typedef typename entropy_codec_t::int_context_t    int_context_t;
    typedef typename entropy_codec_t::symbol_context_t symbol_context_t;
    
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };

    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    typedef fixed_size_square_matrix<half_sz, int32_t> int32_half_size_matrix_t;

    typedef hidwt_matrix_transform<log2_sz>                                        matrix_transform_t;
    typedef hidwt_matrix_transform<log2_sz-1>                                      half_size_matrix_transform_t;
    typedef hidwt_quantized_matrix_codec<entropy_codec, log2_sz>                   matrix_codec_t;
    typedef hidwt_quantized_matrix_codec<entropy_codec, log2_sz-1>                 half_size_matrix_codec_t;
    
  public:
     
    static inline void encode(entropy_codec_t& codec,
                              const int32_matrix_t& QM) {
      // Quadtree split
      bool is_null_matrix = true;
      for (uint32_t i=0; (i<sz) && is_null_matrix; ++i) {
        for (uint32_t j=0; (j<sz) && is_null_matrix; ++j) {
          is_null_matrix = (QM(i, j) == 0);
        }
      }

      codec.encode_bit(is_null_matrix);
      if (is_null_matrix) {
        // Nothing to do
        //std::cerr << "Encode NULL " << sz << "x" << sz << std::endl;
     } else {
        // Nonzero matrix, quadtree split
        //std::cerr << "Encode Recurse " << sz << "x" << sz << std::endl;
        int32_half_size_matrix_t QLL, QLH, QHL, QHH;
        matrix_transform_t::decompose(QM, QLL, QLH, QHL, QHH);
        half_size_matrix_codec_t::encode(codec, QLL);
        half_size_matrix_codec_t::encode(codec, QLH);
        half_size_matrix_codec_t::encode(codec, QHL);
        half_size_matrix_codec_t::encode(codec, QHH);
      }
    }
        
    static inline void decode(entropy_codec_t& codec,
                              int32_matrix_t& QM) {
      bool is_null_matrix = codec.decode_bit();
      if (is_null_matrix) {
        //std::cerr << "Decode NULL " << sz << "x" << sz << std::endl;
	
        // Fill with zeros
        for (uint32_t i=0; i<sz; ++i) {
          for (uint32_t j=0; j<sz; ++j) {
            QM(i, j) = 0;
          }
        }
      } else {
        //std::cerr << "Decode Recurse " << sz << "x" << sz << std::endl;

        // Nonzero matrix, quadtree split
        int32_half_size_matrix_t QLL, QLH, QHL, QHH;
        half_size_matrix_codec_t::decode(codec, QLL);
        half_size_matrix_codec_t::decode(codec, QLH);
        half_size_matrix_codec_t::decode(codec, QHL);
        half_size_matrix_codec_t::decode(codec, QHH);
        matrix_transform_t::recompose(QM, QLL, QLH, QHL, QHH);
      }
    }

  }; // hidwt_quantized_matrix_codec

  // -----------------------------------------------------------------

  /**
   * Codec specialization for small size (2x2)
   */
  template <class entropy_codec>
  class hidwt_quantized_matrix_codec<entropy_codec, 1> {
  public:
    enum { log2_sz = 1 };
    typedef entropy_codec entropy_codec_t;
    typedef typename entropy_codec_t::bit_context_t    bit_context_t;
    typedef typename entropy_codec_t::int_context_t    int_context_t;
    typedef typename entropy_codec_t::symbol_context_t symbol_context_t;
    
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };

    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    typedef fixed_size_square_matrix<half_sz, int32_t> int32_half_size_matrix_t;

    typedef hidwt_matrix_transform<log2_sz>                                  matrix_transform_t;
    typedef hidwt_matrix_transform<log2_sz-1>                                half_size_matrix_transform_t;
    typedef hidwt_quantized_matrix_codec<entropy_codec, log2_sz>             matrix_codec_t;
    typedef hidwt_quantized_matrix_codec<entropy_codec, log2_sz-1>           half_size_matrix_codec_t;
    
  public:
    
    static inline void encode(entropy_codec_t& codec,
                              const int32_matrix_t& QM) {
      bool is_null_matrix = true;
      for (uint32_t i=0; (i<sz) && is_null_matrix; ++i) {
        for (uint32_t j=0; (j<sz) && is_null_matrix; ++j) {
          is_null_matrix = (QM(i, j) == 0);
        }
      }
      
      codec.encode_bit(is_null_matrix);
      if (is_null_matrix) {
        // Nothing to do
        //std::cerr << "  Encode NULL " << sz << "x" << sz << std::endl;
      } else {
        //std::cerr << "  Encode: " << sz << "x" << sz << ": ";
        for (uint32_t i=0; i<sz; ++i) {
          for (uint32_t j=0; j<sz; ++j) {
            codec.encode_int(QM(i, j));
	    //std::cerr << QM(i, j) << " ";
          }
        }
	//std::cerr << std::endl;
      }
    }
        
    static inline void decode(entropy_codec_t& codec,
                              int32_matrix_t& QM) {
      bool is_null_matrix = codec.decode_bit();
      if (is_null_matrix) {
        //std::cerr << "Decode NULL " << sz << "x" << sz << std::endl;
        // Fill with zeros
        for (uint32_t i=0; i<sz; ++i) {
          for (uint32_t j=0; j<sz; ++j) {
            QM(i, j) = 0;
          }
        }
      } else {
        //std::cerr << "Decode: " << sz << "x" << sz << ": ";
        for (uint32_t i=0; i<sz; ++i) {
          for (uint32_t j=0; j<sz; ++j) {
            QM(i, j) = codec.decode_int();
	    //std::cerr << QM(i, j) << " ";
          }
        }
	//std::cerr << std::endl;
      }
    }

  }; // hidwt_quantized_matrix_codec
  
  // -----------------------------------------------------------------

  /**
   * Compressor dispatcher on size
   */
  template <class entropy_codec, class int_t, uint32_t log2_sz>
  class hidwt_c22_matrix_codec {
  public:
    typedef entropy_codec entropy_codec_t;
    typedef typename entropy_codec_t::bit_context_t    bit_context_t;
    typedef typename entropy_codec_t::int_context_t    int_context_t;
    typedef typename entropy_codec_t::symbol_context_t symbol_context_t;
    
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };
    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    typedef fixed_size_square_matrix<half_sz, int32_t> int32_half_size_matrix_t;
#if 1
    typedef hidwt_c22_matrix_transform<log2_sz>  matrix_transform_t;
#else
    typedef hidwt_haar_matrix_transform<log2_sz> matrix_transform_t;
#endif
    
    typedef hidwt_quantized_matrix_codec<entropy_codec, log2_sz> quantized_matrix_codec_t;
    typedef hidwt_c22_matrix_codec<entropy_codec, int_t, log2_sz+1>     double_size_matrix_codec_t;

    typedef dense_array<int_t, 2, void> int_matrix_t;
    
  public:

    static inline int32_t average(const int_matrix_t& M) {
      double M_sum = 0.0;
	  uint32_t isz = static_cast<uint32_t>(M.extent()[0]);
	  uint32_t jsz = static_cast<uint32_t>(M.extent()[1]);
      for (uint32_t i=0; i<isz; ++i) {
	for (uint32_t j=0; j<jsz; ++j) {
	  M_sum += double(M(i,j));
	}
      }
      return int32_t(M_sum > 0.0f ? (M_sum/double(isz*jsz)+0.5) : (M_sum/double(isz*jsz)-0.5));
    }
         
    static inline void compress(entropy_codec_t& codec,
                                const int_matrix_t& M,
                                uint32_t log2_threshold) {
	  uint32_t isz = static_cast<uint32_t>(M.extent()[0]);
	  uint32_t jsz = static_cast<uint32_t>(M.extent()[1]);
      if (sz>=isz && sz>=jsz) {
        // Copy and recenter
        int32_t offset = average(M);
        
        int32_matrix_t  MM; //MM= tags::not_initialized();
	int32_t         MM_amax = 0;
        for (uint32_t i=0; i<sz; ++i) {
          const uint32_t ii = (i<isz) ? i : (isz-2-(i-isz));
          for (uint32_t j=0; j<sz; ++j) {
            const uint32_t jj = (j<jsz) ? j : (jsz-2-(j-jsz));
            int32_t M_ij = int32_t(M(ii,jj))-offset;
	    MM(i, j) = M_ij;
	    MM_amax  =  sl::max(MM_amax, abs(M_ij));
          }
        }
	
        codec.encode_int(offset);  // Center of original data
	codec.encode_int(MM_amax); // Amax of recentered data
	//std::cerr << "ENCODED: offset=" << offset << " amax = " << MM_amax << std::endl;
        compress(codec, MM, log2_threshold);
      } else {
        double_size_matrix_codec_t::compress(codec, M, log2_threshold);
      }
    }
    
    static inline double error(const int_matrix_t& M,
                               uint32_t log2_threshold,
                               hidwt_error_kind_t error_kind) {
      uint32_t isz = static_cast<uint32_t>(M.extent()[0]);
	  uint32_t jsz = static_cast<uint32_t>(M.extent()[1]);
      if (sz>=isz && sz>=jsz) {
        int32_t offset = average(M);

        int32_matrix_t  MM; //MM= tags::not_initialized();
        for (uint32_t i=0; i<sz; ++i) {
          const uint32_t ii = (i<isz) ? i : (isz-2-(i-isz));
          for (uint32_t j=0; j<sz; ++j) {
            const uint32_t jj = (j<jsz) ? j : (jsz-2-(j-jsz));
            MM(i, j) = int32_t(M(ii,jj))-offset;
          }
        }
        return pow2quantized_xform_error(matrix_transform_t(), MM, log2_threshold, error_kind);
      } else {
        return double_size_matrix_codec_t::error(M, log2_threshold, error_kind);
      }
    }
    
    static uint32_t log2_threshold_from_target_error(const int_matrix_t& M,
                                                     double target_error,
                                                     hidwt_error_kind_t error_kind) {
      
	  uint32_t isz = static_cast<uint32_t>(M.extent()[0]);
	  uint32_t jsz = static_cast<uint32_t>(M.extent()[1]);
      if (sz>=isz && sz>=jsz) {
        int32_t offset = average(M);
        int32_matrix_t  MM; //MM= tags::not_initialized();
        for (uint32_t i=0; i<sz; ++i) {
          const uint32_t ii = (i<isz) ? i : (isz-2-(i-isz));
          for (uint32_t j=0; j<sz; ++j) {
            const uint32_t jj = (j<jsz) ? j : (jsz-2-(j-jsz));
            MM(i, j) = int32_t(M(ii,jj))-offset;
          }
        }
        return log2_threshold_from_target_error(MM, target_error, error_kind);
      } else {
        return double_size_matrix_codec_t::log2_threshold_from_target_error(M, target_error, error_kind);
      }
    }

    static uint32_t log2_threshold_from_target_error(const int32_matrix_t& M,
                                                     double target_error,
                                                     hidwt_error_kind_t error_kind) {
      // Dmax
#if !FIXME_NO_CLAMPING
      int32_t M_amax =iamax(M);
#endif
      // Wavelet xform
      int32_matrix_t Mw= M, QMw, M_prime; //QMw= tags::not_initialized(), M_prime= tags::not_initialized();
      matrix_transform_t::forward(Mw);

      // Search quantization space
      std::size_t steps = 0;
      // Coarsest possible value
      uint32_t log2_threshold_coarse = sl::median(uint32_t(0), uint32_t(sl::bitops<int32_t>::log2(iamax(Mw))+log2_sz), uint32_t(30));
      //std::cerr << steps << ": " << "T=" << log2_threshold_coarse << std::endl;
      if (log2_threshold_coarse == 0) return log2_threshold_coarse; // null matrix

      rescaled_pow2quantize_in(QMw, Mw, log2_threshold_coarse);
      rescaled_pow2dequantize_in(M_prime, QMw, log2_threshold_coarse);
      matrix_transform_t::backward(M_prime);
#if !FIXME_NO_CLAMPING
      clamp(M_prime, -M_amax, M_amax);
#endif
      double   error_coarse = ierror(M, M_prime, error_kind);
      ++steps;
      //std::cerr << steps << ": " << "T=" << log2_threshold_coarse << " => E= " << error_coarse << std::endl;
      if (error_coarse<=target_error) return log2_threshold_coarse;

      // Finest possible value
      uint32_t log2_threshold_fine   = 0; // 2^0 = 1 
      rescaled_pow2quantize_in(QMw, Mw, log2_threshold_fine);
      rescaled_pow2dequantize_in(M_prime, QMw, log2_threshold_fine);
      matrix_transform_t::backward(M_prime);
#if !FIXME_NO_CLAMPING
      clamp(M_prime, -M_amax, M_amax);
#endif
      double   error_fine = ierror(M, M_prime, error_kind);
      ++steps;
      //std::cerr << steps << ": " << "T=" << log2_threshold_fine << " => E= " << error_fine << std::endl;
      if (error_fine  >=target_error) return log2_threshold_fine;

      // Binary search
      assert(error_fine <= target_error);
      assert(target_error <= error_coarse);
      while (log2_threshold_coarse>log2_threshold_fine+1) {
        uint32_t log2_threshold_mid = (log2_threshold_fine+log2_threshold_coarse)/2;
        rescaled_pow2quantize_in(QMw, Mw, log2_threshold_mid);
        rescaled_pow2dequantize_in(M_prime, QMw, log2_threshold_mid);
        matrix_transform_t::backward(M_prime);
#if !FIXME_NO_CLAMPING
	clamp(M_prime, -M_amax, M_amax);
#endif
        double   error_mid = ierror(M, M_prime, error_kind);
        ++steps;
        //std::cerr << steps << ": " << "T=" << log2_threshold_mid << " => E= " << error_mid << std::endl;
        if (error_mid<=target_error) {
          log2_threshold_fine = log2_threshold_mid;
          error_fine = error_mid;
        } else {
          log2_threshold_coarse = log2_threshold_mid;
          error_coarse = error_mid;
        }
      }
#if 1
      // Conservative
      if (error_coarse <= 1.1*target_error) {
	return log2_threshold_coarse;
      } else {
	return log2_threshold_fine;
      }

#else
      // Now choose the 'best'
      if (sl::abs(error_coarse-target_error) <= sl::abs(error_fine-target_error)) {
        //std::cerr << steps << ": [CONVERGED] " << "T=" << log2_threshold_coarse << " => E= " << error_coarse << std::endl;
        return log2_threshold_coarse;
      } else {
        //std::cerr << steps << ": [CONVERGED] " << "T=" << log2_threshold_fine << " => E= " << error_fine << std::endl;
        return log2_threshold_fine;
      }
#endif
    }
      
    static inline void compress(entropy_codec_t& codec,
                                const int32_matrix_t& M,
                                uint32_t log2_threshold) {
      int32_matrix_t Mw=M, QMw; //QMw= tags::not_initialized();
      matrix_transform_t::forward(Mw);
      rescaled_pow2quantize_in(QMw, Mw, log2_threshold);
      quantized_matrix_codec_t::encode(codec, QMw);
    }
    
    static inline void decompress(entropy_codec_t& codec,
                                  int_matrix_t& M,
                                  uint32_t log2_threshold) {
	  uint32_t isz = static_cast<uint32_t>(M.extent()[0]);
	  uint32_t jsz = static_cast<uint32_t>(M.extent()[1]);
      if (sz>=isz && sz>=jsz) {
        int32_t offset  = codec.decode_int();
	int32_t MM_amax = codec.decode_int(); // amax of original data
	//std::cerr << "DECODED: offset=" << offset << " amax = " << MM_amax << std::endl;
 
        const int32_t int_min = int32_t(std::numeric_limits<int_t>::min());
        const int32_t int_max = int32_t(std::numeric_limits<int_t>::max());

        int32_matrix_t MM; //MM= tags::not_initialized();
        decompress(codec, MM, log2_threshold);
        for (uint32_t i=0; i<isz; ++i) {
          for (uint32_t j=0; j<jsz; ++j) {
#if FIXME_NO_CLAMPING
	    const int32_t MM_ij = MM(i, j);
#else
	    const int32_t MM_ij = median(-MM_amax, MM(i, j),  MM_amax);
#endif
            M(i,j) = int_t(median(MM_ij+offset, int_min, int_max));
          }
        }
      } else {
        double_size_matrix_codec_t::decompress(codec, M, log2_threshold);
      }
    }

    static inline void decompress(entropy_codec_t& codec,
                                  int32_matrix_t& M,
                                  uint32_t log2_threshold) {
      int32_matrix_t QMw; //QMw= tags::not_initialized();
      quantized_matrix_codec_t::decode(codec, QMw);
      rescaled_pow2dequantize_in(M, QMw, log2_threshold);
      matrix_transform_t::backward(M);
    }

  };
  
  /**
   * Compressor dispatcher on size - specialization
   * for maximum size matrices on which a full transform is performed.
   *
   * FIXME: Currently returns an error, should actually decompose the
   * matrix into blocks and compress the blocks separately.
   */
  template <class entropy_codec, class int_t>
  class hidwt_c22_matrix_codec<entropy_codec, int_t, 8> {
  public:
    enum { log2_sz = 8 };
    
    typedef entropy_codec entropy_codec_t;
    typedef typename entropy_codec_t::bit_context_t    bit_context_t;
    typedef typename entropy_codec_t::int_context_t    int_context_t;
    typedef typename entropy_codec_t::symbol_context_t symbol_context_t;
    
    enum { sz  = (uint32_t(1)<<log2_sz) };
    enum { half_sz = (sz>>1) };
    typedef fixed_size_square_matrix<sz, int32_t> int32_matrix_t;
    typedef fixed_size_square_matrix<half_sz, int32_t> int32_half_size_matrix_t;
        
    typedef dense_array<int_t, 2, void> int_matrix_t;
    
  public:

    static inline double error(const int_matrix_t& /*M*/,
                               uint32_t /*log2_threshold*/,
                               hidwt_error_kind_t /*error_kind*/) {
      SL_FAIL("Unsupported matrix size!");
      return 0;
    }
        
    static uint32_t log2_threshold_from_target_error(const int_matrix_t& /*array*/,
                                                     double /*target_error*/,
                                                     hidwt_error_kind_t /*error_kind*/) {
      SL_FAIL("Unsupported matrix size!");
      return 0;
    }
    
    static uint32_t log2_threshold_from_target_error(const int32_matrix_t& /*M*/,
                                                     double /*target_error*/,
                                                     hidwt_error_kind_t /*error_kind*/) {
      SL_FAIL("Unsupported matrix size!");
      return 0;
    }
      
    static inline void compress(entropy_codec_t& /*codec*/,
                                const int_matrix_t& /*M*/,
                                uint32_t /*threshold*/) {
      SL_FAIL("Unsupported matrix size!");
    }
    
    static inline void decompress(entropy_codec_t& /*codec*/,
                                  int_matrix_t& /*M*/,
                                  uint32_t /*threshold*/) {
      SL_FAIL("Unsupported matrix size!");
    }

  };
} // namespace sl

// ============================================================================
// hidwt_array_codec implementation
// ============================================================================

namespace sl {
  
  /**
   * Compress array to current codec buffer
   */
  template <class int_t>
  void hidwt_array_codec::int_quantized_compress(const dense_array<int_t,2,void>& array,
                                                 uint32_t log2_threshold,
                                                 std::size_t *actual_size, 
                                                 double *actual_error,
                                                 error_kind_t error_kind) {
    typedef hidwt_c22_matrix_codec<entropy_codec_t, int_t, 2> transform_codec_t;

    // Encode array shape
	const uint32_t h = static_cast<uint32_t>(array.extent()[0]);
	const uint32_t w = static_cast<uint32_t>(array.extent()[1]);

    entropy_codec_.start_encoder();
    {
      entropy_codec_.encode_int(h);
      entropy_codec_.encode_int(w);
      entropy_codec_.encode_int(log2_threshold);
      //std::cerr << "ENCODED: h: " << h << " w= " << w << " T= " << log2_threshold << std::endl;
      transform_codec_t::compress(entropy_codec_, array, log2_threshold);
    }
    entropy_codec_.stop_encoder();

    if (actual_size) {
      *actual_size = entropy_codec_.current_byte_count(); 
    }

    if (actual_error) {
      *actual_error = transform_codec_t::error(array, log2_threshold, error_kind);
    }
  }


  
  template <class int_t>
  void hidwt_array_codec::int_compress(const dense_array<int_t,2,void>& array,
                                       std::size_t target_size,
                                       double target_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t *actual_size, 
                                       double *actual_error,
                                       error_kind_t error_kind) {
#if 0
    //////////
    std::cerr << "----- compressing ------" << std::endl;
    for (std::size_t i=0; i< array.extent()[0]; ++i) {
      for (std::size_t j=0; j<array.extent()[1]; ++j) {
	std::cerr << std::setw(8) << array(i,j);
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
    /////////
#endif
    typedef hidwt_c22_matrix_codec<entropy_codec_t, int_t, 2> transform_codec_t;

    // Allocate buffer for compression trials
	const uint32_t max_buffer_size = static_cast<uint32_t>(1024 + 32 * array.count());
    if (entropy_codec_buffer_.size() < max_buffer_size) {
      entropy_codec_buffer_.resize(max_buffer_size);
      entropy_codec_.set_buffer(&(entropy_codec_buffer_[0]), max_buffer_size);
    }

    // Guess threshold from target error
    uint32_t log2_threshold_fine = transform_codec_t::log2_threshold_from_target_error(array, target_error, error_kind);

    // Compress and check size
    int_quantized_compress(array, log2_threshold_fine, actual_size, actual_error, error_kind);
    if (*actual_size <= target_size) {
      // Converged on error
    } else {
      // Buffer overflow!!
      std::cerr << "BUFFER OVERFLOW!!" << std::endl;
      
      uint32_t log2_threshold_coarse = std::max(log2_threshold_fine, uint32_t(30));
      if (log2_threshold_coarse>log2_threshold_fine) {
        std::size_t size_fine = *actual_size;
        std::size_t size_coarse;
        int_quantized_compress(array, log2_threshold_coarse, &size_coarse, NULL, error_kind);
        if (size_coarse>target_size) {
          // Buffer overflow!!
          SL_FAIL("Codec buffer overflow");
        } else {
          // Binary search: 
          while (log2_threshold_coarse>log2_threshold_fine+1) {
            uint32_t log2_threshold_mid = (log2_threshold_fine+log2_threshold_coarse)/2;
            std::size_t size_mid;
            int_quantized_compress(array, log2_threshold_fine, &size_mid, NULL, error_kind);
            if (size_mid<=target_size) {
              log2_threshold_coarse = log2_threshold_mid;
              size_coarse = size_mid;
            } else {
              log2_threshold_fine = log2_threshold_mid;
              size_fine = size_mid;
            }
          }
          // Choose largest size passing threshold
          if (size_fine<=target_size) {
            int_quantized_compress(array, log2_threshold_fine, actual_size, actual_error, error_kind);
          } else {
            int_quantized_compress(array, log2_threshold_coarse, actual_size, actual_error, error_kind);
          }
        }
      }
    }
    
    // Here, we have the result in the codec buffer, copy to output buffer
    if (*actual_size > buf_size) SL_FAIL("Codec buffer overflow");
    
    uint8_t* obuf = static_cast<uint8_t*>(buf);
    for (std::size_t i=0; i<*actual_size; ++i) {
      obuf[i] = entropy_codec_buffer_[i];
    }
  }
  
  /**
   *  Decompress buffer and store result into array.
   */
  template <class int_t>
  void hidwt_array_codec::int_decompress(dense_array<int_t,2,void>& array,
                                         const void* buf,
                                         std::size_t buf_size) {
    typedef hidwt_c22_matrix_codec<entropy_codec_t, int_t, 2> transform_codec_t;

    assert(buf);
    assert(buf_size >= 2);
    
    uint8_t* ibuf = static_cast<uint8_t*>(const_cast<void*>(buf));
    entropy_codec_.set_buffer(ibuf, buf_size);
    entropy_codec_.start_decoder();
    {
      uint32_t h              = entropy_codec_.decode_int();
      uint32_t w              = entropy_codec_.decode_int();
      uint32_t log2_threshold = entropy_codec_.decode_int();
      //std::cerr << "DECODED: h: " << h << " w= " << w << " T= " << log2_threshold << std::endl;
      array.resize(index<2>(h,w));
      transform_codec_t::decompress(entropy_codec_, array, log2_threshold);
    }
    entropy_codec_.stop_decoder();
    entropy_codec_.set_buffer(&(entropy_codec_buffer_[0]), entropy_codec_buffer_.size());

#if 0
    //////////
    std::cerr << "----- decompressed ------" << std::endl;
    for (std::size_t i=0; i< array.extent()[0]; ++i) {
      for (std::size_t j=0; j<array.extent()[1]; ++j) {
	std::cerr << std::setw(8) << array(i,j);
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
    /////////
#endif
  }


  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::int8_compress(const int8_array2_t& array,
                                        std::size_t target_size,
                                        double target_error,
                                        void* buf,
                                        std::size_t buf_size,
                                        std::size_t *actual_size, 
                                        double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::int8_compress_amax(const int8_array2_t& array,
                                             std::size_t target_size,
                                             double target_error,
                                             void* buf,
                                             std::size_t buf_size,
                                             std::size_t *actual_size, 
                                             double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::int8_decompress(int8_array2_t& array,
                                          const void* buf,
                                          std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::int16_compress(const int16_array2_t& array,
                                         std::size_t target_size,
                                         double target_error,
                                         void* buf,
                                         std::size_t buf_size,
                                         std::size_t *actual_size, 
                                         double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::int16_compress_amax(const int16_array2_t& array,
                                              std::size_t target_size,
                                              double target_error,
                                              void* buf,
                                              std::size_t buf_size,
                                              std::size_t *actual_size, 
                                              double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::int16_decompress(int16_array2_t& array,
                                           const void* buf,
                                           std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::int32_compress(const int32_array2_t& array,
                                         std::size_t target_size,
                                         double target_error,
                                         void* buf,
                                         std::size_t buf_size,
                                         std::size_t *actual_size, 
                                         double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::int32_compress_amax(const int32_array2_t& array,
                                              std::size_t target_size,
                                              double target_error,
                                              void* buf,
                                              std::size_t buf_size,
                                              std::size_t *actual_size, 
                                              double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::int32_decompress(int32_array2_t& array,
                                           const void* buf,
                                           std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::uint8_compress(const uint8_array2_t& array,
                                         std::size_t target_size,
                                         double target_error,
                                         void* buf,
                                         std::size_t buf_size,
                                         std::size_t *actual_size, 
                                         double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::uint8_compress_amax(const uint8_array2_t& array,
                                              std::size_t target_size,
                                              double target_error,
                                              void* buf,
                                              std::size_t buf_size,
                                              std::size_t *actual_size, 
                                              double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::uint8_decompress(uint8_array2_t& array,
                                           const void* buf,
                                           std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::uint16_compress(const uint16_array2_t& array,
                                          std::size_t target_size,
                                          double target_error,
                                          void* buf,
                                          std::size_t buf_size,
                                          std::size_t *actual_size, 
                                          double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }
        
  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::uint16_compress_amax(const uint16_array2_t& array,
                                               std::size_t target_size,
                                               double target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t *actual_size, 
                                               double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::uint16_decompress(uint16_array2_t& array,
                                            const void* buf,
                                            std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::uint32_compress(const uint32_array2_t& array,
                                          std::size_t target_size,
                                          double target_error,
                                          void* buf,
                                          std::size_t buf_size,
                                          std::size_t *actual_size, 
                                          double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   * Compress array to buffer
   */
  void hidwt_array_codec::uint32_compress_amax(const uint32_array2_t& array,
                                               std::size_t target_size,
                                               double target_error,
                                               void* buf,
                                               std::size_t buf_size,
                                               std::size_t *actual_size, 
                                               double *actual_error) {
    int_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }
        
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::uint32_decompress(uint32_array2_t& array,
                                            const void* buf,
                                            std::size_t buf_size) {
    int_decompress(array, buf, buf_size);
  }   

  
    
  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::float_compress(const float_array2_t& array,
                                         std::size_t  target_size,
                                         double       target_error,
                                         void* buf,
                                         std::size_t buf_size,
                                         std::size_t *actual_size,
                                         double      *actual_error,
					 error_kind_t error_kind) {
    assert(buf);
    assert(buf_size >= 8);

    enum { Int32_quantize_bits = 20 };
    const float Int32_quantize_scale =(float)(1<<Int32_quantize_bits);
    
    // Quantize to int32    
	const uint8_t h = static_cast<uint8_t>(array.extent()[0]);
	const uint8_t w = static_cast<uint8_t>(array.extent()[1]);
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
  void hidwt_array_codec::float_compress(const float_array2_t& array,
					 std::size_t  target_size,
					 double       target_error,
					 void* buf,
					 std::size_t buf_size,
					 std::size_t *actual_size,
					 double      *actual_error) {
    float_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_rms);
  }

  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::float_compress_amax(const float_array2_t& array,
					      std::size_t  target_size,
					      double       target_error,
					      void* buf,
					      std::size_t buf_size,
					      std::size_t *actual_size,
					      double      *actual_error) {
    float_compress(array, target_size, target_error, buf, buf_size, actual_size, actual_error, Error_kind_amax);
  }

  /**
   *  Decompress buffer and store result into array.
   */
  void hidwt_array_codec::float_decompress(float_array2_t& array,
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
	const uint8_t h = static_cast<uint8_t>(iarray.extent()[0]);
	const uint8_t w = static_cast<uint8_t>(iarray.extent()[1]);
    array.resize(subscript_t(h,w));
    for (uint32_t i=0; i<h; ++i) {
      for (uint32_t j=0; j<w; ++j) {
        array(i,j) = float(iarray(i,j))*qthr;
      }
    }
    
  }

}
