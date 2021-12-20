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
#ifndef SL_WAVELET_TRANSFORM_HPP
#define SL_WAVELET_TRANSFORM_HPP

#include <sl/utility.hpp>
#include <sl/numeric_traits.hpp>
#include <vector>
#include <cassert>

namespace sl {

  template <class T_VALUE>
  class wavelet_transform {
  public:
    typedef T_VALUE value_t;
  protected:
    bool        is_orthonormal_;
    bool        is_symmetric_;
    std::size_t analysis_low_size_;
    
  public:
    inline wavelet_transform()
        :
        is_orthonormal_(false),
        is_symmetric_(false),
        analysis_low_size_(0)
    {
    }

    virtual ~wavelet_transform() {
    }

    bool is_orthonormal() const { return is_orthonormal_; }

    bool is_symmetric() const { return is_symmetric_; }

    std::size_t analysis_low_size() const { return analysis_low_size_; }

    int analysis_low_first_index() const { return is_symmetric() ? (-(((int)analysis_low_size())-1)/2) : (0); }

    int analysis_low_last_index() const { return analysis_low_first_index() + (int)analysis_low_size() -1; }

  public: // Default implementation for orthogonal wavelets with no separate synthesis coeffs
    
    virtual std::size_t analysis_high_size() const { return analysis_low_size(); }

    virtual int analysis_high_first_index() const { return
                                                      2 -
                                                      (int)analysis_low_size()-
                                                      (int)analysis_low_first_index(); }

    int analysis_high_last_index() const { return analysis_high_first_index() + (int)analysis_high_size() -1; }
    
    virtual std::size_t synthesis_low_size() const { return analysis_low_size(); }

    virtual int synthesis_low_first_index() const { return analysis_low_first_index(); }

    int synthesis_low_last_index() const { return synthesis_low_first_index() + (int)synthesis_low_size() -1; }

    virtual std::size_t synthesis_high_size() const { return analysis_high_size(); }

    virtual int synthesis_high_first_index() const { return analysis_high_first_index(); }

    virtual int synthesis_high_last_index() const { return synthesis_high_first_index() + (int)synthesis_high_size() -1; }
    
  protected: // wavelet lifting support

    /**
     *  Split s into even and odd elements,
     *  where the even elements are in the first half
     *  of the vector and the odd elements are in the
     *  second half.
     */
    static void lifting_split(value_t* s, std::size_t N) {
      std::size_t start = 1;
      std::size_t end = N-1;

      while (start < end) {
        for (std::size_t i = start; i < end; i+=2) {
          std::swap(s[i],s[i+1]);
	}
        ++start; --end;
      }
    }

    /**
     * Merge the odd elements from the second half of the N element
     * region in the array with the even elements in the first
     * half of the N element region.  The result will be the
     * combination of the odd and even elements in a region
     * of length N.
     */
    static void lifting_merge(value_t* s, std::size_t N) {
      std::size_t half = N >> 1;
      std::size_t start = half-1;
      std::size_t end = half;
      
      while (start > 0) {
        for (std::size_t i = start; i < end; i += 2) {
          std::swap(s[i],s[i+1]);
        }
        --start; ++end;
      }
    }

    virtual void lifting_forward_predict(value_t* s,
                                         std::size_t l) const = 0;

    virtual void lifting_backward_predict(value_t* s,
                                          std::size_t l) const = 0;
    
    virtual void lifting_forward_update(value_t* s,
                                        std::size_t l) const = 0;
    
    virtual void lifting_backward_update(value_t* s,
                                        std::size_t l) const = 0;

    virtual void lifting_forward_normalize(value_t* /*s*/,
                                           std::size_t /*l*/) const {
      // Nothing by default
    }
    
    virtual void lifting_backward_normalize(value_t* /*s*/,
                                            std::size_t /*l*/) const {
      // Nothing by default
    }

  public:
    
    virtual void forward1d_step(value_t* s,
                                std::size_t l) const {
      lifting_split(s,l);
      lifting_forward_predict(s,l);
      lifting_forward_update(s,l);      
      lifting_forward_normalize(s,l);      
    }

    virtual void backward1d_step(value_t* s,
                                 std::size_t l) const {
      lifting_backward_normalize(s,l);      
      lifting_backward_update(s,l);      
      lifting_backward_predict(s,l);      
      lifting_merge(s,l);
    }

  public:

    template <class ARRAY1_T>
    void forward1d(ARRAY1_T& s, std::size_t N) const {
      std::size_t nn=1; while (nn<N) nn<<=1;
      if (nn==N) {
        // Assume array is a value_t* !!!
        for (std::size_t i=nn; i>=2; i>>=1) {
          forward1d_step(&s[0], i);
        }        
      } else {
        // Extent
        std::vector<value_t> ss(nn);
        for (std::size_t i=0; i<N; ++i) {
          ss[i] = s[i];
        }
        for (std::size_t i=N; i<nn; ++i) {
          ss[i] = s[N-2-(i-N)]; // fake values for non power of two signals
        } 
        for (std::size_t i=nn; i>=2; i>>=1) {
          forward1d_step(&ss[0], i);
        }
        for (std::size_t i=0; i<N; ++i) {
          s[i] = ss[i];
        }
      }
    }
      
    template <class ARRAY1_T>
    void forward1d(const ARRAY1_T& s_in, ARRAY1_T& s_out, std::size_t N) const {
      for (std::size_t i=0; i<N; ++i) {
        s_out[i] = s_in[i];
      }
      forward1d(s_out, N);
    }
    
    template <class ARRAY1_T>
    void backward1d(ARRAY1_T& s, std::size_t N) const {
      std::size_t nn=1; while (nn<N) nn<<=1;
      if (nn==N) {
        // Assume array is a value_t* !!!
        for (std::size_t i=2; i<=nn; i<<=1) {
          backward1d_step(&s[0], i);
        }
      } else {
        std::vector<value_t> ss(nn);
        for (std::size_t i=0; i<N; ++i) {
          ss[i] = s[i];
        }
        for (std::size_t i=N; i<nn; ++i) {
          ss[i] = 0; // missing detail for non power of two signals
        } 
        for (std::size_t i=2; i<=nn; i<<=1) {
          backward1d_step(&ss[0], i);
        }
        for (std::size_t i=0; i<N; ++i) {
          s[i] = ss[i];
        }
      }
    }

    template <class ARRAY1_T>
    void backward1d(const ARRAY1_T& s_in, ARRAY1_T& s_out, std::size_t N) const {
      for (std::size_t i=0; i<N; ++i) {
        s_out[i] = s_in[i];
      }
      backward1d(s_out, N);
    }
    
  public:

    template <class ARRAY2_T>
    void forward2d(ARRAY2_T& s,
                   std::size_t row_count,
                   std::size_t column_count) const {
      std::size_t nr=1; while (nr<row_count) nr<<=1;
      std::size_t nc=1; while (nc<column_count) nc<<=1;
      
      const std::size_t dmax = std::max(nr,nc);
      std::vector<value_t> ss(dmax);

      /* levels */
      for (std::size_t l=std::max(nr,nc); l>=2; l>>=1) {
        std::size_t lc = std::min(l,column_count);
        std::size_t lr = std::min(l,row_count);
        if (l<=nc) {
          /* rows */
          for(std::size_t j=0; j<std::min(l,row_count); ++j) {
            for (std::size_t i=0; i<lc; ++i) ss[i] = s(j,i);
            for (std::size_t i=lc; i<l; ++i) ss[i] = s(j,lc-2-(i-lc)); // fake values for non power of two signals
            forward1d_step(&ss[0],l);
            for (std::size_t i=0; i<lc; ++i) s(j,i) = ss[i];
          }
        }
        if (l<=nr) {
          /* cols */
          for(std::size_t i=0; i<std::min(l,column_count); ++i) {
            for (std::size_t j=0; j<lr; ++j) ss[j] = s(j,i);
            for (std::size_t j=lr; j<l; ++j) ss[j] = s(lr-2-(j-lr),i); // fake values for non power of two signals
            forward1d_step(&ss[0],l);
            for (std::size_t j=0; j<lr; ++j) s(j,i) = ss[j];
          }
        }
      }
    }
    
    template <class ARRAY2_T>
    void forward2d(const ARRAY2_T& s_in,
                   ARRAY2_T& s_out,
                   std::size_t row_count,
                   std::size_t column_count) const {
      for (std::size_t i=0; i<row_count; ++i) {
        for (std::size_t j=0; j<column_count; ++j) {
          s_out(i,j) = s_in(i,j);
        }
      }
      forward2d(s_out, row_count, column_count);
    }

    template <class ARRAY2_T>
    void backward2d(ARRAY2_T& s,
                   std::size_t row_count,
                   std::size_t column_count) const {
      std::size_t nr=1; while (nr<row_count) nr*=2;
      std::size_t nc=1; while (nc<column_count) nc*=2;
      
      const std::size_t dmax = std::max(nr,nc);
      std::vector<value_t> ss(dmax);

      /* levels */
      for (std::size_t l=2; l<=std::max(nr,nc); l<<=1) {
        std::size_t lc = std::min(l,column_count);
        std::size_t lr = std::min(l,row_count);
        if (l<=nc) {
          /* rows */
          for(std::size_t j=0; j<std::min(l,row_count); ++j) {
            for (std::size_t i=0; i<lc; ++i) ss[i] = s(j,i);
            for (std::size_t i=lc; i<l; ++i) ss[i] = 0; // missing detail for non power of two signals
            backward1d_step(&ss[0],l);
            for (std::size_t i=0; i<lc; ++i) s(j,i) = ss[i];
          }
        }
        if (l<=nr) {
          /* cols */
          for(std::size_t i=0; i<std::min(l,column_count); ++i) {
            for (std::size_t j=0; j<lr; ++j) ss[j] = s(j,i);
            for (std::size_t j=lr; j<l; ++j) ss[j] = 0; // missing detail for non power of two signals
            backward1d_step(&ss[0],l);
            for (std::size_t j=0; j<lr; ++j) s(j,i) = ss[j];
          }
        }
      }
    }
    
    template <class ARRAY2_T>
    void backward2d(const ARRAY2_T& s_in,
                    ARRAY2_T& s_out,
                    std::size_t row_count,
                    std::size_t column_count) const {
      for (std::size_t i=0; i<row_count; ++i) {
        for (std::size_t j=0; j<column_count; ++j) {
          s_out(i,j) = s_in(i,j);
        }
      }
      backward2d(s_out, row_count, column_count);
    }

  };

  /**
   *  Identity wavelet transform
   */
  template <class T_VALUE>
  class identity_wavelet_transform: public wavelet_transform<T_VALUE> {
  public:
    typedef T_VALUE value_t;
  public:
    inline identity_wavelet_transform() {
      this->is_orthonormal_ = true;
      this->is_symmetric_ = true;
      this->analysis_low_size_ = 1;
    }

    virtual ~identity_wavelet_transform() {
    }

  public:
    
    virtual void lifting_forward_predict(value_t* /*s*/,
                                         std::size_t /*l*/) const {
    }

    virtual void lifting_backward_predict(value_t* /*s*/,
                                         std::size_t /*l*/) const {
    }
    
    virtual void lifting_forward_update(value_t* /*s*/,
                                        std::size_t /*l*/) const {
    }
    
    virtual void lifting_backward_update(value_t* /*s*/,
                                        std::size_t /*l*/) const {
    }

    virtual void forward1d_step(value_t* /*s*/,
                                std::size_t /*l*/) const {
    }

    virtual void backward1d_step(value_t* /*s*/,
                                 std::size_t /*l*/) const {
    }
    
  };

  /**
   *  Haar wavelet transform
   */
  template <class T_VALUE>
  class haar_wavelet_transform: public wavelet_transform<T_VALUE> {
  public:
    typedef T_VALUE value_t;
  public:
    inline haar_wavelet_transform() {
      this->is_orthonormal_ = false;
      this->is_symmetric_ = false;
      this->analysis_low_size_   = 2;
    }

    virtual ~haar_wavelet_transform() {
    }
  public:

    static inline void haar_lifting_forward_predict(value_t* s,
                                                    std::size_t l) {
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i+half] -= s[i];
      }
    }

    static inline void haar_lifting_backward_predict(value_t* s,
                                                     std::size_t l) {
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i+half] += s[i];
      }
    }
    
    static inline void haar_lifting_forward_update(value_t* s,
                                                   std::size_t l) {
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i] += s[i+half]/2;
      }
    }
    
    static inline void haar_lifting_backward_update(value_t* s,
                                                    std::size_t l) {
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i] -= s[i+half]/2;
      }
    }

    virtual void lifting_forward_predict(value_t* s,
                                         std::size_t l) const {
      haar_lifting_forward_predict(s,l);
    }

    virtual void lifting_backward_predict(value_t* s,
                                         std::size_t l) const {
      haar_lifting_backward_predict(s,l);
    }
    
    virtual void lifting_forward_update(value_t* s,
                                        std::size_t l) const {
      haar_lifting_forward_update(s,l);
    }
    
    virtual void lifting_backward_update(value_t* s,
                                        std::size_t l) const {
      haar_lifting_backward_update(s,l);
    }

#if 0
    virtual void forward1d_step(value_t* s,
                                std::size_t l) const {
      lifting_split(s,l);
      haar_lifting_forward_predict(s,l);
      haar_lifting_forward_update(s,l); 
      lifting_forward_normalize(s,l);      
    }

    virtual void backward1d_step(value_t* s,
                                 std::size_t l) const {
      lifting_backward_normalize(s,l);      
      haar_lifting_backward_update(s,l);      
      haar_lifting_backward_predict(s,l);      
      lifting_merge(s,l);
    }
#endif
    
  };

  /**
   *  Normalized Haar wavelet transform
   */
  template <class T_VALUE>
  class normalized_haar_wavelet_transform: public haar_wavelet_transform<T_VALUE> {
  public:
    typedef T_VALUE value_t;
  public:
    inline normalized_haar_wavelet_transform() {
      this->is_orthonormal_ = true;
      this->is_symmetric_ = false;
      this->analysis_low_size_   = 2;
    }
    virtual ~normalized_haar_wavelet_transform() {}
  public:
    
    static void haar_lifting_forward_normalize(value_t* s,
                                                      std::size_t l) {
      const SL_FLOATTYPENAME(value_t) one_over_sqrt2 = 0.70710678118654752440;
      const SL_FLOATTYPENAME(value_t) sqrt2 = 1.41421356237309504880;
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i] *= sqrt2;
        s[i+half] *= one_over_sqrt2;
      }
    }

    static void haar_lifting_backward_normalize(value_t* s,
                                                       std::size_t l) {
      const SL_FLOATTYPENAME(value_t) one_over_sqrt2 = 0.70710678118654752440;
      const SL_FLOATTYPENAME(value_t) sqrt2 = 1.41421356237309504880;
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i] *= one_over_sqrt2;
        s[i+half] *= sqrt2;
      }
    }

    virtual void  lifting_forward_normalize(value_t* s,
                                            std::size_t l) const {
      return haar_lifting_forward_normalize(s,l);
    }
    
    virtual void  lifting_backward_normalize(value_t* s,
                                             std::size_t l) const {
      return haar_lifting_backward_normalize(s,l);
    }

#if 0
    virtual void forward1d_step(value_t* s,
                                std::size_t l) const {
      lifting_split(s,l);
      haar_lifting_forward_predict(s,l);
      haar_lifting_forward_update(s,l);      
      haar_lifting_forward_normalize(s,l);      
    }

    virtual void backward1d_step(value_t* s,
                                 std::size_t l) const {
      haar_lifting_backward_normalize(s,l);      
      haar_lifting_backward_update(s,l);      
      haar_lifting_backward_predict(s,l);      
      lifting_merge(s,l);
    }
#endif

  };

  /**
   *  TS wavelet transform
   */
  template <class T_VALUE>
  class ts_wavelet_transform: public normalized_haar_wavelet_transform<T_VALUE> {
  public:
    typedef T_VALUE value_t;
    typedef ts_wavelet_transform<value_t>   this_t;
    typedef haar_wavelet_transform<value_t> super_t;
    
  public:
    inline ts_wavelet_transform() {
      this->is_orthonormal_ = false;
      this->is_symmetric_ = false;
      this->analysis_low_size_   = 4;
    }
    virtual ~ts_wavelet_transform() {}
  public:

    virtual void lifting_forward_predict_lerp(value_t* s,
                                              std::size_t l) const {
      std::size_t half = l>>1;
      for (std::size_t i=0; i<half; ++i) {
        value_t s_i_plus1;
        value_t s_i_minus1;
        if (l==2) {
          s_i_minus1 = s[0];
          s_i_plus1= s[0];
        } else if (i==0) {
          s_i_minus1 = 2 * s[0] - s[1];
          s_i_plus1 = s[1];
        } else if (i==half-1) {
          s_i_minus1 = s[i-1];
          s_i_plus1 = 2 * s[i] - s[i-1];
        } else {
          s_i_minus1 = s[i-1];
          s_i_plus1 = s[i+1];
        }
        s[i+half] += (s_i_minus1 - s_i_plus1)/4;
      }
    }
    
    virtual void lifting_backward_predict_lerp(value_t* s,
                                              std::size_t l) const {
      std::size_t half = l>>1;
      for (std::size_t i=0; i<half; ++i) {
        value_t s_i_plus1;
        value_t s_i_minus1;
        if (l==2) {
          s_i_minus1 = s[0];
          s_i_plus1= s[0];
        } else if (i==0) {
          s_i_minus1 = 2 * s[0] - s[1];
          s_i_plus1 = s[1];
        } else if (i==half-1) {
          s_i_minus1 = s[i-1];
          s_i_plus1 = 2 * s[i] - s[i-1];
        } else {
          s_i_minus1 = s[i-1];
          s_i_plus1 = s[i+1];
        }
        s[i+half] -= (s_i_minus1 - s_i_plus1)/4;
      }
    }

    virtual void lifting_forward_update(value_t* s,
                                        std::size_t l) const {
      super_t::lifting_forward_update(s,l);
      lifting_forward_predict_lerp(s,l);
    }

    virtual void lifting_backward_update(value_t* s,
                                         std::size_t l) const {
      lifting_backward_predict_lerp(s,l);
      super_t::lifting_backward_update(s,l);
    }
    
  };
  
  /**
   *  Daubechies-4 wavelet transform
   */
  template <class T_VALUE>
  class normalized_daubechies4_wavelet_transform: public normalized_haar_wavelet_transform<T_VALUE> {
  public:
    typedef T_VALUE value_t;
  public:
    inline normalized_daubechies4_wavelet_transform() {
      this->is_orthonormal_ = true;
      this->is_symmetric_ = false;
      this->analysis_low_size_   = 4;
    }

    virtual ~normalized_daubechies4_wavelet_transform() {}
  public:

    virtual void lifting_forward_update_one(value_t* s,
                                            std::size_t l) const {
      assert(l>=4);

      const SL_FLOATTYPENAME(value_t) sqrt3 = 1.73205080756887729352;
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i] += sqrt3 * s[i+half];
      }
    }

    virtual void lifting_backward_update_one(value_t* s,
                                             std::size_t l) const {
      assert(l>=4);

      const SL_FLOATTYPENAME(value_t) sqrt3 = 1.73205080756887729352;
      const std::size_t half = l >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        s[i] -= sqrt3 * s[i+half];
      }
    }
    
    virtual void lifting_forward_predict(value_t* s,
                                         std::size_t l) const {
      if (l<4) {
        normalized_haar_wavelet_transform<value_t>::lifting_forward_predict(s,l);
      } else {
        lifting_forward_update_one(s,l);
        
        const SL_FLOATTYPENAME(value_t) sqrt3_div4        =  0.43301270189221932338;
        const SL_FLOATTYPENAME(value_t) sqrt3_minus2_div4 = -0.06698729810778067662;
        const std::size_t half = l >> 1;
        s[half] -= sqrt3_div4 * s[0] + sqrt3_minus2_div4 * s[half-1];
        for (std::size_t i = 1; i < half; ++i) {
          s[i+half] -= sqrt3_div4 * s[i] + sqrt3_minus2_div4 * s[i-1];
        }
      }
    }

    virtual void lifting_backward_predict(value_t* s,
                                         std::size_t l) const {
      if (l<4) {
        normalized_haar_wavelet_transform<value_t>::lifting_backward_predict(s,l);
      } else {
        const SL_FLOATTYPENAME(value_t) sqrt3_div4        =  0.43301270189221932338;
        const SL_FLOATTYPENAME(value_t) sqrt3_minus2_div4 = -0.06698729810778067662;
        const std::size_t half = l >> 1;
        s[half] += sqrt3_div4 * s[0] + sqrt3_minus2_div4 * s[half-1];
        for (std::size_t i = 1; i < half; ++i) {
          s[i+half] += sqrt3_div4 * s[i] + sqrt3_minus2_div4 * s[i-1];
        }
        
        lifting_backward_update_one(s,l);
      }
    }
    
    virtual void lifting_forward_update(value_t* s,
                                        std::size_t l) const {
      if (l<4) {
        normalized_haar_wavelet_transform<value_t>::lifting_forward_update(s,l);
      } else {
        const std::size_t half = l >> 1;
        for (std::size_t i = 0; i < half-1; ++i) {
          s[i] -= s[i+half+1];
        }
        s[half-1] -=s[half];
      }
    }
    
    virtual void lifting_backward_update(value_t* s,
                                        std::size_t l) const {
      if (l<4) {
        normalized_haar_wavelet_transform<value_t>::lifting_backward_update(s,l);
      } else {
        const std::size_t half = l >> 1;
        for (std::size_t i = 0; i < half-1; ++i) {
          s[i] += s[i+half+1];
        }
        s[half-1] +=s[half];
      }
    }

    virtual void lifting_forward_normalize(value_t* s,
                                           std::size_t l) const {
      if (l<4) {
        normalized_haar_wavelet_transform<value_t>::lifting_forward_normalize(s,l);
      } else {
        const SL_FLOATTYPENAME(value_t) sqrt3_minus_1_div_sqrt2= 0.51763809020504152469;
        const SL_FLOATTYPENAME(value_t) sqrt3_plus_1_div_sqrt2 = 1.93185165257813657349;
        
        const std::size_t half = l >> 1;
        for (std::size_t i = 0; i < half; ++i) {
          s[i     ] *= sqrt3_minus_1_div_sqrt2;
          s[i+half] *= sqrt3_plus_1_div_sqrt2;
        }
      }
    }
    
    virtual void lifting_backward_normalize(value_t* s,
                                            std::size_t l) const {
      if (l<4) {
        normalized_haar_wavelet_transform<value_t>::lifting_backward_normalize(s,l);
      } else {
        const SL_FLOATTYPENAME(value_t) sqrt3_minus_1_div_sqrt2= 0.51763809020504152469;
        const SL_FLOATTYPENAME(value_t) sqrt3_plus_1_div_sqrt2 = 1.93185165257813657349;
        
        const std::size_t half = l >> 1;
        for (std::size_t i = 0; i < half; ++i) {
          s[i     ] *= sqrt3_plus_1_div_sqrt2;
          s[i+half] *= sqrt3_minus_1_div_sqrt2;
        }
      }
    }

  };

  /**
   *  Daubechies-4 wavelet transform.
   *
   *  Adapted from FastSymmlet8.c, (c) 1998-2002 Daniel Lemire, http://www.ondelette.com/
   */
  template <class T_VALUE>
  class normalized_symmlet8_wavelet_transform: public normalized_daubechies4_wavelet_transform<T_VALUE> {
  public:
    typedef T_VALUE value_t;
    typedef normalized_daubechies4_wavelet_transform<T_VALUE> super_t;
  public:
    inline normalized_symmlet8_wavelet_transform() {
      this->is_orthonormal_ = true;
      this->is_symmetric_ = false;
      this->analysis_low_size_   = 4;
    }

    virtual ~normalized_symmlet8_wavelet_transform() {}

  public:// FIXME implement lifting version

    virtual void forward1d_step(value_t* s,
                                std::size_t l) const {
      if (l<8) {
        super_t::forward1d_step(s,l);
      } else {
        static const value_t scale[8] = {0.0322231006040782, -0.0126039672622638, -0.0992195435769564,  0.297857795605605,
                                         0.803738751805386,   0.497618667632563,  -0.0296355276459604, -0.0757657147893567};
        static const float wavelet[8] = {0.0757657147893567, -0.0296355276459604, -0.497618667632563,   0.803738751805386,
                                         -0.297857795605605,  -0.0992195435769564,  0.0126039672622638,  0.0322231006040782};
      
        value_t * ans=new float[l];
	for(std::size_t k=0;k<l;k++) {
          ans[k]=0;
	}
        std::size_t half=l/2;
        for(std::size_t k=0;k<half-3;++k) {
          ans[k+half]=s[(2*k+0)]*wavelet[0]+s[(2*k+1)]*wavelet[1]+s[(2*k+2)]*wavelet[2]+s[(2*k+3)]*wavelet[3]+s[(2*k+4)]*wavelet[4]+s[(2*k+5)]*wavelet[5]+s[(2*k+6)]*wavelet[6]+s[(2*k+7)]*wavelet[7];
          ans[k]=s[(2*k+0)]*scale[0]+s[(2*k+1)]*scale[1]+s[(2*k+2)]*scale[2]+s[(2*k+3)]*scale[3]+s[(2*k+4)]*scale[4]+s[(2*k+5)]*scale[5]+s[(2*k+6)]*scale[6]+s[(2*k+7)]*scale[7];
	}
        ans[l-3]=s[l-6]*wavelet[0]+s[l-5]*wavelet[1]+s[l-4]*wavelet[2]+s[l-3]*wavelet[3]+s[l-2]*wavelet[4]+s[l-1]*wavelet[5]+s[0]*wavelet[6]+s[1]*wavelet[7];
        ans[half-3]=s[l-6]*scale[0]+s[l-5]*scale[1]+s[l-4]*scale[2]+s[l-3]*scale[3]+s[l-2]*scale[4]+s[l-1]*scale[5]+s[0]*scale[6]+s[1]*scale[7];
        ans[l-2]=s[l-4]*wavelet[0]+s[l-3]*wavelet[1]+s[l-2]*wavelet[2]+s[l-1]*wavelet[3]+s[0]*wavelet[4]+s[1]*wavelet[5]+s[2]*wavelet[6]+s[3]*wavelet[7];
        ans[half-2]=s[l-4]*scale[0]+s[l-3]*scale[1]+s[l-2]*scale[2]+s[l-1]*scale[3]+s[0]*scale[4]+s[1]*scale[5]+s[2]*scale[6]+s[3]*scale[7];
        ans[l-1]=s[l-2]*wavelet[0]+s[l-1]*wavelet[1]+s[0]*wavelet[2]+s[1]*wavelet[3]+s[2]*wavelet[4]+s[3]*wavelet[5]+s[4]*wavelet[6]+s[5]*wavelet[7];
        ans[half-1]=s[l-2]*scale[0]+s[l-1]*scale[1]+s[0]*scale[2]+s[1]*scale[3]+s[2]*scale[4]+s[3]*scale[5]+s[4]*scale[6]+s[5]*scale[7];
	for(std::size_t k=0;k<l;++k) {
          s[k]=ans[k];
	}
	delete[] ans;
      }
    }

    virtual void backward1d_step(value_t* s,
                                 std::size_t l) const {
      if (l<8) {
        super_t::backward1d_step(s,l);
      } else {
        static const value_t scale[8] = {0.0322231006040782, -0.0126039672622638, -0.0992195435769564,  0.297857795605605,
                                         0.803738751805386,   0.497618667632563,  -0.0296355276459604, -0.0757657147893567};
        static const float wavelet[8] = {0.0757657147893567, -0.0296355276459604, -0.497618667632563,   0.803738751805386,
                                         -0.297857795605605,  -0.0992195435769564,  0.0126039672622638,  0.0322231006040782};
        
        const std::size_t half = l >> 1;
        float * ans=new float[l];
	for(std::size_t k=0;k<l;k++) {
		ans[k]=0;
	}
        for(std::size_t k=0;2*k+7<l;k++) {
          ans[(2*k+7)]+=scale[7]*s[k]+wavelet[7]*s[k+half];
          ans[(2*k+6)]+=scale[6]*s[k]+wavelet[6]*s[k+half];
          ans[(2*k+5)]+=scale[5]*s[k]+wavelet[5]*s[k+half];
          ans[(2*k+4)]+=scale[4]*s[k]+wavelet[4]*s[k+half];
          ans[(2*k+3)]+=scale[3]*s[k]+wavelet[3]*s[k+half];
          ans[(2*k+2)]+=scale[2]*s[k]+wavelet[2]*s[k+half];
          ans[(2*k+1)]+=scale[1]*s[k]+wavelet[1]*s[k+half];
          ans[(2*k+0)]+=scale[0]*s[k]+wavelet[0]*s[k+half];
        }
        ans[l-6]+=scale[0]*s[half-3]+wavelet[0]*s[l-3];
        ans[l-5]+=scale[1]*s[half-3]+wavelet[1]*s[l-3];
        ans[l-4]+=scale[2]*s[half-3]+wavelet[2]*s[l-3];
        ans[l-3]+=scale[3]*s[half-3]+wavelet[3]*s[l-3];
        ans[l-2]+=scale[4]*s[half-3]+wavelet[4]*s[l-3];
        ans[l-1]+=scale[5]*s[half-3]+wavelet[5]*s[l-3];
        ans[0]+=scale[6]*s[half-3]+wavelet[6]*s[l-3];
        ans[1]+=scale[7]*s[half-3]+wavelet[7]*s[l-3];
        ans[l-4]+=scale[0]*s[half-2]+wavelet[0]*s[l-2];
        ans[l-3]+=scale[1]*s[half-2]+wavelet[1]*s[l-2];
        ans[l-2]+=scale[2]*s[half-2]+wavelet[2]*s[l-2];
        ans[l-1]+=scale[3]*s[half-2]+wavelet[3]*s[l-2];
        ans[0]+=scale[4]*s[half-2]+wavelet[4]*s[l-2];
        ans[1]+=scale[5]*s[half-2]+wavelet[5]*s[l-2];
        ans[2]+=scale[6]*s[half-2]+wavelet[6]*s[l-2];
        ans[3]+=scale[7]*s[half-2]+wavelet[7]*s[l-2];
        ans[l-2]+=scale[0]*s[half-1]+wavelet[0]*s[l-1];
        ans[l-1]+=scale[1]*s[half-1]+wavelet[1]*s[l-1];
        ans[0]+=scale[2]*s[half-1]+wavelet[2]*s[l-1];
        ans[1]+=scale[3]*s[half-1]+wavelet[3]*s[l-1];
        ans[2]+=scale[4]*s[half-1]+wavelet[4]*s[l-1];
        ans[3]+=scale[5]*s[half-1]+wavelet[5]*s[l-1];
        ans[4]+=scale[6]*s[half-1]+wavelet[6]*s[l-1];
        ans[5]+=scale[7]*s[half-1]+wavelet[7]*s[l-1];
	for(std::size_t k=0;k<l;++k) {
          s[k]=ans[k];
	}
	delete[] ans;
      }
    }

  };
  
} // namespace sl

#endif
