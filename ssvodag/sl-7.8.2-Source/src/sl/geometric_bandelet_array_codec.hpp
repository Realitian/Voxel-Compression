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
#ifndef SL_GEOMETRIC_BANDELET_ARRAY_CODEC_HPP
#define SL_GEOMETRIC_BANDELET_ARRAY_CODEC_HPP

#include <sl/wavelet_array_codec.hpp>
#include <sl/fixed_ac_int_codec.hpp>
#include <sl/fixed_huffman_rle_codec.hpp>
#include <sl/fixed_rc_int_codec.hpp>


#define BANDELET_TRACE 0

namespace sl {

  /**
   *  Lossy compression of arbitrary scalar arrays using a geometric
   *  bandelet coder.
   *
   *  Based on "Surface Compression With Geometric Bandelets"
   *  Gabriel Peyré and Stéphane Mallat,
   *  Proceedings of SIGGRAPH'05.
   */
  class geometric_bandelet_array_codec: public wavelet_array_codec {
  public:
    typedef geometric_bandelet_array_codec this_t;
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
    
    typedef fixed_rc_int_codec             entropy_codec_t;
    typedef entropy_codec_t::bit_context_t bit_context_t;
    typedef entropy_codec_t::int_context_t int_context_t;
    
  protected: // Compressed header

    class compressed_data_header {
    protected:
      std::size_t      h_;
      std::size_t      w_;
      transform_kind_t transform_kind_;
      float            signal_mean_;
      float            xform_threshold_;
    public:
      compressed_data_header(std::size_t h,
                             std::size_t w,
                             const transform_kind_t& transform_kind,
                             float signal_mean,
                             float xform_threshold)
          :
          h_(h),
          w_(w),
          transform_kind_(transform_kind),
          signal_mean_(signal_mean),
          xform_threshold_(xform_threshold) {
      }

      compressed_data_header(const void* header_buffer) {
        const uint8_t* buf = (const uint8_t*)header_buffer;
        uint8_t hdr = be_decode<uint8_t>(buf); buf += 1;
        const bool square_matrix  = (hdr & 0x80)!=0;
        const bool h_16bit        = (hdr & 0x40)!=0;
        const bool w_16bit        = (hdr & 0x20)!=0;
        const bool zero_mean      = (hdr & 0x10)!=0;
        transform_kind_ = transform_kind_t(hdr & (~(0x80|0x40|0x20|0x10|0x08)));
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
        xform_threshold_     = be_decode<float>(buf); buf+=4;
        assert(std::size_t(buf-(const uint8_t*)header_buffer) == compressed_size());
      }

      void encode(void* header_buffer) {
        uint8_t* buf = (uint8_t*)header_buffer;
        const bool square_matrix = h() == w();
        const bool h_16bit = h() > 0xff;
        const bool w_16bit = w() > 0xff;
        const bool zero_mean = (signal_mean_ == 0.0f);
        uint8_t hdr = uint8_t(transform_kind());
        assert(hdr <= 0x08);
        
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
        buf += be_encode(xform_threshold(), buf);
        assert(std::size_t(buf-(uint8_t*)header_buffer) == compressed_size());
      }

      inline std::size_t h() const { return h_; }
      inline std::size_t w() const { return w_; }
      inline transform_kind_t transform_kind() const { return transform_kind_; }
      inline float signal_mean() const { return signal_mean_; }
      inline float xform_threshold() const { return xform_threshold_; }

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
      
  public: // Construction and destruction

    /// Construct codec
    inline geometric_bandelet_array_codec() {
      bandelet_init_idxmap();
    }

    /// Destruct codec
    inline ~geometric_bandelet_array_codec() {
    }
    
    virtual std::string description() const { return std::string("Bandelet") + "/" + this->transform_name_[transform_kind_]; }

  protected: // Bandelet compression helpers

    entropy_codec_t entropy_codec_;

    typedef enum {
      Bandelet_size_min_quadtree_root = 16,  // Start quadtree decomposition here or above
      Bandelet_size_min_geometry = 4,        // No geometry search for matrices smaller than this
      Bandelet_size_max_geometry = 32        // No geometry search for matrices larger than this 
    } bandelet_limits_t;
    
    inline std::size_t bandelet_direction_count(std::size_t sz) const {
      return (sz<Bandelet_size_min_geometry || sz>Bandelet_size_max_geometry) ? 0 : 2*int32_t(sz)-2;
    }

    inline float bandelet_direction_theta(int32_t i, std::size_t sz) const {
      const int32_t N_dir = bandelet_direction_count(sz);
      return (N_dir==0 || i<0) ? 0.0f : scalar_math<float>::Pi()/float(N_dir)*float(i);
    }

  protected: // Index sorting
    
    typedef std::pair<uint8_t, uint8_t>            small_subscript_t;
    typedef std::pair<std::size_t, std::size_t>    size_direction_pair_t;
    
    std::vector<small_subscript_t>                 bandelet_idxmap_data_; // Make it shared
    std::map<size_direction_pair_t, std::size_t>   bandelet_idxmap_offset_;

    struct bandelet_cmp {
      float cos_theta_;
      float sin_theta_;
      bandelet_cmp(float theta) {
        cos_theta_ = std::cos(theta);
        sin_theta_ = std::sin(theta);
      }
      bool operator()(const small_subscript_t& i0,
                      const small_subscript_t& i1) const {
        const float t0 = -sin_theta_*i0.first + cos_theta_*i0.second;
        const float t1 = -sin_theta_*i1.first + cos_theta_*i1.second;
        return t0<t1;
      }
    };
    
    void bandelet_init_idxmap() {
      if (bandelet_idxmap_offset_.empty()) {
        for (std::size_t sz=Bandelet_size_min_geometry; sz<=Bandelet_size_max_geometry; ++sz) {
          const std::size_t sz2 = sz*sz;
          for (std::size_t direction=0; direction<bandelet_direction_count(sz); ++direction) {
            std::vector<small_subscript_t> idxmap(sz2);

            for (std::size_t i=0; i<sz; ++i) {
              for (std::size_t j=0; j<sz; ++j) {
                idxmap[i*sz+j] = small_subscript_t(i,j);
              }
            }
            const float theta = bandelet_direction_theta(direction,sz);
            std::stable_sort(idxmap.begin(), idxmap.begin()+sz2, bandelet_cmp(theta));

            bandelet_idxmap_offset_[size_direction_pair_t(sz,direction)] = bandelet_idxmap_data_.size();
            for (std::size_t k=0; k<sz2; ++k) {
              bandelet_idxmap_data_.push_back(idxmap[k]);
            }
          }
        }
      }
    }
    
    void bandelet_sort_forward(int32_t direction,
                               const float_array2_t& wm,
                               std::size_t i0,
                               std::size_t j0,
                               std::size_t sz,
                               float* ss) {
      if (direction<0 ||
          sz < Bandelet_size_min_geometry ||
          sz > Bandelet_size_max_geometry) {
        // No geometry, Morton order
        typedef morton_bitops<std::size_t,2> morton2d_t;
        for (std::size_t i=0; i<sz; ++i) {
          for (std::size_t j=0; j<sz; ++j) {
            const std::size_t idx_k =
              morton2d_t::encoded(i,0) |
              morton2d_t::encoded(j,1);
            ss[idx_k] = wm(i0+i, j0+j);
          }
        }
      } else {
        // Geometry, orthogonal projection
        const small_subscript_t *idxmap =
          &(bandelet_idxmap_data_[bandelet_idxmap_offset_[size_direction_pair_t(sz,direction)]]);
        const std::size_t sz2 = sz*sz;
        for (std::size_t k=0; k<sz2; ++k) {
          const small_subscript_t& idx_k = idxmap[k];
          ss[k] = wm(i0+idx_k.first,
                     j0+idx_k.second);
        }
      }
    }
    
    void bandelet_sort_backward(int32_t direction,
                                float_array2_t& wm,
                                std::size_t i0,
                                std::size_t j0,
                                std::size_t sz,
                                const float* ss) {
      if (direction<0 ||
          sz < Bandelet_size_min_geometry ||
          sz > Bandelet_size_max_geometry) {
        // No geometry, Morton order
        typedef morton_bitops<std::size_t,2> morton2d_t;
        for (std::size_t i=0; i<sz; ++i) {
          for (std::size_t j=0; j<sz; ++j) {
            const std::size_t idx_k =
              morton2d_t::encoded(i,0) |
              morton2d_t::encoded(j,1);
            wm(i0+i, j0+j) = ss[idx_k];
          }
        }
      } else {
        // Geometry, orthogonal projection
        const small_subscript_t *idxmap =
          &(bandelet_idxmap_data_[bandelet_idxmap_offset_[size_direction_pair_t(sz,direction)]]);
        const std::size_t sz2 = sz*sz;
        for (std::size_t k=0; k<sz2; ++k) {
          const small_subscript_t& idx_k = idxmap[k];
          wm(i0+idx_k.first,j0+idx_k.second) = ss[k];
        }
      }
    }
    
    typedef enum { Bandelet_tag_split = -3, Bandelet_tag_zero_block = -2, Bandelet_tag_untransformed = -1 } bandelet_tag_t;

    class bandelet_quad_descriptor {
    protected:
      int32_t   tag_;
      uint32_t  bit_count_;
      float     sqr_err_;
    public:
      inline bandelet_quad_descriptor(int32_t  tag = Bandelet_tag_untransformed,
                                      uint32_t bit_count = 0,
                                      float    sqr_err = 0) :
          tag_(tag), bit_count_(bit_count), sqr_err_(sqr_err) {
        
      }
      inline int32_t   tag() const { return tag_; }
      inline uint32_t  bit_count() const { return bit_count_; }
      inline float     sqr_err() const { return sqr_err_; }
      inline int32_t&  tag() { return tag_; }
      inline uint32_t& bit_count(){ return bit_count_; }
      inline float&    sqr_err(){ return sqr_err_; }
      inline float     lagrangian(float threshold) const {
        const float One_Over_Lambda = 4.0/3.0f; 
        return One_Over_Lambda * sqr_err() / (threshold * threshold) + bit_count();
      }
    };
    
    typedef bandelet_quad_descriptor                    bandelet_quad_descriptor_t;
    typedef std::vector<bandelet_quad_descriptor_t>     bandelet_quadtree_t;
    
    bandelet_quadtree_t                                 bandelet_quadtree_;
    std::vector<float>                                  bandelet_transformed_data_;
    
    void bandelet_begin_encoding(const float_array2_t& wm) {
      bandelet_quadtree_.resize(2*wm.extent()[0]*wm.extent()[1]);
      bandelet_transformed_data_.resize(2*wm.extent()[0]*wm.extent()[1]);
    }

    void bandelet_end_encoding() {
    }

    void bandelet_begin_decoding(float_array2_t& wm) {
      //bandelet_quadtree_.resize(2*wm.extent()[0]*wm.extent()[1]);
      bandelet_transformed_data_.resize(2*wm.extent()[0]*wm.extent()[1]);
    }

    void bandelet_end_decoding() {
    }
    
  protected:
    
    // Forward bandelet transform of a square subimage of wm
    void bandelet_forward(int32_t direction,
                          const float_array2_t& wm,
                          std::size_t i0,
                          std::size_t j0,
                          std::size_t sz,
                          float* ss) {
      bandelet_sort_forward(direction, wm, i0, j0, sz, ss);
      if (direction>=0) {
        transform_table_[TRANSFORM_KIND_HAAR]->forward1d(ss, sz*sz);
      }
    }

    // Backward bandelet transform of a square subimage of wm
    void bandelet_backward(int32_t direction,
                           float_array2_t& wm,
                           std::size_t i0,
                           std::size_t j0,
                           std::size_t sz,
                           float* ss) {
      if (direction>=0) {
        transform_table_[TRANSFORM_KIND_HAAR]->backward1d(ss, sz*sz);
      }
      bandelet_sort_backward(direction, wm, i0, j0, sz, ss);
    }

  protected:

    void bandelet_encode(entropy_codec_t& entropy_codec,
                         const float_array2_t& wm,
                         bool is_wavelet_xform,
                         float threshold,
                         std::size_t* actual_size,
                         double* actual_rms) {
#if BANDELET_TRACE
      std::cerr << "BEGIN ENCODE: " << wm.extent()[0] << "x" << wm.extent()[1] << " thr=" << threshold << std::endl;
#endif
      bandelet_begin_encoding(wm);
      if (is_wavelet_xform) {
        bandelet_encode_subtree(entropy_codec,
                                wm,
                                0, 0, wm.extent()[0],
                                threshold,
                                actual_size,
                                actual_rms);
      } else {
        bandelet_encode_subimage(entropy_codec,
                                 wm,
                                 0, 0, wm.extent()[0],
                                 threshold,
                                 actual_size,
                                 actual_rms);
      }
      bandelet_end_encoding();
#if BANDELET_TRACE
      std::cerr << "END ENCODE: " << wm.extent()[0] << "x" << wm.extent()[1] << " thr=" << threshold << std::endl;
#endif
    }
    
    void bandelet_decode(entropy_codec_t& entropy_codec,
                         float_array2_t& wm,
                         bool is_wavelet_xform,
                         float threshold) {
#if BANDELET_TRACE
      std::cerr << "BEGIN DECODE: " << wm.extent()[0] << "x" << wm.extent()[1] << " thr=" << threshold << std::endl;
#endif
      bandelet_begin_decoding(wm);
      if (is_wavelet_xform) {
        bandelet_decode_subtree(entropy_codec,
                                wm,
                                0, 0, wm.extent()[0],
                                threshold);
      } else {
        bandelet_decode_subimage(entropy_codec,
                                 wm,
                                 0, 0, wm.extent()[0],
                                 threshold);
      }        
      bandelet_end_decoding();
#if BANDELET_TRACE
      std::cerr << "END DECODE: " << wm.extent()[0] << "x" << wm.extent()[1] << " thr=" << threshold << std::endl;
#endif
    }

    void bandelet_encode_subtree(entropy_codec_t& entropy_codec,
                                 const float_array2_t& wm,
                                 std::size_t i0, std::size_t j0,
                                 std::size_t sz,
                                 float threshold,
                                 std::size_t* actual_size,
                                 double* actual_rms) {
      if (sz<=Bandelet_size_min_quadtree_root) {
        // Root
        bandelet_encode_subimage(entropy_codec, wm, i0, j0, sz, threshold, actual_size, actual_rms);
      } else {
        // Inner node of wavelet tree
        const std::size_t child_sz = sz/2;
        std::size_t sz0, sz1, sz2, sz3;
        double rms0, rms1, rms2, rms3;
        bandelet_encode_subtree (entropy_codec, wm, i0,          j0,          child_sz, threshold, &sz0, &rms0); // LL
        bandelet_encode_subimage(entropy_codec, wm, i0+child_sz, j0,          child_sz, threshold, &sz1, &rms1); // HL
        bandelet_encode_subimage(entropy_codec, wm, i0,          j0+child_sz, child_sz, threshold, &sz2, &rms2); // LH
        bandelet_encode_subimage(entropy_codec, wm, i0+child_sz, j0+child_sz, child_sz, threshold, &sz3, &rms3); // HH
        *actual_size = sz0+sz1+sz2+sz3;
        *actual_rms = std::sqrt(0.25f*(rms0*rms0+rms1*rms1+rms2*rms2+rms3*rms3));
      }
    }

    void bandelet_decode_subtree(entropy_codec_t& entropy_codec,
                                 float_array2_t& wm,
                                 std::size_t i0, std::size_t j0,
                                 std::size_t sz,
                                 float threshold) {
      if (sz<=Bandelet_size_min_quadtree_root) {
        // Root
        bandelet_decode_subimage(entropy_codec, wm, i0, j0, sz, threshold);
      } else {
        // Inner node of wavelet tree
        const std::size_t child_sz = sz/2;
        bandelet_decode_subtree (entropy_codec, wm, i0,          j0,          child_sz, threshold); // LL
        bandelet_decode_subimage(entropy_codec, wm, i0+child_sz, j0,          child_sz, threshold); // HL
        bandelet_decode_subimage(entropy_codec, wm, i0,          j0+child_sz, child_sz, threshold); // LH
        bandelet_decode_subimage(entropy_codec, wm, i0+child_sz, j0+child_sz, child_sz, threshold); // HH
      }
    }
    
    void bandelet_encode_subimage(entropy_codec_t& entropy_codec,
                                  const float_array2_t& wm,
                                  std::size_t i0, std::size_t j0,
                                  std::size_t sz,
                                  float threshold,
                                  std::size_t* actual_size,
                                  double* actual_rms) {
      bandelet_compute_optimal_quadtree(0, wm, i0, j0, sz,
                                        threshold);
      *actual_rms = std::sqrt(bandelet_quadtree_[0].sqr_err()/double(sz*sz));
      
      std::size_t sz0 = entropy_codec.current_bit_count();

      bit_context_t quadtree_context;
      int_context_t tag_context;
      int_context_t coeff_context;
      bandelet_encode_optimal_quadtree(entropy_codec,
                                       quadtree_context,
                                       tag_context,
                                       coeff_context,
                                       0, wm, i0, j0, sz,
                                       threshold);
      *actual_size = (entropy_codec.current_bit_count() - sz0)/8;

      //std::cerr << "SZ[" << sz << "x" << sz << "]: " << "tag = " << bandelet_quadtree_[0].tag() << " " << bandelet_quadtree_[0].bit_count() << " vs. " << (*actual_size)*8 << std::endl;
    }

    void bandelet_decode_subimage(entropy_codec_t& entropy_codec,
                                  float_array2_t& wm,
                                  std::size_t i0, std::size_t j0,
                                  std::size_t sz,
                                  float threshold) {
      bit_context_t quadtree_context;
      int_context_t tag_context;
      int_context_t coeff_context;
      bandelet_decode_optimal_quadtree(entropy_codec,
                                       quadtree_context,
                                       tag_context,
                                       coeff_context,
                                       wm, i0, j0, sz,
                                       threshold);
    }
    
    void bandelet_compute_optimal_quadtree(std::size_t root_idx,
                                           const float_array2_t& wm,
                                           std::size_t i0,
                                           std::size_t j0,
                                           std::size_t sz,
                                           float threshold) {
      if (sz <= Bandelet_size_min_geometry) {
        // Leaf
        bandelet_quadtree_[root_idx] = bandelet_optimal_encoding_descriptor(wm, i0, j0, sz, threshold);
      } else {
        // Compute children, then check wether to merge or split
        const std::size_t child_sz = sz/2;
        bandelet_compute_optimal_quadtree(root_idx*4+1, wm, i0         , j0         , child_sz, threshold);
        bandelet_compute_optimal_quadtree(root_idx*4+2, wm, i0+child_sz, j0         , child_sz, threshold);
        bandelet_compute_optimal_quadtree(root_idx*4+3, wm, i0         , j0+child_sz, child_sz, threshold);
        bandelet_compute_optimal_quadtree(root_idx*4+4, wm, i0+child_sz, j0+child_sz, child_sz, threshold);
        
        bandelet_quad_descriptor_t when_split = bandelet_quad_descriptor_t(Bandelet_tag_split,
                                                                           1 + 
                                                                           bandelet_quadtree_[root_idx*4+1].bit_count() +
                                                                           bandelet_quadtree_[root_idx*4+2].bit_count() +
                                                                           bandelet_quadtree_[root_idx*4+3].bit_count() +
                                                                           bandelet_quadtree_[root_idx*4+4].bit_count(),
                                                                           bandelet_quadtree_[root_idx*4+1].sqr_err() +
                                                                           bandelet_quadtree_[root_idx*4+2].sqr_err() +
                                                                           bandelet_quadtree_[root_idx*4+3].sqr_err() +
                                                                           bandelet_quadtree_[root_idx*4+4].sqr_err());
        bandelet_quad_descriptor_t when_merged = bandelet_optimal_encoding_descriptor(wm, i0, j0, sz, threshold);
        when_merged.bit_count() += 1; // cost of merge bit
          
        bandelet_quadtree_[root_idx] = (when_split.lagrangian(threshold) < when_merged.lagrangian(threshold)) ? when_split : when_merged;
      }
    }

    int32_t bandelet_dominant_direction(const float_array2_t& wm,
                                        std::size_t i0,
                                        std::size_t j0,
                                        std::size_t sz) const {
      const int32_t N_dir = bandelet_direction_count(sz);
      if (N_dir <= 0) {
        return Bandelet_tag_untransformed;
      } else if (N_dir == 1) {
        return 0;
      } else {
        // Find eigenbasis
        typedef fixed_size_square_matrix<2,float>             covariance_matrix_t;
        typedef fixed_size_square_matrix<2,float>             eigen_matrix_t;
        typedef fixed_size_vector<column_orientation,2,float> eigen_vector_t;
        eigen_vector_t mean;
        float wt_sum = 0.0f;
        for (std::size_t i=0; i<sz; ++i) {
          for (std::size_t j=0; j<sz; ++j) {
            float    wt = abs(wm(i0+i,j0+j));
            eigen_vector_t xp = wt * eigen_vector_t(float(i), float(j));
            mean += xp;
            wt_sum += wt;
          }
        }
        if (!wt_sum) {
          // Null matrix
          return Bandelet_tag_zero_block;
        } else {
          mean /= wt_sum;
          
          covariance_matrix_t cov;
          for (std::size_t i=0; i<sz; ++i) {
            for (std::size_t j=0; j<sz; ++j) {
              float    wt = abs(wm(i0+i,j0+j));
              eigen_vector_t xp = wt * eigen_vector_t(float(i), float(j));
              xp -= mean;
              xp /= wt_sum;
              cov(0,0) += xp[0]*xp[0];
              cov(0,1) += xp[0]*xp[1];
              cov(1,0) += xp[1]*xp[0];
              cov(1,1) += xp[1]*xp[1];
            }
          }
          
          eigen_vector_t eigen_values;
          eigen_matrix_t eigen_vectors;
          bool ok = false;
          cov.symmetric_sorted_eigen_in(eigen_values,
                                        eigen_vectors,
                                        &ok);
          if (!ok) {
            // Choose random basis
            eigen_vectors.to_identity();
            std::cerr << "BAD EIGENBASIS" << std::endl;
          }
          eigen_vectors.make_right_handed();
          
          eigen_vector_t eigen_axis = eigen_vectors.axis(1); //std::cerr << "EIGEN_AXIS=" << eigen_axis[0] << " " << eigen_axis[1] << std::endl;
          float   best_dot = -1.0f;
          int32_t best_direction = -1; 
          for (int32_t direction=0; direction<N_dir; ++direction) {
            const float theta = bandelet_direction_theta(direction,sz);
            const float cos_theta = std::cos(theta);
            const float sin_theta = std::sin(theta);
            const float dot = abs(cos_theta * eigen_axis[0] + sin_theta * eigen_axis[1]);
            if (dot > best_dot) {
              best_direction = direction;
              best_dot = dot;
            }
          }
          return best_direction;
        }
      }
    }
    
    void JUNK() {
      float_array2_t wm(32,32);

      bandelet_quadtree_.resize(2*wm.extent()[0]*wm.extent()[1]);
      bandelet_transformed_data_.resize(2*wm.extent()[0]*wm.extent()[1]);

      for (std::size_t i=1; i<31; ++i) {
        wm(i,i) = 1.0;
        wm(i+1,i) = 0.5;
        wm(i-1,i) = 0.5;
      }
      std::cerr << "WM" << std::endl;
      bandelet_quad_descriptor_t result = bandelet_encoding_descriptor(Bandelet_tag_untransformed, wm, 0, 0, 32, 0.1);
      std::cerr << " START= " << 32 << "x" << 32 << ": dir = " << result.tag() << " L= " << result.lagrangian(0.1) << std::endl;
      const int32_t N_dir = bandelet_direction_count(32);
      for (int32_t direction=0; direction<N_dir; ++direction) {
        bandelet_quad_descriptor_t result_i = bandelet_encoding_descriptor(direction, wm, 0, 0, 32, 0.1);
        std::cerr << 32 << "x" << 32 << ": dir = " << result_i.tag() << " L= " << result_i.lagrangian(0.1) << std::endl;
        if (result_i.lagrangian(0.1)<result.lagrangian(0.1)) {
          result = result_i;
        }
      }
      {
        int32_t direction = bandelet_dominant_direction(wm, 0, 0, 32);
        if (direction != Bandelet_tag_untransformed) {
          bandelet_quad_descriptor_t result_i = bandelet_encoding_descriptor(direction, wm, 0, 0, 32, 0.1);
          std::cerr << "GUESS: " << 32 << "x" << 32 << ": dir = " << result_i.tag() << " L= " << result_i.lagrangian(0.1) << std::endl;
          if (result_i.lagrangian(0.1)<result.lagrangian(0.1)) {
            result = result_i;
          }
        }
        std::cerr << " BEST= " << 32 << "x" << 32 << ": dir = " << result.tag() << " L= " << result.lagrangian(0.1) << std::endl;
        exit(1);
      }
    }
    
    bool is_zero_block(const float_array2_t& wm,
                       std::size_t i0,
                       std::size_t j0,
                       std::size_t sz,
                       float threshold) const {
      for (std::size_t i=0; i<sz; ++i) {
        for (std::size_t j=0; j<sz; ++j) {
          if (int32_quantize(wm(i0+i,j0+j),threshold)!=0) return false;
        }
      }
      return true;
    }
    
    bandelet_quad_descriptor_t bandelet_optimal_encoding_descriptor(const float_array2_t& wm,
                                                                    std::size_t i0,
                                                                    std::size_t j0,
                                                                    std::size_t sz,
                                                                    float threshold) {
      if (is_zero_block(wm,i0,j0,sz,threshold)) {
        bandelet_quad_descriptor_t result = bandelet_encoding_descriptor(Bandelet_tag_untransformed, wm, i0, j0, sz, threshold);
        result.tag() = Bandelet_tag_zero_block;
        result.bit_count() = 1;
        return result;
      } else {
        bandelet_quad_descriptor_t result = bandelet_encoding_descriptor(Bandelet_tag_untransformed, wm, i0, j0, sz, threshold);
        int32_t dominant_direction = bandelet_dominant_direction(wm, i0, j0, sz);
        if (dominant_direction >= 0) {
          const int32_t N_dir = bandelet_direction_count(sz);
          int32_t min_direction = std::max(dominant_direction-1, 0);
          int32_t max_direction = std::min(dominant_direction, N_dir-1);
          for (int32_t direction = min_direction; direction<=max_direction; ++direction) {
            bandelet_quad_descriptor_t result_i = bandelet_encoding_descriptor(direction, wm, i0, j0, sz, threshold);
            if (result_i.lagrangian(threshold)<result.lagrangian(threshold)) {
              result = result_i;
            }
          }
        }
        return result;      
      }
    }

    inline int32_t ac_estimated_bit_count(int32_t x) const {
      int32_t result = 1+bitops<int32_t>::log2(x>0?x+1:-x+1);
      return result;
    }

    bandelet_quad_descriptor_t bandelet_encoding_descriptor(int32_t direction,
                                                            const float_array2_t& wm,
                                                            std::size_t i0,
                                                            std::size_t j0,
                                                            std::size_t sz,
                                                            float threshold) {
      float       sqr_err     = 0.0f;
      std::size_t bit_count = 0;
      
      bandelet_forward(direction, wm, i0,j0,sz, &(bandelet_transformed_data_[0]));
      
      bit_count += ac_estimated_bit_count(direction);
      const std::size_t sz2 = sz*sz;
      for (std::size_t i=0; i<sz2; ++i) {
        const float   v_i = bandelet_transformed_data_[i];
        const int32_t q_i = int32_quantize(v_i, threshold);
        const float   v_i_tilde = int32_dequantize(q_i, threshold);
        const float   delta = v_i-v_i_tilde;
        sqr_err += delta*delta;
        bit_count += ac_estimated_bit_count(q_i);
      }
      
      //std::cerr << "   " << sz << "x" << sz << ": dir = " << direction << " thr = " << threshold << " err= " << sqr_err << " bits = " << bit_count << std::endl;
      return bandelet_quad_descriptor_t(direction,
                                        bit_count,
                                        sqr_err);
    }
    
    void bandelet_encode_optimal_quadtree(entropy_codec_t& entropy_codec,
                                          bit_context_t& quadtree_context,
                                          int_context_t& tag_context,
                                          int_context_t& coeff_context,
                                          std::size_t root_idx,
                                          const float_array2_t& wm,
                                          std::size_t i0,
                                          std::size_t j0,
                                          std::size_t sz,
                                          float threshold) {
      const int32_t tag = bandelet_quadtree_[root_idx].tag();
#if BANDELET_TRACE
      std::cerr << "  QT: " << i0 << " " << j0 << ";" << sz;
#endif
      if (tag == Bandelet_tag_split) {
        // Split
#if BANDELET_TRACE
        std::cerr << "  ENC: split " << std::endl;
#endif
        entropy_codec.encode_bit(quadtree_context, true);
        const std::size_t child_sz = sz/2;
        bandelet_encode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, root_idx*4+1, wm, i0         , j0         , child_sz, threshold);
        bandelet_encode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, root_idx*4+2, wm, i0+child_sz, j0         , child_sz, threshold);
        bandelet_encode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, root_idx*4+3, wm, i0         , j0+child_sz, child_sz, threshold);
        bandelet_encode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, root_idx*4+4, wm, i0+child_sz, j0+child_sz, child_sz, threshold);
      } else {
        // Leaf of segmentation
#if BANDELET_TRACE
        std::cerr << "  ENC: leaf " << " dir = " << tag << " data = ";
#endif
        entropy_codec.encode_bit(quadtree_context, false);
        entropy_codec.encode_int(tag_context, tag);
        if (tag != Bandelet_tag_zero_block) {
          bandelet_forward(tag, wm, i0, j0, sz, &(bandelet_transformed_data_[0]));
          const std::size_t sz2 = sz*sz;
          for (std::size_t i=0; i<sz2; ++i) {
            const float   v_i = bandelet_transformed_data_[i];
            const int32_t q_i = int32_quantize(v_i,threshold);
#if BANDELET_TRACE
            std::cerr << " " << q_i;
#endif
            entropy_codec.encode_int(coeff_context, q_i);
          }
        }
#if BANDELET_TRACE
        std::cerr << std::endl;
#endif
      }
    }

    void bandelet_decode_optimal_quadtree(entropy_codec_t& entropy_codec,
                                          bit_context_t& quadtree_context,
                                          int_context_t& tag_context,
                                          int_context_t& coeff_context,
                                          float_array2_t& wm,
                                          std::size_t i0,
                                          std::size_t j0,
                                          std::size_t sz,
                                          float threshold) {
#if BANDELET_TRACE
      std::cerr << "  QT: " << i0 << " " << j0 << ";" << sz;
#endif
      const bool is_split = entropy_codec.decode_bit(quadtree_context);
      if (is_split) {
        // Split
#if BANDELET_TRACE
        std::cerr << "  DEC: split " << std::endl;
#endif
        const std::size_t child_sz = sz/2;
        bandelet_decode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, wm, i0         , j0         , child_sz, threshold);
        bandelet_decode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, wm, i0+child_sz, j0         , child_sz, threshold);
        bandelet_decode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, wm, i0         , j0+child_sz, child_sz, threshold);
        bandelet_decode_optimal_quadtree(entropy_codec, quadtree_context, tag_context, coeff_context, wm, i0+child_sz, j0+child_sz, child_sz, threshold);
      } else {
        const int32_t theta = entropy_codec.decode_int(tag_context);
#if BANDELET_TRACE
        std::cerr << "  DEC: leaf " << " dir = " << theta << " data = ";
#endif
        if (theta == Bandelet_tag_zero_block) {
          for (std::size_t i=0; i<sz; ++i) {
            for (std::size_t j=0; j<sz; ++j) {
              wm(i0+i,j0+j) = 0;
            }
          }
        } else {
          const std::size_t sz2 = sz*sz;
          for (std::size_t i=0; i<sz2; ++i) {
            const int32_t q_i = entropy_codec.decode_int(coeff_context);
#if BANDELET_TRACE
            std::cerr << " " << q_i;
#endif
            bandelet_transformed_data_[i] = int32_dequantize(q_i, threshold);
          }
          bandelet_backward(theta, wm, i0, j0, sz, &(bandelet_transformed_data_[0]));
        }
#if BANDELET_TRACE
        std::cerr << std::endl;
#endif
      }
    }

  protected:

    std::pair<std::size_t,double> bandelet_size_and_rms(const float_array2_t& wm,
                                                        bool is_wavelet_xform,
                                                        float threshold) {
      std::vector<uint8_t> tmpbuf(8*sizeof(float)*wm.count());
      static entropy_codec_t tmpcodec; // FIXME
      tmpcodec.set_buffer(&(tmpbuf[0]), tmpbuf.size());
      tmpcodec.start_encoder();
      double actual_rms;
      std::size_t actual_size;
      bandelet_encode(tmpcodec, wm, is_wavelet_xform, threshold, &actual_size, &actual_rms);
      tmpcodec.stop_encoder();
      actual_size = tmpcodec.current_byte_count();

      return std::make_pair(actual_size, actual_rms);
    }
    
    void bandelet_select_compression_parameters(const float_array2_t& wm,
                                                bool         is_wavelet_xform,
                                                std::size_t  target_size,
                                                double       target_rms,
                                                float       *actual_qt,
                                                std::size_t *actual_size,
                                                double      *actual_rms) {
      //std::cerr << " begin select_compression(sz=" << target_size << ", rms=" << target_rms << ")" << std::endl;
      const std::size_t max_bisections = 20; // FIXME
      const std::size_t max_steps      = 1<<max_bisections;
      float stepsize_coarse =  std::max(1e-6f, 
					std::min(float(4*wm.extent()[0]*target_rms),
                                        std::max(1e-6f,
                                                 0.5f*float(amax(wm)))));
      float stepsize_fine = std::max(float(target_rms),
                                     std::max(1e-6f, stepsize_coarse/float(max_steps)));
      std::pair<std::size_t,double> sz_rms_coarse = bandelet_size_and_rms(wm, is_wavelet_xform, stepsize_coarse);
      if (sz_rms_coarse.first>=target_size) {
        // Cannot meet size constraints!
        *actual_qt   = stepsize_coarse;
        *actual_size = sz_rms_coarse.first;
        *actual_rms  = sz_rms_coarse.second;
      } else {
        std::pair<std::size_t,double> sz_rms_fine = bandelet_size_and_rms(wm, is_wavelet_xform, stepsize_fine);
        if (sz_rms_fine.first<=target_size && sz_rms_fine.second >= target_rms) {
          *actual_qt   = stepsize_fine;
          *actual_size = sz_rms_fine.first;
          *actual_rms  = sz_rms_fine.second;
        } else {
          std::size_t N_steps = max_bisections;
          while (N_steps>0 &&
                 (stepsize_coarse-stepsize_fine)>1e-6f &&
                 (sz_rms_coarse.first-sz_rms_fine.first)>4 &&
                 (sz_rms_coarse.second-sz_rms_fine.second)>0.1*target_rms) {
            // FIXME
            float t = 0.5f;
            if (sz_rms_fine.second < target_rms && target_rms < sz_rms_coarse.second) {
              t = (float)(median((target_rms - sz_rms_fine.second) / (sz_rms_coarse.second-sz_rms_fine.second), 0.1, 0.9));
            }
            // FIXME
            float stepsize_mid = stepsize_fine + t * (stepsize_coarse-stepsize_fine);
            std::pair<std::size_t,double> sz_rms_mid = bandelet_size_and_rms(wm, is_wavelet_xform, stepsize_mid);
            if (sz_rms_mid.first>=target_size || sz_rms_mid.second < target_rms) {
              // Too big or too precise, coarsen
              sz_rms_fine = sz_rms_mid;
              stepsize_fine = stepsize_mid;
            } else {
              sz_rms_coarse = sz_rms_mid;
              stepsize_coarse = stepsize_mid;
            }
            //std::cerr << "STEP " << N_steps << " RANGE: [" << stepsize_coarse << " sz = " << sz_rms_coarse.first << " rms = " << sz_rms_coarse.second << "] [" << stepsize_fine << " sz = " << sz_rms_fine.first << " rms = " << sz_rms_fine.second << "] " << std::endl;
            --N_steps;
          }
          if (sz_rms_fine.first<=target_size) {
            // Choose fine
            *actual_qt   = stepsize_fine;
            *actual_size = sz_rms_fine.first;
            *actual_rms  = sz_rms_fine.second;
          } else {
            // Choose coarse
            *actual_qt   = stepsize_coarse;
            *actual_size = sz_rms_coarse.first;
            *actual_rms  = sz_rms_coarse.second;
          }
        }
        //std::cerr << " end selec_compression(sz=" << target_size << ", rms=" << target_rms << ") => T=" << *actual_qt << ", sz=" << *actual_size << ", rms=" << *actual_rms << std::endl; 
      }
    }
    
  protected:

    /**
     *  Compress array and encode it in given buffer using given transform kind and
     *  quantization threshold.
     *  Compression terminates when the size equals or exceeds target size or when the root mean
     *  square error falls below the target value. The array is transformed by the current
     *  transform kind and then compressed by a the geometric bandelet method.
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
      const uint16_t h_wavelet = pow2size(h);
      const uint16_t w_wavelet = pow2size(w);
      
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

      // Wavelet transform
      bool is_wavelet_xform = tkind != TRANSFORM_KIND_IDENTITY;
      if (is_wavelet_xform) {
        // Extend to power of 2 dimension
        symmetric_extend_top_left(xform, h, w);
        transform_table_[tkind]->forward2d(xform,h_wavelet,w_wavelet);
        // Clear coefficients referring purely to virtual non pow2 part
        // FIXME clear_extra_wavelet_coefficients(xform, h, w, transform_table_[tkind]);
      }

      // Select bandelet stepsize
      float bandelet_stepsize = 0.0f;
      compressed_data_header hdr(h, w, tkind, signal_mean, bandelet_stepsize); // To estimate size
      const std::size_t max_header_bytes     = hdr.compressed_size();
      const std::size_t max_code_bytes       = buf_size - max_header_bytes; 
      const std::size_t target_code_size     = target_size > max_header_bytes ? target_size-max_header_bytes : 0;
      std::size_t       bandelet_code_size;
      double            bandelet_rms;
      bandelet_select_compression_parameters(xform,
                                             is_wavelet_xform,
                                             target_code_size,
                                             target_rms,
                                             &bandelet_stepsize,
                                             &bandelet_code_size,
                                             &bandelet_rms);
      // Encode header
      uint8_t* header_buffer = (uint8_t*)buf;
      hdr = compressed_data_header(h, w, tkind, signal_mean, bandelet_stepsize);
      hdr.encode(header_buffer);
      const std::size_t header_bytes   = hdr.compressed_size();
      assert(header_bytes == max_header_bytes);

      // Encode bandelet transform
      uint8_t* code_buffer = header_buffer + header_bytes;
      entropy_codec_.set_buffer(code_buffer, max_code_bytes);
      entropy_codec_.start_encoder();
      std::size_t encode_size;
      double encode_rms;
      bandelet_encode(entropy_codec_, xform, is_wavelet_xform, bandelet_stepsize, &encode_size, &encode_rms);
      entropy_codec_.stop_encoder();
      std::size_t code_size = entropy_codec_.current_byte_count();
      assert(code_size == bandelet_code_size);

      // Return results
      *actual_size = header_bytes + code_size;
      *actual_rms  = bandelet_rms;
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

      const transform_kind_t tkind = hdr.transform_kind();
      const std::size_t h_wavelet = pow2size(hdr.h());
      const std::size_t w_wavelet = pow2size(hdr.w());
      const bool is_wavelet_xform = tkind != TRANSFORM_KIND_IDENTITY;
      
      float_array2_t xform(h_wavelet,w_wavelet);
      entropy_codec_.set_buffer(const_cast<uint8_t*>(code_buffer), code_bytes);
      entropy_codec_.start_decoder();
      bandelet_decode(entropy_codec_, xform, is_wavelet_xform, hdr.xform_threshold());
      entropy_codec_.stop_decoder();
      
      // Invert wavelet transform
      if (is_wavelet_xform) {
        transform_table_[tkind]->backward2d(xform,h_wavelet,w_wavelet);
      }
      
      // Back-transform from zero-centered range
      array.resize(subscript_t(hdr.h(),hdr.w()));
      for (std::size_t i=0; i<hdr.h(); ++i) {
        for (std::size_t j=0; j<hdr.w(); ++j) {
          array(i,j) = xform(i,j) + hdr.signal_mean();
        }
      }
    }


  };

} // namespace sl


#endif
