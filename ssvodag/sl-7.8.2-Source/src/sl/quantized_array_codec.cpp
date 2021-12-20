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
#include <sl/quantized_array_codec.hpp>

namespace sl {

  namespace quantized_array_codec_detail {

    // ============= Types
    
    typedef array_codec::int8_array2_t   int8_array2_t;
    typedef array_codec::uint8_array2_t  uint8_array2_t;
    typedef array_codec::int16_array2_t  int16_array2_t;
    typedef array_codec::uint16_array2_t uint16_array2_t;
    typedef array_codec::int32_array2_t  int32_array2_t;
    typedef array_codec::uint32_array2_t uint32_array2_t;
    typedef array_codec::float_array2_t  float_array2_t;
    typedef float_array2_t::subscript_t  subscript_t;

    typedef quantized_array_codec::bitio_t bitio_t;

    // ============== Encoding/decoding
    
    static inline int32_t int32_delta_predict(const int32_array2_t& x, std::size_t i, std::size_t j) {
      int32_t result = 0;
      if (i>0 && j>0) {
	const int16_t a = x(i,j-1);
	const int16_t b = x(i-1,j);
	const int16_t c = x(i-1,j-1);
	
	const int16_t max_b_a = std::max(b, a);
	const int16_t min_b_a = std::min(b, a);
	if (c >= max_b_a) {
	  result = min_b_a;
	} else if (c <= min_b_a) {
	  result = max_b_a;
	} else {        
	  result = a+b-c;
	}
      } else if (j>0) {
	result = x(i,j-1);
      } else if (i>0) {
	result = x(i-1,j);
      } else {
	result = 0;
      }

      return result;
    }
    
    static void int32_delta_encode_in(int32_array2_t& delta_x, const int32_array2_t& x) {
      delta_x.resize(x.extent());
      const std::size_t h = x.extent()[0];
      const std::size_t w = x.extent()[1];
      
      for (std::size_t i=0; i<h; ++i) {
	for (std::size_t j=0; j<w; ++j) {
	  delta_x(i,j) = x(i,j)-int32_delta_predict(x,i,j);
	}
      }
    }

    static void int32_delta_decode_in(int32_array2_t& x, const int32_array2_t& delta_x) {
      x.resize(delta_x.extent());
      const std::size_t h = x.extent()[0];
      const std::size_t w = x.extent()[1];
      
      for (std::size_t i=0; i<h; ++i) {
	for (std::size_t j=0; j<w; ++j) {
	  x(i,j) = delta_x(i,j)+int32_delta_predict(x,i,j);
	}
      }
    }

    static bool int32_array_is_zero(const int32_array2_t& A,
                                    uint32_t i0,
                                    uint32_t j0,
                                    uint32_t isz,
                                    uint32_t jsz) {
      bool result = true;
      for (uint32_t i=i0; result && (i<i0+isz); ++i) {
        for (uint32_t j=j0; result && (j<j0+jsz); ++j) {
          result = (A(i,j) == 0);
        }
      }
      return result;
    }
    
    static void int32_array_quadtree_structure_compute(std::size_t idx,
						       std::vector<bool>& is_zero_block,
						       const int32_array2_t& A,
						       uint32_t i0,
						       uint32_t j0,
						       uint32_t isz,
						       uint32_t jsz) {
      assert(idx < is_zero_block.size());
      
      uint32_t isz_lo = isz/2;
      uint32_t isz_hi = isz-isz_lo;
      uint32_t jsz_lo = jsz/2;
      uint32_t jsz_hi = jsz-jsz_lo;

      if (isz_lo<2 || isz_hi<2) {
	if (jsz_lo<2 || jsz_hi<2) {
	  // Block is small, becomes leaf
	  is_zero_block[idx] = int32_array_is_zero(A, i0, j0, isz, jsz);
	} else {
	  isz_lo = isz;
	  int32_array_quadtree_structure_compute(4*idx+1+0, is_zero_block, A, i0,        j0,        isz_lo, jsz_lo);
	  int32_array_quadtree_structure_compute(4*idx+1+1, is_zero_block, A, i0,        j0+jsz_lo, isz_lo, jsz_hi);
	  is_zero_block[idx] = (is_zero_block[4*idx+1+0] && is_zero_block[4*idx+1+1]);
	}
      } else if(jsz_lo<2 || jsz_hi<2) {
	jsz_lo = jsz;
	int32_array_quadtree_structure_compute(4*idx+1+0, is_zero_block, A, i0,        j0,        isz_lo, jsz_lo);
	int32_array_quadtree_structure_compute(4*idx+1+1, is_zero_block, A, i0+isz_lo, j0,        isz_hi, jsz_lo);
	is_zero_block[idx] = (is_zero_block[4*idx+1+0] && is_zero_block[4*idx+1+1]);
      } else {
	int32_array_quadtree_structure_compute(4*idx+1+0, is_zero_block, A, i0,        j0,        isz_lo, jsz_lo);
	int32_array_quadtree_structure_compute(4*idx+1+1, is_zero_block, A, i0+isz_lo, j0,        isz_hi, jsz_lo);
	int32_array_quadtree_structure_compute(4*idx+1+2, is_zero_block, A, i0,        j0+jsz_lo, isz_lo, jsz_hi);
	int32_array_quadtree_structure_compute(4*idx+1+3, is_zero_block, A, i0+isz_lo, j0+jsz_lo, isz_hi, jsz_hi);
	is_zero_block[idx] = (is_zero_block[4*idx+1+0] && 
			      is_zero_block[4*idx+1+1] &&
			      is_zero_block[4*idx+1+2] &&
			      is_zero_block[4*idx+1+3]);
      }
    }

					    
    static void int32_array_quadtree_encode(std::size_t idx,
					    std::vector<bool>& is_zero_block,
					    bitio_t& ec,
                                            const int32_array2_t& A,
                                            uint32_t i0,
                                            uint32_t j0,
                                            uint32_t isz,
                                            uint32_t jsz) {
      assert(idx < is_zero_block.size());
      bool is_zero_block_bit = is_zero_block[idx];

      assert(is_zero_block_bit == int32_array_is_zero(A, i0, j0, isz, jsz));

      ec.encode_bit(is_zero_block_bit);
      if (is_zero_block_bit) {
        //std::cerr << "ENC: ZB" << " [ " << i0 << " " << j0 << "] " <<" [ " << i0+isz << " " << j0+jsz << "] " << std::endl;
      } else {
        uint32_t isz_lo = isz/2;
        uint32_t isz_hi = isz-isz_lo;
        uint32_t jsz_lo = jsz/2;
        uint32_t jsz_hi = jsz-jsz_lo;

        if (isz_lo<2 || isz_hi<2) {
          if (jsz_lo<2 || jsz_hi<2) {
            //std::cerr << "ENC: B" << " [ " << i0 << " " << j0 << "] " <<" [ " << i0+isz << " " << j0+jsz << "] " << std::endl;
            for (uint32_t i=i0; i<i0+isz; ++i) {
              for (uint32_t j=j0; j<j0+jsz; ++j) {
                ec.encode_int(A(i,j)); // FIXME
              }
            }
          } else {
            isz_lo = isz;
            int32_array_quadtree_encode(4*idx+1+0, is_zero_block, ec, A, i0,        j0,        isz_lo, jsz_lo);
            int32_array_quadtree_encode(4*idx+1+1, is_zero_block, ec, A, i0,        j0+jsz_lo, isz_lo, jsz_hi);
          }
        } else if(jsz_lo<2 || jsz_hi<2) {
          jsz_lo = jsz;
          int32_array_quadtree_encode(4*idx+1+0, is_zero_block, ec, A, i0,        j0,        isz_lo, jsz_lo);
          int32_array_quadtree_encode(4*idx+1+1, is_zero_block, ec, A, i0+isz_lo, j0,        isz_hi, jsz_lo);
        } else {
          int32_array_quadtree_encode(4*idx+1+0, is_zero_block, ec, A, i0,        j0,        isz_lo, jsz_lo);
          int32_array_quadtree_encode(4*idx+1+1, is_zero_block, ec, A, i0+isz_lo, j0,        isz_hi, jsz_lo);
          int32_array_quadtree_encode(4*idx+1+2, is_zero_block, ec, A, i0,        j0+jsz_lo, isz_lo, jsz_hi);
          int32_array_quadtree_encode(4*idx+1+3, is_zero_block, ec, A, i0+isz_lo, j0+jsz_lo, isz_hi, jsz_hi);
        }
      }
    }

    static void int32_array_quadtree_encode(bitio_t& ec,
                                            const int32_array2_t& A) {
      std::vector<bool> is_zero_block(A.count()*8); // Too large than needed...
      int32_array_quadtree_structure_compute(0, is_zero_block, A, 0, 0, A.extent()[0], A.extent()[1]);
      int32_array_quadtree_encode(0, is_zero_block, ec, A, 0, 0, A.extent()[0], A.extent()[1]);
    }

    static void int32_array_quadtree_decode(bitio_t& ec,
                                            int32_array2_t& A,
                                            uint32_t i0,
                                            uint32_t j0,
                                            uint32_t isz,
                                            uint32_t jsz) {
      bool is_zero_block = ec.decode_bit();
      if (is_zero_block) {
        //std::cerr << "DEC: ZB" << " [ " << i0 << " " << j0 << "] " <<" [ " << i0+isz << " " << j0+jsz << "] " << std::endl;
        for (uint32_t i=i0; i<i0+isz; ++i) {
          for (uint32_t j=j0; j<j0+jsz; ++j) {
            A(i,j) = 0; // FIXME
          }
        }
      } else {
        uint32_t isz_lo = isz/2;
        uint32_t isz_hi = isz-isz_lo;
        uint32_t jsz_lo = jsz/2;
        uint32_t jsz_hi = jsz-jsz_lo;

        if (isz_lo<2 || isz_hi<2) {
          if (jsz_lo<2 || jsz_hi<2) {
            //std::cerr << "DEC: B" << " [ " << i0 << " " << j0 << "] " <<" [ " << i0+isz << " " << j0+jsz << "] " << std::endl;
            for (uint32_t i=i0; i<i0+isz; ++i) {
              for (uint32_t j=j0; j<j0+jsz; ++j) {
                A(i,j) = ec.decode_int(); // FIXME
              }
            }
          } else {
            isz_lo = isz;
            int32_array_quadtree_decode(ec, A, i0,        j0,        isz_lo, jsz_lo);
            int32_array_quadtree_decode(ec, A, i0,        j0+jsz_lo, isz_lo, jsz_hi);
          }
        } else if (jsz_lo<2 || jsz_hi<2) {
          jsz_lo = jsz;
          int32_array_quadtree_decode(ec, A, i0,        j0,        isz_lo, jsz_lo);
          int32_array_quadtree_decode(ec, A, i0+isz_lo, j0,        isz_hi, jsz_lo);
        } else {
          int32_array_quadtree_decode(ec, A, i0,        j0,        isz_lo, jsz_lo);
          int32_array_quadtree_decode(ec, A, i0+isz_lo, j0,        isz_hi, jsz_lo);
          int32_array_quadtree_decode(ec, A, i0,        j0+jsz_lo, isz_lo, jsz_hi);
          int32_array_quadtree_decode(ec, A, i0+isz_lo, j0+jsz_lo, isz_hi, jsz_hi);
        }
      }
    }

    // ================= Quantizing/Dequantizing
    
    template <class int_t>
    static void int_array_quantize_in(int32_array2_t& Aq,
                                      const sl::dense_array<int_t, 2, void>& A,
                                      int32_t eps) {
      const int32_t two_eps_plus_one = 2*eps+1;
      const uint32_t isz = A.extent()[0];
      const uint32_t jsz = A.extent()[1];
      for (uint32_t i=0; i<isz; ++i) {
        for (uint32_t j=0; j<jsz; ++j) {
          int32_t A_ij = int32_t(A(i,j));
          int32_t Aq_ij = (A_ij>0) ? ((A_ij+eps)/two_eps_plus_one) : ((A_ij-eps)/two_eps_plus_one);
          Aq(i,j) = Aq_ij;
        }
      }
    }
    
    template <class int_t>
    static void int_array_dequantize_in(sl::dense_array<int_t, 2, void>& A,
                                        const int32_array2_t& Aq,
                                        int32_t eps) {
      const int_t imin = std::numeric_limits<int_t>::min();
      const int_t imax = std::numeric_limits<int_t>::max();

      const int32_t two_eps_plus_one = 2*eps+1;
      const uint32_t isz = A.extent()[0];
      const uint32_t jsz = A.extent()[1];
      for (uint32_t i=0; i<isz; ++i) {
        for (uint32_t j=0; j<jsz; ++j) {
          int32_t Aq_ij = Aq(i,j);
          int_t A_ij = int_t(sl::median(Aq_ij*two_eps_plus_one,
                                        int32_t(imin),
                                        int32_t(imax)));                                        
          A(i,j) = A_ij;
        }
      }
    }

    static void float_array_quantize_in(int32_array2_t& Aq,
                                        const float_array2_t& A,
                                        float eps) {
      const uint32_t isz = A.extent()[0];
      const uint32_t jsz = A.extent()[1];

      const float one_over_two_eps = 0.5f/eps;
      for (uint32_t i=0; i<isz; ++i) {
        for (uint32_t j=0; j<jsz; ++j) {
          float A_ij = A(i,j)*one_over_two_eps;
          int32_t Aq_ij = (A_ij>0) ? int32_t(A_ij+0.5f) : int32_t(A_ij-0.5f);
          Aq(i,j) = Aq_ij;
        }
      }
    }
    
    static void float_array_dequantize_in(float_array2_t& A,
                                          const int32_array2_t& Aq,                                          
                                          float eps) {
      const uint32_t isz = A.extent()[0];
      const uint32_t jsz = A.extent()[1];

      const float two_eps = 2.0f*eps;
      for (uint32_t i=0; i<isz; ++i) {
        for (uint32_t j=0; j<jsz; ++j) {
          int32_t Aq_ij = Aq(i,j);
          float   A_ij = float(Aq_ij)*two_eps;
          A(i,j) = A_ij;
        }
      }
    }

    
    // ================= Compresing/Decompresing: integers
    
    template <class int_t>
    static void int_array_compress_amax(bitio_t& ec,
                                        const sl::dense_array<int_t, 2, void>& A,
                                        std::size_t target_size,
                                        double target_error,
                                        void* buf,
                                        std::size_t buf_size,
                                        std::size_t *actual_size, 
                                        double *actual_error,
					bool is_actual_error_rms,
					bool is_compressing_header,
					bool is_delta_encoding) {
      typedef sl::dense_array<int_t, 2, void> int_array2_t;

      const uint32_t h = A.extent()[0];
      const uint32_t w = A.extent()[1];

      int32_t tol    = int32_t(target_error+0.1);
      int32_t tolmax = amax(A);
      if (tol>tolmax) tol = tolmax;
      
      int32_array2_t       Aq(h,w);
      std::vector<uint8_t> tmpbufdata(h*w*16); // Large
      
      bool done = false;
      while (!done) {
        //std::cerr << "Compress try with " << tol << std::endl;
        // Quantize with given amax tolerance
        int_array_quantize_in(Aq, A, tol);

        // Encode to large tmp buffer
        uint8_t* tmpbuf = &(tmpbufdata[0]); 
	if (is_compressing_header) {
	  const bool is_square_matrix = h == w;
	  const bool h_gt_8bit = h > 0xff;
	  const bool w_gt_8bit = w > 0xff;
	  const bool tol_gt_8bit = (uint32_t(tol) > 0xff);
	  const bool tol_gt_16bit = (uint32_t(tol) > 0xffff);
	  uint8_t hdr = 0;
	  if (is_square_matrix) hdr |= 0x80;
	  if (h_gt_8bit)        hdr |= 0x40;
	  if (w_gt_8bit)        hdr |= 0x20;
	  if (tol_gt_8bit)      hdr |= 0x10;
	  if (tol_gt_16bit)     hdr |= 0x08;

	  tmpbuf += be_encode(hdr, tmpbuf);
 	  if (h_gt_8bit) { 
	    tmpbuf += be_encode(uint16_t(h), tmpbuf); 
	  } else { 
	    tmpbuf += be_encode(uint8_t(h), tmpbuf); 
	  }
	  if (!is_square_matrix) {
	    if (w_gt_8bit) { 
	      tmpbuf += be_encode(uint16_t(w), tmpbuf); 
	    } else { 
	      tmpbuf += be_encode(uint8_t(w), tmpbuf); 
	    }
	  } 
	  if (tol_gt_16bit) { 
	    tmpbuf += be_encode(uint32_t(tol), tmpbuf); 
	  } else if (tol_gt_8bit) { 
	    tmpbuf += be_encode(uint16_t(tol), tmpbuf); 
	  } else {
	    tmpbuf += be_encode(uint8_t(tol), tmpbuf); 
	  }	    
	} else {
	  tmpbuf += be_encode(uint16_t(h),   tmpbuf);
	  tmpbuf += be_encode(uint16_t(w),   tmpbuf);
	  tmpbuf += be_encode(uint32_t(tol), tmpbuf);
	}
        std::size_t header_size = tmpbuf - &(tmpbufdata[0]);
        ec.set_buffer(tmpbuf, tmpbufdata.size()-header_size);
        ec.start_encoder();
        {
	  if (is_delta_encoding) {
	    int32_array2_t       delta_Aq(h,w);
	    int32_delta_encode_in(delta_Aq, Aq);
	    int32_array_quadtree_encode(ec, delta_Aq);
	  } else {
	    int32_array_quadtree_encode(ec, Aq);
	  }
        }
        ec.stop_encoder();

        // Check whether we fit into buf
        std::size_t total_size = header_size + ec.current_byte_count();
        done = ((tol == tolmax) || total_size <= target_size);
        if (done) {
          //std::cerr << "Converged with " << tol << std::endl;
          if (total_size > buf_size) {
            SL_FAIL("Buffer overflow!");
          }
          // Copy to output buf
          uint8_t* obuf = static_cast<uint8_t*>(buf);
          for (std::size_t i=0; i<total_size; ++i) {
            obuf[i] = tmpbufdata[i];
          }
          if (actual_size) *actual_size = total_size;
          if (actual_error) {
            // FIXME - could be faster
            int_array2_t A_prime(h,w);
            int_array_dequantize_in(A_prime, Aq, tol);
            *actual_error = (is_actual_error_rms) ? rms(A, A_prime) : amax_diff(A, A_prime);
          }
        } else {
          // Try next tol
          tol *= 2;
          if (tol>tolmax) tol = tolmax;
        }
      }
    }
    
    template <class int_t>
    static void int_array_compress_rms(bitio_t& ec,
                                       const sl::dense_array<int_t, 2, void>& A,
                                       std::size_t target_size,
                                       double target_error,
                                       void* buf,
                                       std::size_t buf_size,
                                       std::size_t *actual_size, 
                                       double *actual_error,
				       bool is_compressing_header,
				       bool is_delta_encoding) {
      typedef sl::dense_array<int_t, 2, void> int_array2_t;
      
      const uint32_t h = A.extent()[0];
      const uint32_t w = A.extent()[1];

      int32_t tolmin = int32_t(target_error+0.1);
      int32_t tolmax = amax(A);
      if (tolmin>tolmax) tolmin = tolmax;
      if (tolmin == 0) {
        // Lossless!
        int_array_compress_amax(ec,
                                A,
                                target_size,
                                target_error,
                                buf, buf_size,
                                actual_size, actual_error, true, 
				is_compressing_header,
				is_delta_encoding);
      } else {
        // Lossy
        
        int32_array2_t Aq(h,w);
        int_array2_t   A_prime(h,w);

        int_array_quantize_in(Aq, A, tolmin);
        int_array_dequantize_in(A_prime, Aq, tolmin);
        double rmse_tolmin = rms(A, A_prime);

        int_array_quantize_in(Aq, A, tolmax);
        int_array_dequantize_in(A_prime, Aq, tolmax);
        double rmse_tolmax = rms(A, A_prime);

        if (rmse_tolmax<target_error) {
          tolmin = tolmax;
          rmse_tolmin = rmse_tolmax;
        }
        while (tolmax>tolmin+1) {
          //std::cerr << "tolmin = " << tolmin << " tolmax = " << tolmax << std::endl;
          int32_t tolmid = (tolmax+tolmin)/2;
          int_array_quantize_in(Aq, A, tolmid);
          int_array_dequantize_in(A_prime, Aq, tolmid);
          double rmse_tolmid = rms(A, A_prime);
          if (rmse_tolmid>=target_error) {
            tolmax = tolmid;
            rmse_tolmax = rmse_tolmid;
          } else {
            tolmin = tolmid;
            rmse_tolmin = rmse_tolmid;
          }
        }
	
	SL_USEVAR(rmse_tolmin);
	
        // Compress using found tolerance
        int_array_compress_amax(ec,
                                A,
                                target_size,
                                double(tolmin),
                                buf, buf_size,
                                actual_size, actual_error, true,
				is_compressing_header,
				is_delta_encoding);
      }
    }

    
    template <class int_t>
    static void int_array_decompress(bitio_t& ec,
                                     sl::dense_array<int_t, 2, void>& A,
                                     const void* buf,
                                     std::size_t buf_size,
				     bool is_compressing_header,
				     bool is_delta_encoding) {
      // Decode array shape
      uint8_t* ibuf = static_cast<uint8_t*>(const_cast<void*>(buf));
      std::size_t h;
      std::size_t w;
      uint32_t tol;
      if (is_compressing_header) {
        const uint8_t hdr = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
        const bool is_square_matrix  = (hdr & 0x80)!=0;
	const bool h_gt_8bit =  (hdr & 0x40)!=0;
	const bool w_gt_8bit = (hdr & 0x20)!=0;
	const bool tol_gt_8bit = (hdr & 0x10)!=0;
	const bool tol_gt_16bit = (hdr & 0x08)!=0;

	if (h_gt_8bit) { 
	  h = be_decode<uint16_t>(ibuf); ibuf += sizeof(uint16_t);
	} else { 
	  h = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
	}
	if (is_square_matrix) {
	  w = h;
	} else if (w_gt_8bit) { 
	  w = be_decode<uint16_t>(ibuf); ibuf += sizeof(uint16_t);
	} else {
	  w = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
	} 
	if (tol_gt_16bit) { 
	  tol = be_decode<uint32_t>(ibuf); ibuf += sizeof(uint32_t);
	} else if (tol_gt_8bit) { 
	  tol = be_decode<uint16_t>(ibuf); ibuf += sizeof(uint16_t);
	} else {
	  tol = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
	}	    
      } else {
	h   = be_decode<uint16_t>(ibuf); ibuf+= sizeof(uint16_t);
	w   = be_decode<uint16_t>(ibuf); ibuf+= sizeof(uint16_t);
	tol = be_decode<uint32_t>(ibuf); ibuf+= sizeof(uint32_t);
      }

      std::size_t header_size = ibuf - static_cast<const uint8_t*>(buf);

      // Decode quantized
      int32_array2_t Aq(h,w);
      ec.set_buffer(ibuf, buf_size-header_size);
      ec.start_decoder();
      {
	if (is_delta_encoding) {
	  int32_array2_t delta_Aq(h,w);
	  int32_array_quadtree_decode(ec, delta_Aq, 0, 0, h, w);
	  int32_delta_decode_in(Aq, delta_Aq);
	} else {
	  int32_array_quadtree_decode(ec, Aq, 0, 0, h, w);
	}
      }
      ec.stop_decoder();

      // Prepare array
      A.resize(subscript_t(h,w));
      int_array_dequantize_in(A, Aq, tol);
    }

    // ================= Compresing/Decompresing: floats

    static void float_array_compress_amax(bitio_t& ec,
                                          const float_array2_t& A,
                                          std::size_t target_size,
                                          double target_error,
                                          void* buf,
                                          std::size_t buf_size,
                                          std::size_t *actual_size, 
                                          double *actual_error,
					  bool is_actual_error_rms,
					  bool is_compressing_header,
					  bool is_delta_encoding) {
      //std::cerr << "Float compress amax enter" << std::endl;
      enum { Int32_quantize_bits = 20 };
      const float Int32_quantize_scale =(float)(1<<Int32_quantize_bits);

      const uint32_t h = A.extent()[0];
      const uint32_t w = A.extent()[1];

      float   tolmax = std::max(1e-6f, float(amax(A)));
      float   tolmin = std::max(float(target_error),
                                std::max(1e-6f, tolmax/Int32_quantize_scale));
      if (tolmax<tolmin) { tolmax=tolmin; }

      float    tol = tolmin;

      int32_array2_t Aq(h,w);
      std::vector<uint8_t> tmpbufdata(h*w*16); // Large

      bool done = false;
      while (!done) {
        //std::cerr << "Float Compress try with " << tol << std::endl;
        // Quantize with given amax tolerance
        float_array_quantize_in(Aq, A, tol);

        // Encode to large tmp buffer
        uint8_t* tmpbuf = &(tmpbufdata[0]);
        if (is_compressing_header) {
	  const bool is_square_matrix = h == w;
	  const bool h_gt_8bit = h > 0xff;
	  const bool w_gt_8bit = w > 0xff;
	  uint8_t hdr = 0;
	  if (is_square_matrix) hdr |= 0x80;
	  if (h_gt_8bit)        hdr |= 0x40;
	  if (w_gt_8bit)        hdr |= 0x20;

	  tmpbuf += be_encode(hdr, tmpbuf);
 	  if (h_gt_8bit) { 
	    tmpbuf += be_encode(uint16_t(h), tmpbuf); 
	  } else { 
	    tmpbuf += be_encode(uint8_t(h), tmpbuf); 
	  }
	  if (!is_square_matrix) {
	    if (w_gt_8bit) { 
	      tmpbuf += be_encode(uint16_t(w), tmpbuf); 
	    } else { 
	      tmpbuf += be_encode(uint8_t(w), tmpbuf); 
	    }
	  } 
	  tmpbuf += be_encode(float(tol), tmpbuf);
	} else {
	  tmpbuf += be_encode(uint16_t(h),tmpbuf);
	  tmpbuf += be_encode(uint16_t(w),tmpbuf);
	  tmpbuf += be_encode(float(tol), tmpbuf);
	}
        std::size_t header_size = tmpbuf - &(tmpbufdata[0]);
        ec.set_buffer(tmpbuf, tmpbufdata.size()-header_size); 
        ec.start_encoder();
        {
	  if (is_delta_encoding) {
	    int32_array2_t       delta_Aq(h,w);
	    int32_delta_encode_in(delta_Aq, Aq);
	    int32_array_quadtree_encode(ec, delta_Aq);
	  } else {
	    int32_array_quadtree_encode(ec, Aq);
	  }
        }
        ec.stop_encoder();

        // Check whether we fit into buf
        std::size_t total_size = header_size + ec.current_byte_count();
        done = ((tol >= tolmax) ||
                (total_size <= target_size));
        if (done) {
          //std::cerr << "Converged with " << tol << std::endl;
          if (total_size > buf_size) {
            SL_FAIL("Buffer overflow!");
          }
          // Copy to output buf
          uint8_t* obuf = static_cast<uint8_t*>(buf);
          for (std::size_t i=0; i<total_size; ++i) {
            obuf[i] = tmpbufdata[i];
          }
          if (actual_size) *actual_size = total_size;
          if (actual_error) {
            // FIXME
            float_array2_t A_prime(h,w);
            float_array_dequantize_in(A_prime, Aq, tol);
            *actual_error = (is_actual_error_rms ? rms(A, A_prime) : amax_diff(A, A_prime));
          }
        } else {
          // Try next tol
          tol *= 2.0f; // FIXME
          if (tol>tolmax) tol = tolmax;
        }
      }
      //std::cerr << "Float compress amax exit" << std::endl;
    }

    static void float_array_compress_rms(bitio_t& ec,
                                         const float_array2_t& A,
                                         std::size_t target_size,
                                         double target_error,
                                         void* buf,
                                         std::size_t buf_size,
                                         std::size_t *actual_size, 
                                         double *actual_error,
					 bool is_compressing_header,
					 bool is_delta_encoding) {
      //std::cerr << "Float compress rms enter" << std::endl;
      enum { Int32_quantize_bits = 20 };
      const float Int32_quantize_scale =(float)(1<<Int32_quantize_bits);

      const uint32_t h = A.extent()[0];
      const uint32_t w = A.extent()[1];

      float   tolmax = std::max(1e-6f, float(amax(A)));
      float   tolmin = std::max(1e-6f, tolmax/Int32_quantize_scale);

      float   tol    = (float)(target_error);
      if (tol<tolmin) { tol=tolmin; }
      if (tol>tolmax) { tol=tolmax; tolmin=tolmax; }
      if (tolmin<tol) { tolmin = tol; }

      if (tolmax <= 1e-6f) {
        // Lossless!
        float_array_compress_amax(ec,
                                  A,
                                  target_size,
                                  target_error,
                                  buf, buf_size,
                                  actual_size, actual_error, true,
				  is_compressing_header,
				  is_delta_encoding);
      } else {
        // Lossy
        
        int32_array2_t Aq(h,w);
        float_array2_t A_prime(h,w);

        float_array_quantize_in(Aq, A, tolmin);
        float_array_dequantize_in(A_prime, Aq, tolmin);
        double rmse_tolmin = rms(A, A_prime);

        float_array_quantize_in(Aq, A, tolmax);
        float_array_dequantize_in(A_prime, Aq, tolmax);
        double rmse_tolmax = rms(A, A_prime);

        if (rmse_tolmax<target_error) {
          tolmin = tolmax;
          rmse_tolmin = rmse_tolmax;
        }
	
        while ((tolmax-tolmin)>0.1f*(1e-6+target_error)) {
          //std::cerr << "tolmin = " << tolmin << " tolmax = " << tolmax << std::endl;
          float tolmid = 0.5f*(tolmax+tolmin);
          float_array_quantize_in(Aq, A, tolmid);
          float_array_dequantize_in(A_prime, Aq, tolmid);
          double rmse_tolmid = rms(A, A_prime);
          if (rmse_tolmid>=target_error) {
            tolmax = tolmid;
            rmse_tolmax = rmse_tolmid;
          } else {
            tolmin = tolmid;
            rmse_tolmin = rmse_tolmid;
          }
        }

	SL_USEVAR(rmse_tolmin);

        // Compress using found tolerance
        float_array_compress_amax(ec,
                                  A,
                                  target_size,
                                  double(tolmin),
                                  buf, buf_size,
                                  actual_size, actual_error, true,
				  is_compressing_header,
				  is_delta_encoding);
      }
      //std::cerr << "Float compress rms exit" << std::endl;
    }
    
    static void float_array_decompress(bitio_t& ec,
                                       float_array2_t& A,
                                       const void* buf,
                                       std::size_t buf_size,
				       bool is_compressing_header,
				       bool is_delta_encoding) {
      // Decode array shape
      uint8_t* ibuf = static_cast<uint8_t*>(const_cast<void*>(buf));
      std::size_t h;
      std::size_t w;
      float tol;
      if (is_compressing_header) {
        const uint8_t hdr = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
        const bool is_square_matrix  = (hdr & 0x80)!=0;
	const bool h_gt_8bit =  (hdr & 0x40)!=0;
	const bool w_gt_8bit = (hdr & 0x20)!=0;

	if (h_gt_8bit) { 
	  h = be_decode<uint16_t>(ibuf); ibuf += sizeof(uint16_t);
	} else { 
	  h = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
	}
	if (is_square_matrix) {
	  w = h;
	} else if (w_gt_8bit) { 
	  w = be_decode<uint16_t>(ibuf); ibuf += sizeof(uint16_t);
	} else {
	  w = be_decode<uint8_t>(ibuf); ibuf += sizeof(uint8_t);
	} 
	tol = be_decode<float>(ibuf);    ibuf+= sizeof(float);
      } else {
	h   = be_decode<uint16_t>(ibuf); ibuf+= sizeof(uint16_t);
	w   = be_decode<uint16_t>(ibuf); ibuf+= sizeof(uint16_t);
	tol = be_decode<float>(ibuf);    ibuf+= sizeof(float);
      }
      std::size_t header_size = ibuf - static_cast<const uint8_t*>(buf);

      //std::cerr << "Float decode (1) - (" << h << " " << w << ")" << std::endl;
      // Decode quantized
      int32_array2_t Aq(h,w);
      ec.set_buffer(ibuf, buf_size-header_size);
      ec.start_decoder();
      {
	if (is_delta_encoding) {
	  int32_array2_t delta_Aq(h,w);
	  int32_array_quadtree_decode(ec, delta_Aq, 0, 0, h, w);
	  int32_delta_decode_in(Aq, delta_Aq);
	} else {
	  int32_array_quadtree_decode(ec, Aq, 0, 0, h, w);
	}
      }
      ec.stop_decoder();

      //std::cerr << "Float decode (2)" << std::endl;
      // Prepare array
      A.resize(subscript_t(h,w));
      float_array_dequantize_in(A, Aq, tol);
      //std::cerr << "Float decode exit" << std::endl;
    }

  } // namespace


  // ====================================== Class implementation

  
  void quantized_array_codec::float_compress(const float_array2_t& array,
                                             std::size_t target_size,
                                             double target_rms_error,
                                             void* buf,
                                             std::size_t buf_size,
                                             std::size_t* actual_size,
                                             double* actual_rms_error) {
    quantized_array_codec_detail::float_array_compress_rms(bitio_,
                                                           array,
                                                           target_size,
                                                           target_rms_error,
                                                           buf, buf_size,
                                                           actual_size,
                                                           actual_rms_error,
							   is_compressing_header_, is_delta_encoding_);
  }

  void quantized_array_codec::float_compress_amax(const float_array2_t& array,
                                                  std::size_t target_size,
                                                  double target_amax_error,
                                                  void* buf,
                                                  std::size_t buf_size,
                                                  std::size_t* actual_size,
                                                  double* actual_amax_error) {
    quantized_array_codec_detail::float_array_compress_amax(bitio_,
                                                            array,
                                                            target_size,
                                                            target_amax_error,
                                                            buf, buf_size,
                                                            actual_size,
                                                            actual_amax_error, false, 
							    is_compressing_header_, is_delta_encoding_);
  }

  void quantized_array_codec::float_decompress(float_array2_t& array,
                                               const void* buf,
                                               std::size_t buf_size) {
    quantized_array_codec_detail::float_array_decompress(bitio_,
                                                         array,
                                                         buf, buf_size,
							 is_compressing_header_, is_delta_encoding_);
  }

    
  void quantized_array_codec::int8_compress(const int8_array2_t& array,
                                            std::size_t target_size,
                                            double target_rms_error,
                                            void* buf,
                                            std::size_t buf_size,
                                            std::size_t* actual_size,
                                            double* actual_rms_error) {
    quantized_array_codec_detail::int_array_compress_rms(bitio_,
                                                         array,
                                                         target_size,
                                                         target_rms_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_rms_error,
							 is_compressing_header_, is_delta_encoding_);
  }

  void quantized_array_codec::int8_compress_amax(const int8_array2_t& array,
                                                 std::size_t target_size,
                                                 double target_amax_error,
                                                 void* buf,
                                                 std::size_t buf_size,
                                                 std::size_t* actual_size,
                                                 double* actual_amax_error) {
    quantized_array_codec_detail::int_array_compress_amax(bitio_,
							  array,
							  target_size,
							  target_amax_error,
							  buf, buf_size,
							  actual_size,
							  actual_amax_error, false,
							  is_compressing_header_, is_delta_encoding_);
  }

  void quantized_array_codec::int8_decompress(int8_array2_t& array,
                                              const void* buf,
                                              std::size_t buf_size) {
    quantized_array_codec_detail::int_array_decompress(bitio_,
                                                       array,
                                                       buf, buf_size,
						       is_compressing_header_, is_delta_encoding_);
  }

  void quantized_array_codec::uint8_compress(const uint8_array2_t& array,
                                             std::size_t target_size,
                                             double target_rms_error,
                                             void* buf,
                                             std::size_t buf_size,
                                             std::size_t* actual_size,
                                             double* actual_rms_error) {
    quantized_array_codec_detail::int_array_compress_rms(bitio_,
                                                         array,
                                                         target_size,
                                                         target_rms_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_rms_error,
							 is_compressing_header_, is_delta_encoding_);

  }

  void quantized_array_codec::uint8_compress_amax(const uint8_array2_t& array,
                                                  std::size_t target_size,
                                                  double target_amax_error,
                                                  void* buf,
                                                  std::size_t buf_size,
                                                  std::size_t* actual_size,
                                                  double* actual_amax_error) {
    quantized_array_codec_detail::int_array_compress_amax(bitio_,
                                                         array,
                                                         target_size,
                                                         target_amax_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_amax_error, false,
							  is_compressing_header_, is_delta_encoding_);
  }
  
  void quantized_array_codec::uint8_decompress(uint8_array2_t& array,
                                               const void* buf,
                                               std::size_t buf_size) {
    quantized_array_codec_detail::int_array_decompress(bitio_,
                                                       array,
                                                       buf, buf_size,
						       is_compressing_header_, is_delta_encoding_);
  }

    
  void quantized_array_codec::int16_compress(const int16_array2_t& array,
                                             std::size_t target_size,
                                             double target_rms_error,
                                             void* buf,
                                             std::size_t buf_size,
                                             std::size_t* actual_size,
                                             double* actual_rms_error)  {
    quantized_array_codec_detail::int_array_compress_rms(bitio_,
                                                         array,
                                                         target_size,
                                                         target_rms_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_rms_error,
							 is_compressing_header_, is_delta_encoding_);
  }

  void quantized_array_codec::int16_compress_amax(const int16_array2_t& array,
                                                  std::size_t target_size,
                                                  double target_amax_error,
                                                  void* buf,
                                                  std::size_t buf_size,
                                                  std::size_t* actual_size,
                                                  double* actual_amax_error) {
    quantized_array_codec_detail::int_array_compress_amax(bitio_,
                                                         array,
                                                         target_size,
                                                         target_amax_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_amax_error, false,
							  is_compressing_header_, is_delta_encoding_);
  }
  
  
  void quantized_array_codec::int16_decompress(int16_array2_t& array,
                                               const void* buf,
                                               std::size_t buf_size) {
    quantized_array_codec_detail::int_array_decompress(bitio_,
                                                       array,
                                                       buf, buf_size,
						       is_compressing_header_, is_delta_encoding_);
  }
  
  void quantized_array_codec::uint16_compress(const uint16_array2_t& array,
                                              std::size_t target_size,
                                              double target_rms_error,
                                              void* buf,
                                              std::size_t buf_size,
                                              std::size_t* actual_size,
                                              double* actual_rms_error) {
    quantized_array_codec_detail::int_array_compress_rms(bitio_,
                                                         array,
                                                         target_size,
                                                         target_rms_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_rms_error,
							 is_compressing_header_, is_delta_encoding_);
  }
  

    
  void quantized_array_codec::uint16_compress_amax(const uint16_array2_t& array,
                                                   std::size_t target_size,
                                                   double target_amax_error,
                                                   void* buf,
                                                   std::size_t buf_size,
                                                   std::size_t* actual_size,
                                                   double* actual_amax_error)  {
    quantized_array_codec_detail::int_array_compress_amax(bitio_,
                                                         array,
                                                         target_size,
                                                         target_amax_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_amax_error, false,
							  is_compressing_header_, is_delta_encoding_);
  }
  

  void quantized_array_codec::uint16_decompress(uint16_array2_t& array,
                                                const void* buf,
                                                std::size_t buf_size) {
    quantized_array_codec_detail::int_array_decompress(bitio_,
                                                       array,
                                                       buf, buf_size,
						       is_compressing_header_, is_delta_encoding_);
  }
  
  void quantized_array_codec::int32_compress(const int32_array2_t& array,
                                             std::size_t target_size,
                                             double target_rms_error,
                                             void* buf,
                                             std::size_t buf_size,
                                             std::size_t* actual_size,
                                             double* actual_rms_error) {
    quantized_array_codec_detail::int_array_compress_rms(bitio_,
                                                         array,
                                                         target_size,
                                                         target_rms_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_rms_error,
							 is_compressing_header_, is_delta_encoding_);
  }
  
  void quantized_array_codec::int32_compress_amax(const int32_array2_t& array,
                                                  std::size_t target_size,
                                                  double target_amax_error,
                                                  void* buf,
                                                  std::size_t buf_size,
                                                  std::size_t* actual_size,
                                                  double* actual_amax_error) {
    quantized_array_codec_detail::int_array_compress_amax(bitio_,
                                                          array,
                                                          target_size,
                                                          target_amax_error,
                                                          buf, buf_size,
                                                          actual_size,
                                                          actual_amax_error, false,
							  is_compressing_header_, is_delta_encoding_);
  }
  
  
  void quantized_array_codec::int32_decompress(int32_array2_t& array,
                                               const void* buf,
                                               std::size_t buf_size)  {
    quantized_array_codec_detail::int_array_decompress(bitio_,
                                                       array,
                                                       buf, buf_size,
						       is_compressing_header_, is_delta_encoding_);
  }
  
  void quantized_array_codec::uint32_compress(const uint32_array2_t& array,
                                              std::size_t target_size,
                                              double target_rms_error,
                                              void* buf,
                                              std::size_t buf_size,
                                              std::size_t* actual_size,
                                              double* actual_rms_error) {
    quantized_array_codec_detail::int_array_compress_rms(bitio_,
                                                         array,
                                                         target_size,
                                                         target_rms_error,
                                                         buf, buf_size,
                                                         actual_size,
                                                         actual_rms_error,
							 is_compressing_header_, is_delta_encoding_);
  }
  
  void quantized_array_codec::uint32_compress_amax(const uint32_array2_t& array,
                                                   std::size_t target_size,
                                                   double target_amax_error,
                                                   void* buf,
                                                   std::size_t buf_size,
                                                   std::size_t* actual_size,
                                                   double* actual_amax_error)  {
    quantized_array_codec_detail::int_array_compress_amax(bitio_,
                                                          array,
                                                          target_size,
                                                          target_amax_error,
                                                          buf, buf_size,
                                                          actual_size,
                                                          actual_amax_error, false,
							  is_compressing_header_, is_delta_encoding_);
  }
  
  
  void quantized_array_codec::uint32_decompress(uint32_array2_t& array,
                                                const void* buf,
                                                std::size_t buf_size)  {
    quantized_array_codec_detail::int_array_decompress(bitio_,
                                                       array,
                                                       buf, buf_size,
						       is_compressing_header_, is_delta_encoding_);
  }

}
