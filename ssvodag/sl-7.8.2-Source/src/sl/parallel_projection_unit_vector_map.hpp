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
#ifndef SL_PARALLEL_PROJECTION_UNIT_VECTOR_MAP_HPP
#define SL_PARALLEL_PROJECTION_UNIT_VECTOR_MAP_HPP

#include <sl/fixed_size_vector.hpp>
#include <sl/math.hpp>

namespace sl {

  class parallel_projection_unit_vector_map {

  public:
    sl::uint32_t n_bits_;
    sl::uint32_t i_max_;
    float        eps_;

  protected: // Uniformely quantized numbers

    static inline float sigma(float r) {
      return r>=0.0f ? 1.0f : -1.0f;
    }

    inline sl::uint32_t packed(float r) const {
      return sl::median(0, int((r*0.5f+0.5f) * i_max_ + 0.5f), int(i_max_));
    }

    inline float unpacked(sl::uint32_t i) const {
      return sl::median(-1.0f, 1.0f,
			(float(i) / float(i_max_)) * 2.0f - 1.0f);
    }
    
    inline float quantized(float r) const {
      return unpacked(packed(r));
    }
    
  public: // Construction

    parallel_projection_unit_vector_map(sl::uint32_t n_bits = 8) {
      set_bit_count_per_component(n_bits);
    }

    ~parallel_projection_unit_vector_map() {
    }
    
    void set_bit_count_per_component(sl::uint32_t n_bits) {
      n_bits_ = n_bits;
      i_max_ =  (1<<(n_bits-1))+((1<<(n_bits-1))-1);
      eps_ = (1.0f-(-1.0f))/float(i_max_);
#if 0
      std::cerr << "nb: " << n_bits_ << " imax: " << i_max_ << " eps: " << eps_ << std::endl;
      std::cerr << 1.0f << " -> " << packed(1.0f) << ": " << quantized(1.0f) << std::endl;
      std::cerr << 0.0f << " -> " << packed(0.0f) << ": " << quantized(0.0f) << std::endl;
      std::cerr << 0.5f << " -> " << packed(0.5f) << ": " << quantized(0.5f) << std::endl;
      std::cerr << -1.0f << " -> " << packed(-1.0f) << ": " << quantized(-1.0f) << std::endl;
      std::cerr << -0.5f << " -> " << packed(-0.5f) << ": " << quantized(-0.5f) << std::endl;
#endif
    }

  public: // Mapping
    
    inline void map_in(float& u, float& v, bool& signz, const sl::vector3f& n) const {
      u = n[0];
      v = n[1];
      signz = (n[2]>=0);
    }
    
    inline void map_in(float& u, float& v, bool& signz, const sl::row_vector3f& n) const {
      map_in(u,v,signz, as_dual(n));
    }
  
    inline void unmap_in(sl::vector3f& n, float u, float v, bool signz) const {
      n[0]=u;
      n[1]=v;
      n[2]=std::sqrt(max(0.0f, 1.0f-u*u-v*v));
      if (!signz) n[2]*=-1.0f;
    }
    
    inline void unmap_in(sl::row_vector3f& n, float u, float v, bool signz) const {
      sl::vector3f dn; unmap_in(dn,u,v,signz); n=as_dual(dn);
    }

  public: // Packing


    // domain quantization
    inline void domain_map_in(sl::uint32_t& iu, sl::uint32_t& iv, bool& signz, const sl::vector3f& n) const {
      float u,v;
      map_in(u,v,signz, n);
      iu=packed(u); iv=packed(v);
    }
    
    inline void domain_map_in(sl::uint32_t& iu, sl::uint32_t& iv, bool& signz, const sl::row_vector3f& n) const {
      domain_map_in(iu,iv,signz, as_dual(n));
    }

    // range quantization
    inline void range_map_in(sl::uint32_t& iu, sl::uint32_t& iv, bool& signz, const sl::vector3f& n) const {
      float u0,v0;
      map_in(u0,v0,signz,n);
      sl::uint32_t iu0=packed(u0);
      sl::uint32_t iv0=packed(v0);
      sl::uint32_t du=(u0>=0.0f)?1:-1;
      sl::uint32_t dv=(v0>=0.0f)?1:-1;
      if (iu0==0 || iu0==i_max_) du=0;
      if (iv0==0 || iv0==i_max_) dv=0;
      
      sl::vector3f n_prime;
      unmap_in(n_prime, iu0,iv0,signz);
      iu=iu0; iv=iv0;
      float eps = n.dot(n_prime);
      
      for (std::size_t i=0; i<2; ++i) {
	sl::uint32_t iu1=iu0+i*du;
	for (std::size_t j=0; j<2; ++j) {
	  sl::uint32_t iv1=iv0+j*dv;
	  if (iu1!=iu0 || iv1!=iv0) {
	    unmap_in(n_prime, iu1,iv1,signz);
	    float eps1 = n.dot(n_prime);
	    if (eps1>eps) {
	      //std::cerr<<"XXX(" << iu1 << " " << iv1 << "): " << i << " " << j << std::endl;
	      eps=eps1;
	      iu=iu1; iv=iv1;
	    }
	  }
	}
      }
    }
    
    inline void range_map_in(sl::uint32_t& iu, sl::uint32_t& iv, bool& signz, const sl::row_vector3f& n) const {
      range_map_in(iu,iv,signz,as_dual(n));
    }

    inline void map_in(sl::uint32_t& iu, sl::uint32_t& iv, bool& signz, const sl::vector3f& n) const {
      range_map_in(iu,iv,signz,n);
    }

    inline void map_in(sl::uint32_t& iu, sl::uint32_t& iv, bool& signz, const sl::row_vector3f& n) const {
      range_map_in(iu,iv,signz,n);
    }
    
    inline void unmap_in(sl::vector3f& n, sl::uint32_t iu, sl::uint32_t iv, bool signz) const {
      float u = unpacked(iu);
      float v = unpacked(iv);
      unmap_in(n, u,v,signz);
    }

    inline void unmap_in(sl::row_vector3f& n, sl::uint32_t iu, sl::uint32_t iv, bool signz) const {
      float u = unpacked(iu);
      float v = unpacked(iv);
      unmap_in(n, u,v,signz);
    }

  public: // Accuracy
    
    /// An estimate of quantized representation accuracy (in radians)
    float accuracy() const {
      double dDotMin = 1.0;
      sl::vector3f n_worst;
      
      const std::size_t N=2048;
      for (std::size_t i=0; i<N; ++i) {
	float z = -1.0f + 2.0f*float(i)/float(N-1);
	for (std::size_t j=0; j<N; ++j) {
	  float theta = 0.0f + 2.0f*sl::scalar_math<float>::Pi()*float(i)/float(N-1);

	  sl::vector3f n=sl::vector3f(std::sqrt(1.0f-z*z)*std::cos(theta),
				      std::sqrt(1.0f-z*z)*std::sin(theta),
				      z).ok_normalized();
	  sl::uint32_t iu, iv; bool signz;
	  map_in(iu, iv, signz, n);
	  sl::vector3f n_prime;
	  unmap_in(n_prime, iu, iv, signz);
		  
	  double dDot =
	    double(n[0])*double(n_prime[0]) + 
	    double(n[1])*double(n_prime[1]) + 
	    double(n[2])*double(n_prime[2]);
	  if ( dDot < dDotMin ) {
	    dDotMin = dDot;

	    n_worst = n;
	  }
	}
      }

#if 0
      {
	sl::uint32_t iu, iv; bool signz;
	range_map_in(iu,iv, signz, n_worst);
	sl::vector3f n_prime;
	unmap_in(n_prime, iu, iv, signz);
	
	double dDot =
	  double(n_worst[0])*double(n_prime[0]) + 
	  double(n_worst[1])*double(n_prime[1]) + 
	  double(n_worst[2])*double(n_prime[2]);
	
	std::cerr << "WORST: " << 
	  n_worst[0] << " " << n_worst[1] << " " << n_worst[2] << "->" << iu << "," << iv << " " << signz << "=>" <<
	  n_prime[0] << " " << n_prime[1] << " " << n_prime[2] << "=> E=" 
		     <<
	  std::acos(dDotMin)*180.0/sl::scalar_math<double>::Pi() <<
	  " DOT: " << dDotMin << std::endl;
      }
#endif		    

      return std::acos(dDotMin);
    }
     
  };
}

#endif
