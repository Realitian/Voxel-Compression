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
#include <sl/cstdint.hpp>
#include <algorithm> // std::swap

#ifndef SL_PLHAAR_HPP
#define SL_PLHAAR_HPP

namespace sl {

  /**
   *  An Improved N-Bit to N-Bit Reversible Haar-Like Transform
   *  Joshua G. Senecal, Peter Lindstrom, Mark A. Duchaineau, Kenneth I. Joy
   *  Proc. Pacific Graphics 2004
   *
   *  If the inputs are N-bit values from a domain [0, 2^n-1] the bias parameter 
   *  should be set to 2^(n-1) to keep the high- and low-pass coefficient range equal 
   *  to the domain.
   */
  template <class T, intmax_t G_bias>
  class plhaar {
  public:
    typedef T value_t;
    
    static inline void apply(value_t *l, // low-pass output
			     value_t *h, // high-pass output
			     value_t a, // input #1
			     value_t b // input #2
			     ) {
      const value_t bias = value_t(G_bias);
      const value_t s = (a < bias), t = (b < bias);
      a += s; b += t; // (**) nudge origin
      if (s == t) { // A * B > 0?
	a -= b - bias; // H = A - B
	if ((a < bias) == s) { // |A| > |B|?
	  b += a - bias; // L = A (replaces L = B)
	}
      } else { // A * B < 0
	b += a - bias; // L = A + B
	if ((b < bias) == t) { // |B| > |A|?
	  a -= b - bias; // H = -B (replaces H = A)
	}
      }
      a -= s; b -= t; // (**) restore origin
      *l = b; *h = a; // store result
    }
    
    static void split(value_t* s, std::size_t N) {
      const std::size_t N2 = ((N/2)*2 == N) ? N : N+1;
      const std::size_t half = N2>>1;
      value_t* s_even = new value_t[half+1];
      value_t* s_odd  = new value_t[half+1];
      
      std::size_t n_half = 0;
      for (std::size_t i=0; i<N; i+=2) {
	s_even[n_half] = s[i]; 
	if (i+1<N) { 
	  s_odd[n_half] = s[i+1]; 
	}
	++n_half;
      }
      for (std::size_t i=0; i<half; ++i) {
	s[i] = s_even[i];
	if (i+half<N) { 
	  s[i+half] = s_odd[i]; 
	}
      }
      
      delete[] s_even;
      delete[] s_odd;
    }
      
    static void merge(value_t* s, std::size_t N) {
      const std::size_t N2 = ((N/2)*2 == N) ? N : N+1;
      const std::size_t half = N2>>1;
      value_t* s_even = new value_t[half+1];
      value_t* s_odd  = new value_t[half+1];
      
      for (std::size_t i=0; i<half; ++i) {
	s_even[i] = s[i];
	if (i+half<N) { 
	  s_odd[i] = s[i+half]; 
	}
      }
	
      std::size_t n_half = 0;
      for (std::size_t i=0; i<N; i+=2) {
	s[i] = s_even[n_half]; 
	if (i+1<N) { 
	  s[i+1] = s_odd[n_half]; 
	}
	++n_half;
      }
      
      delete[] s_even;
      delete[] s_odd;
    }
  
    static void apply(value_t* s, std::size_t N) {
      const std::size_t half = N>>1;
      for (std::size_t i = 0; i < half; ++i) {
	apply(&(s[i]), &(s[i+half]),
	      s[i], s[i+half]);
      }
    }

    static void forward1d(value_t* s, std::size_t N) {
      std::size_t Nk = N;
      while ((Nk>>1) > 0) {
	split(s,Nk);
	apply(s, Nk);
	Nk>>=1;
      }
    }
    
    static void backward1d(value_t* s, std::size_t N) {
      std::size_t RNk[64];
      std::size_t RNkN=0;
      std::size_t Nk = N;
      while ((Nk>>1) > 0) {
	RNk[RNkN]=Nk;
	++RNkN;
	Nk>>=1;
      }
      for (std::size_t k=0; k<RNkN; ++k) {
	apply(s,RNk[RNkN-k-1]);
	merge(s,RNk[RNkN-k-1]);
      }
    }

  }; // class plhaar

  typedef plhaar<uint8_t,  1<<7>   plhaar_uint8_t;
  typedef plhaar<uint16_t, 1<<15>  plhaar_uint16_t;
  typedef plhaar<uint32_t, 1<<31>  plhaar_uint32_t;
  typedef plhaar<int8_t, 0>        plhaar_int8_t;
  typedef plhaar<int16_t, 0>       plhaar_int16_t;
  typedef plhaar<int32_t, 0>       plhaar_int32_t;
  
} // namespace sl

#endif
