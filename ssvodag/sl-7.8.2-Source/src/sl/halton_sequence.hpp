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
#ifndef SL_HALTON_SEQUENCE_HPP
#define SL_HALTON_SEQUENCE_HPP

#include <sl/utility.hpp>
#include <sl/fixed_size_point.hpp>

namespace sl {

  /**
   * Halton sequences are quasi random sequences used to generate points in space for numerical
   * methods such as Monte Carlo simulations. Although these sequences are deterministic
   * they are of low discrepancy, that is, appear to be random for many purposes and
   * do a better job of filling in space, avoiding the large gaps that can occur
   * with a random number generator.
   */
  template <size_t DIMENSION, class T>
  class halton_sequence {
  public: // Constants and types
    typedef halton_sequence<DIMENSION,T> self_t;

    enum { dimension = DIMENSION };
    typedef T       value_t;
    
    typedef fixed_size_point<DIMENSION,T>   point_t;
    typedef typename point_t::vector_t      vector_t;
    typedef typename point_t::dual_vector_t dual_vector_t;

  protected: // Storage
    uint64_t n_;
    uint64_t n0_;
    uint64_t base_[dimension];
    double   radical_[dimension];
    double   x_[dimension];
    point_t  p_;
    
  protected:

    SL_COMPILE_TIME_CHECK("Good dimension", DIMENSION > 0 && DIMENSION<256);
    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);

  protected:
    
    inline void update_p() {
      for (int j=0; j<dimension; ++j) {
	p_[j] = value_t(x_[j]);
      }
    }

  public:

    /// Initialize the generator
    halton_sequence(uint64_t seed = 0) {
      n_  = seed;
      n0_ = seed;

      for (int j=0; j<dimension; ++j) {
	base_[j] = small_prime((unsigned char)j);
	radical_[j] = 1.0/double(base_[j]);
	x_[j] = 0.0;
	p_[j] = value_t(0.0);
      }

      restart(n0_);
    }

    ~halton_sequence() {
    }

    /// Restart generation
    void restart(uint64_t seed = 0) {
      n_ = seed;
   
      for (int j = 0; j < dimension; ++j) {
	x_[j] = 0.0;
	double factor = radical_[j];
	uint64_t b = base_[j];
      
	for(uint64_t im = n_; im > 0; im /= b, factor *= radical_[j]) {
	  x_[j] += factor * (double) (im % b);
	}
      }
      
      update_p();
    }

    /// Return current quasi random point in the sequence
    const point_t& current() const {
      return p_;
    }

    /// Advance to next quasi random point in the sequence
    void forth() {
      ++n_;

      if (n_ & 8191) {
	// First 8192 numbers in the sequence
	const double almost_one = 1.0 - 1e-10;

	for (int j = 0; j < dimension; j++) {
	  double remainder = almost_one - x_[j];
	 
	  if (remainder < 0.0) {
	    x_[j] = 0.0;
	  } else if (radical_[j] < remainder) {
	    x_[j] += radical_[j];
	  } else {
	    double h = radical_[j];
	    double hh;
	    do {
	       hh = h;
	       h *= radical_[j];
	    } while(h >= remainder);
	    x_[j] += hh + h - 1.0;
	  }
	}
      } else if (n_ >= 1073741824) {  // == 2^30
	restart(0);
      } else {
	restart(n_);
      }

      update_p();
    }

    /// Advance to next quasi random point in the sequence
    void operator++() {
      forth();
    }
  }; // class halton_sequence

  typedef halton_sequence<1,float> halton_sequence1f;
  typedef halton_sequence<2,float> halton_sequence2f;
  typedef halton_sequence<3,float> halton_sequence3f;
  typedef halton_sequence<4,float> halton_sequence4f;
  typedef halton_sequence<5,float> halton_sequence5f;
  typedef halton_sequence<6,float> halton_sequence6f;
  
  typedef halton_sequence<1,double> halton_sequence1d;
  typedef halton_sequence<2,double> halton_sequence2d;
  typedef halton_sequence<3,double> halton_sequence3d;
  typedef halton_sequence<4,double> halton_sequence4d;
  typedef halton_sequence<5,double> halton_sequence5d;
  typedef halton_sequence<6,double> halton_sequence6d;
   
} // namespace sl

#endif

    
