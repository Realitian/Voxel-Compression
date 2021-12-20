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
#ifndef SL_RANDOM_HPP
#define SL_RANDOM_HPP

#include <sl/cstdint.hpp> // for int32_t
#include <sl/utility.hpp> // for tags
#include <sl/math.hpp>    // for numeric_limits
#include <cassert>

namespace sl {

  /// Random number generators
  namespace random {

    /**
     *  A 32 bit integer pseudo random number uniformly distributed between 0 and 0xffffffff
     *  Implementation based on George Marsaglia method.
     */
    class irng_marsaglia {
    protected:
      volatile uint32_t seed_x_;
      volatile uint32_t seed_y_;
      volatile uint32_t seed_z_;
      volatile uint32_t seed_c_;

      enum {
	reference_seed_x = 123456789,
	reference_seed_y = 987654321,
	reference_seed_z = 43219876,
	reference_seed_c = 6543217
      };

    public:

      inline void set_seed(uint32_t x = reference_seed_x,
			   uint32_t y = reference_seed_y,
			   uint32_t z = reference_seed_z,
			   uint32_t c = reference_seed_c) {
	if (x == 0) x = reference_seed_x;
	if (y == 0) y = reference_seed_y;
	if (z == 0) z = reference_seed_z;
	if (c == 0) c = reference_seed_c;
	seed_x_ = x;
	seed_y_ = y;
	seed_z_ = z;
	seed_c_ = c;
      }


      inline irng_marsaglia() {
	set_seed();
      }

      /// Return a random number uniformely distributed between 0 and 0xffffffff
      inline uint32_t value() {
	// Generator JKISS from Good Practice in (Pseudo) Random Number Generation
	// for Bioinformatics Applications David Jones, UCL Bioinformatics Group
	// www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf
	uint64_t t;
	seed_x_ = 314527869 * seed_x_ + 1234567;
	seed_y_ ^= seed_y_ << 5; seed_y_ ^= seed_y_ >> 7; seed_y_ ^= seed_y_ << 22;
	t = 4294584393ULL * seed_z_ + seed_c_; seed_c_ = t >> 32; seed_z_ = t;
	return seed_x_ + seed_y_ + seed_z_;
      }

      /// Return a random number uniformely distributed between 0 and n included
      inline uint32_t value_leq(uint32_t n) {
	const uint32_t rand_max =  uint32_t(0xffffffff);
	
	uint32_t result;

	if (n==0) {
	  // No randomness
	  result = 0;
	} else if (n==rand_max) {
	  // Full range
	  result = value();
	} else {
#if 0
	  // Fast modulo
	  result = rnd % (n+uint32_t(1));
#else
	  // This implementation avoids modulo bias
	  uint32_t limit = rand_max - rand_max % (n+uint32_t(1));
	  uint32_t rnd = value();
	  while (rnd >= limit) {
	    rnd = value();
	  }
	  result = rnd % (n+uint32_t(1));
#endif
	}

	return result;
      }

      /// Return a random number uniformely distributed between lo and hi included
      inline uint32_t value_within(uint32_t lo, uint32_t hi) {
	assert(lo<=hi);

	return lo + value_leq(hi-lo);
      }

      /**
       *  Chose a random subset of k elements from {0,1,...n-1}, and return
       *  these indices in the k first elements of the a[] array.
       */
      void pick_k_out_of_n_in(uint32_t k, uint32_t *a,
			      uint32_t N);

    }; 

    /// The default random number generator
    typedef irng_marsaglia std_irng_t;

    /**
     *  The default sharing state
     *  WARNING: UP TO sl 7.7.0 THE DEFAULT WAS "SHARED", WHICH MEANS
     *  THAT ALL RANDOM GENERATORS WERE DEPENDENT, AND, FOR INSTANCE,
     *  SETTING THE SEED TO ONE OF THEM CHANGED ALL SEEDS. NOW THE
     *  DEFAULT IS NOT SHARING, WHICH MAKES IT POSSIBLE TO USE THEM
     *  MORE EASILY IN PARALLEL PROGRAMS
     */
    typedef ::sl::tags::not_shared_state std_state_sharing_t;

    namespace detail {
      /**
       *  A wrapper for integer random number generators. If
       *  T_STATE is tags::shared, the integer random number
       *  generator is shared among all instances, otherwise
       *  each instance gets its own.
       */
      template <class T_IRNG, class T_STATE>
      class irng_wrapper {
      };
      
      template <class T_IRNG>
      class irng_wrapper<T_IRNG, ::sl::tags::shared_state> {
      protected:
	static T_IRNG irng_;
      public:
	inline void set_seed(uint32_t x) {
	  (this->irng_).set_seed(x);
	}

      };

      template <class T_IRNG>
      T_IRNG irng_wrapper<T_IRNG, ::sl::tags::shared_state>::irng_;
      
      template <class T_IRNG>
      class irng_wrapper<T_IRNG, ::sl::tags::not_shared_state> {
      protected:
	T_IRNG irng_;
      public:
	inline void set_seed(uint32_t x) {
	  (this->irng_).set_seed(x);
	}
      };
    } // namespace detail

    /**
     *  Uniform random numbers in [0..1)
     */
    template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
    class uniform_closed_open: public detail::irng_wrapper<T_IRNG, T_STATE> {
    public: 
      typedef T value_t;
    public:
      value_t value() {
#       define detail_norm32open value_t(.2328306436538696289062500000000000000000E-09L)
#       define detail_norm64open value_t(.5421010862427522170037264004349708557129E-19L)
#       define detail_norm96open value_t(.1262177448353618888658765704452457967477E-28L)
#       define detail_norm128open value_t(.2938735877055718769921841343055614194547E-38L)

	value_t result =  (this->irng_).value()* detail_norm32open;	 
	if (std::numeric_limits<value_t>::digits > 32) {

	  result += (this->irng_).value() * detail_norm64open;
	} 
	if (std::numeric_limits<value_t>::digits > 64) {

	  result += (this->irng_).value() * detail_norm96open;
	}
	if (std::numeric_limits<value_t>::digits > 128) {
	  result += (this->irng_).value() * detail_norm128open;
	}
	return result;
#       undef detail_norm32open
#       undef detail_norm64open
#       undef detail_norm96open
#       undef detail_norm128open
      }
    };

    /**
     *  Uniform value numbers in [0..1]
     */
    template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
    class uniform_closed: public detail::irng_wrapper<T_IRNG, T_STATE> {
    public: 
      typedef T value_t;

    public:
      value_t value() {
#       define detail_norm32closed   value_t(.2328306437080797375431469961868475648078E-9L)
#       define detail_norm64closed1  value_t(.23283064365386962891887177448353618888727188481031E-9L)
#       define detail_norm64closed2  value_t(.54210108624275221703311375920552804341370213034169E-19L)
#       define detail_norm96closed1  value_t(.2328306436538696289062500000029387358771E-9L)
#       define detail_norm96closed2  value_t(.5421010862427522170037264004418131333707E-19L)
#       define detail_norm96closed3  value_t(.1262177448353618888658765704468388886588E-28L)
#       define detail_norm128closed1 value_t(.2328306436538696289062500000000000000007E-9L)
#       define detail_norm128closed2 value_t(.5421010862427522170037264004349708557145E-19L)
#       define detail_norm128closed3 value_t(.1262177448353618888658765704452457967481E-28L)
#       define detail_norm128closed4 value_t(.2938735877055718769921841343055614194555E-38L)

	value_t result;
	if (std::numeric_limits<value_t>::digits > 96) {
	  uint32_t i1 = (this->irng_).value();
	  uint32_t i2 = (this->irng_).value();
	  uint32_t i3 = (this->irng_).value();
	  uint32_t i4 = (this->irng_).value();

	  result = i1 * detail_norm128closed1 + i2 * detail_norm128closed2
	    + i3 * detail_norm128closed3 + i4 * detail_norm128closed4;
	} else if (std::numeric_limits<value_t>::digits > 64) {
	  uint32_t i1 = (this->irng_).value();
	  uint32_t i2 = (this->irng_).value();
	  uint32_t i3 = (this->irng_).value();

	  result = i1 * detail_norm96closed1 + i2 * detail_norm96closed2
	    + i3 * detail_norm96closed3;
	} else if (std::numeric_limits<value_t>::digits > 32) {
	  uint32_t i1 = (this->irng_).value();
	  uint32_t i2 = (this->irng_).value();

	  result = i1 * detail_norm64closed1 + i2 * detail_norm64closed2;
	} else {
	  uint32_t i = (this->irng_).value();
	  result = i * detail_norm32closed;
	}
	return result;
#       undef detail_norm32closed   
#       undef detail_norm64closed1  
#       undef detail_norm64closed2  
#       undef detail_norm96closed1  
#       undef detail_norm96closed2  
#       undef detail_norm96closed3  
#       undef detail_norm128closed1 
#       undef detail_norm128closed2 
#       undef detail_norm128closed3 
#       undef detail_norm128closed4 
      }
    };

    /**
     *  Uniform random numbers in [0..1]
     */
    template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
    class uniform_open: public uniform_closed_open<T,T_IRNG, T_STATE> {
    public: 
      typedef T value_t;
      typedef uniform_closed_open<T,T_IRNG, T_STATE> super_t;
    public:
      inline value_t value() {
	value_t result;
	do {
	  result = super_t::value();
	} while (::sl::is_zero(result));
	return result;
      }
    };

    /**
     *  Uniform random numbers in (0..1]
     */
    template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
    class uniform_open_closed: public uniform_closed_open<T,T_IRNG, T_STATE> {
    public: 
      typedef T value_t;
      typedef uniform_closed_open<T,T_IRNG, T_STATE> super_t;
    public:
      inline value_t value() {
	return ::sl::scalar_math<value_t>::one() - super_t::value();
      }
    };


    /**
     *  Uniform random numbers in [0..1) (Alias of uniform_closed_open)
     */
    template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
    class uniform: public uniform_closed_open<T, T_IRNG, T_STATE> {
    public: 
      typedef T value_t;
    };

    /**
     *  Uniform integer random number generator in 0..n
     */
    template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
    class discrete_uniform: public detail::irng_wrapper<T_IRNG, T_STATE> {
    public: 
      typedef T value_t;
    protected:
      uint32_t n_;
    public:

      inline discrete_uniform(T n) {
	SL_REQUIRE("Good n", n > 0 && n<=std::numeric_limits<uint32_t>::max());
	n_ = uint32_t(n);
      }

      inline value_t value() {
	return value_t((this->irng_).value() % n_);
      }
    };

    /**
     *  Normal random number generator. Based on:
     *   J.L. Leva, Algorithm 712. A normal random number generator, ACM Trans. 
     *   Math. Softw.  18 (1992) 454--455. 
     *
     *   http://www.acm.org/pubs/citations/journals/toms/1992-18-4/p449-leva/
     */
     template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
     class unit_normal: public uniform_open<T, T_IRNG, T_STATE> {
     public: 
       typedef T value_t;
       typedef uniform_open<T, T_IRNG, T_STATE> super_t;
     public:
       value_t value() {
	 const value_t s = value_t(0.449871), t = value_t(-0.386595), a = value_t(0.19600), b = value_t(0.25472);
	 const value_t r1 = value_t(0.27597), r2 = value_t(0.27846);

	 value_t u, v;

	 for (;;) {
	   // Generate P = (u,v) uniform in rectangle enclosing
	   // acceptance region:
	   //   0 < u < 1
	   // - sqrt(2/e) < v < sqrt(2/e)
	   // The constant below is 2*sqrt(2/e).
	   
	   u = super_t::value();
	   v = value_t(1.715527769921413592960379282557544956242L) 
	     * (super_t::value() - value_t(0.5));
	   
	   // Evaluate the quadratic form
	   value_t x = u - s;
	   value_t y = sl::abs(v) - t;
	   value_t q = x*x + y*(a*y - b*x);
	   
	   // Accept P if inside inner ellipse
	   if (q < r1)
	     break;
	   
	   // Reject P if outside outer ellipse
	   if (q > r2)
	     continue;
	   
	   // Between ellipses: perform exact test
	   if (v*v <= -4.0 * std::log(u)*u*u)
	     break;
	 }
	 
	 return v/u;
       }
	 
     }; // class normal
    
    
    /**
     *  Normal random number generator.
     */
     template <class T, class T_IRNG = std_irng_t, class T_STATE = std_state_sharing_t>
     class normal: public unit_normal<T, T_IRNG, T_STATE> {
     public: 
       typedef T value_t;
       typedef unit_normal<T, T_IRNG, T_STATE> super_t;
     protected:
       T mean_;
       T std_dev_;
     public:

       inline normal(T mean, T std_dev): mean_(mean), std_dev_(std_dev) {
	 SL_REQUIRE("Good std_dev", sl::is_non_negative(std_dev));
       }

       inline value_t random() {
	 return mean_ + std_dev_ * super_t::value();
       }

     };

  } // namespace random

} // namespace sl

#endif
