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
#include <sl/random.hpp>

namespace sl {
  namespace random {
    namespace detail {
      template <>
      std_irng_t irng_wrapper<std_irng_t, ::sl::tags::shared_state>::irng_;
    }

    
    void irng_marsaglia::pick_k_out_of_n_in(uint32_t K,
					    uint32_t *a,
					    uint32_t N) {
      assert(K <= N);

      if (K>0) {
	uint32_t n = N;
	uint32_t k = K;
	
	// while k is close to n, go though all candidates from n-1 to 0, and pick
	// each one with probability k/n
	while((n > k) && (k > n/2)) {
	  // each number has probability k/n of being picked up 
	  if (k > value_leq(n-1)) { 
	    // pick this one up 
	    k--;
	    n--;
	    a[k]= n;
	  } else {
	    // don't pick this one
	    n--;
	  }
	}
	
	if (n == k) {
	  // we've got k numbers to sample from a set of k, easy...
	  for (uint32_t i=0; i<n; ++i) {
	    a[i] = i;
	  }
	  k = 0; n=0;
	}
	if (k > 0) {
	  assert(k <= n/2);

	  // Use ranksb method
	  // pp 38-9, Chapter 4 of Nijenhuis & Wilk (1975) Combinatorial Algorithms, Academic Press.
	  
	  if (k == 1) {
	    a[0] = value_leq(n-1);
	    return;
	  } else {
	    /* Partition [0 : n-1] into k intervals:
	       Store the least element of the I'th interval in a[i] */
	    for (uint32_t i = 0; i < k; ++i) {
	      a[i] = i * n / k;
	    }
	    
	    /* Using a uniformly distributed random variable in the
	       range 0 <= x < n, make k selections of x such that
	       x lands in the remaining portion of some interval.
	       At each successful selection, reduce the remaining
	       portion of the interval by 1. */
	    for(uint32_t i = 0; i < k;  ++i) {
	      uint32_t x, l;
	      do {
		x = value_leq(n-1) + 1;
		l = (x * k - 1) / n;
	      } while (x <= a[l]);
	      ++a[l];
	    }

	    /* Collect the least elements of any interval which
	       absorbed a selection in the previous step into the
	       low-order indices of a. */
	    int32_t p = -1; // FIXME signed!
	    for (uint32_t i = 0; i < k; ++i) {
	      uint32_t m = a[i];
	      a[i] = 0;
	      if (m != i * n / k) {
		/* A non-empty partition */
		++p;
		a[p] = m;
	      }
	    }
	    /* Allocate space for each non-empty partition starting
	       from the high-order indices.  At the last position
	       in each partition's segment, store the interval index
	       of the partitions's successor. */
	    int32_t s = k-1;
	    for(; p >=0; p--) {
	      int32_t l = (a[p] * int32_t(k) - 1) / n;
	      int32_t ds = a[p] - l * n / k;
	      a[p] = 0;
	      a[s] = l + 1;
	      s -= ds;
	    }
	    
	    uint32_t m0 = 0;
	    uint32_t m  = k-1;
	    int32_t r = k-1;
	    for(int32_t l = k-1; l >= 0; l--) {
	      /* ranksb each of the sub-problems */
	      uint32_t x = a[l];
	      if (x != 0) {
		/* Start a new bin */
		r = l;
		m0 = (x - 1) * n / k;
		m = x * n / k - m0;
		/* The order of arithmetic operations is important!
		   The same rounding errors must be produced in each
		   computation of a boundary. */
	      }

	      /* m0 is the least element of the current (l'th)
		 interval.  m is the count of the number of
		 unselected members of the interval. */
	      x = m0 + value_leq(uint32_t(m-1));

	      /* Bubble Merge the (x-base_l)'th unselected member
		 of the current interval into the current interval's
		 segment (a [l..r]). */

	      int32_t i = l;
	      while (i < r && x >= a[i+1]) {
		a[i] = a[i+1];
		++x;
		++i;
	      }
	      a[i] = x;
	      --m;
	    } // for each l
	  } // if k>1
	} // if k>0
      } // if K> 0
    }
    
  }
}
