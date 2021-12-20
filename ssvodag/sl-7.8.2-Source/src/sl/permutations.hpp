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
#ifndef SL_PERMUTATIONS_HPP
#define SL_PERMUTATIONS_HPP

#include <sl/config.hpp>
#include <cassert>

namespace sl {

  /**
   *  Enumerate all possible permutations of N elements by generating a sequence
   *  in which each permutation differs from the previous one by only a single swap
   *  of an element with the previous one. Algorithm due to Gideon Ehrlich:
   *
   *  Alg E. 7.2.1.2 in Donald E. Knuth: The Art of Computer Programming,
   *  pre-fascicles for Volume 4. URL: http://www-cs-staff.stanford.edu/~knuth/
   */
  class star_swap_permutator {
  public:
    typedef star_swap_permutator this_t;
  protected:
    std::size_t n_;     // number of elements
    std::size_t swp_;   // index of element swapped with index 0
    std::size_t b_[32]; // auxiliary array
    std::size_t c_[32]; // auxiliary array: mixed radix number in rising factorial base
  public:
    explicit inline star_swap_permutator(std::size_t n) {
      assert(n<32);
      n_ = n;
      restart();
    }

    inline ~star_swap_permutator() {
    }

    inline void restart() {
      swp_ = 0;
      for (std::size_t k=0; k<n_; ++k)  b_[k] = k;
      for (std::size_t k=0; k<=n_; ++k)  c_[k] = 0;
    }

    inline bool off() const {
      return swp_ >= n_;
    }
    
    inline this_t& operator++() {
      assert(!off());
      std::size_t k = 1;
      while ( c_[k]==k )  { c_[k]=0;  ++k; }

      if ( k == n_ ) {
	// Off
	swp_ = n_;
      } else {
	// Not off
        ++c_[k];
        swp_ = b_[k];

        std::size_t j = 1;
        --k;
        while (j < k) {  // Knuth: < 0.18 iterations/call
	  std::swap(b_[j], b_[k]);
	  ++j; --k;
        }
      }
      return *this;
    }
    
    /// Pre-increment - this is much slower than post-increment, don't use it if not necessary
    inline this_t operator++(int) { 
      this_t tmp(*this); ++(*this); return tmp; 
    }

    // The index that must be swapped with index 0
    inline std::size_t value() const {
      assert(!off());
      return swp_;
    }
    
  }; // class star_swap_permutator

} // namespace sl

#endif
