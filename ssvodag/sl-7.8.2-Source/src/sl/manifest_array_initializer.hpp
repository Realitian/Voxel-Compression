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
#ifndef SL_MANIFEST_ARRAY_INITIALIZER_HPP
#define SL_MANIFEST_ARRAY_INITIALIZER_HPP

#include <sl/assert.hpp>

namespace sl {

  /// Utility class for initializing 1D arrays
  template <class T_ARRAY, class T_VALUE>
  class manifest_array1d_initializer {
  protected:
    T_ARRAY a_;
    size_t  nmax_;
    size_t  i_;
  public:

    inline manifest_array1d_initializer(T_ARRAY a,  
					size_t  nmax,
					T_VALUE v):
      a_(a), nmax_(nmax), i_(0) {
      SL_REQUIRE("Good size", nmax > 0);
      a_[i_] = v;
    };

    inline ~manifest_array1d_initializer() {
      SL_REQUIRE("Matching sizes", i_ == 0 || i_ == nmax_-1);
      if (i_ == 0) { 
	// Fill with constant
	T_VALUE v = a_[0];
	for (size_t i=1; i<nmax_; i++) {
	  a_[i] = v;
	}
      }
    }
  
    inline manifest_array1d_initializer<T_ARRAY, T_VALUE>& operator, (T_VALUE v) {
      SL_REQUIRE("Good index", i_ < nmax_-1);
      i_++;
      a_[i_] = v;
      return *this;
    }

  }; 

  /// Utility class for initializing 2D arrays
  template <class T_ARRAY, class T_VALUE>
  class manifest_array2d_initializer {
  protected:
    T_ARRAY a_;
    size_t  row_count_;
    size_t  col_count_;
    size_t  i_;
    size_t  j_;
  public:

    inline manifest_array2d_initializer(T_ARRAY a,  
					size_t  row_count,
					size_t  col_count,
					T_VALUE v):
      a_(a), row_count_(row_count), col_count_(col_count), i_(0), j_(0) {
      SL_REQUIRE("Good size", row_count > 0 && col_count > 0);
      a_(i_,j_) = v;
    };

    inline ~manifest_array2d_initializer() {
      SL_REQUIRE("Matching sizes", (i_ == 0 && j_ == 0) || (i_ == row_count_-1 && j_ == col_count_-1) );
      if (i_ == 0 && j_ == 0) { 
	// Fill with constant
	T_VALUE v = a_(0,0);
	for (size_t i=0; i<row_count_; i++) {
	  for (size_t j=0; j<col_count_; j++) {
	    a_(i,j) = v;
	  }
	}
      }
    }
  
    inline manifest_array2d_initializer<T_ARRAY, T_VALUE>& operator, (T_VALUE v) {
      SL_REQUIRE("Good index", i_ < row_count_-1 || j_ < col_count_-1);
      j_++; if (j_ == col_count_) { j_ = 0; i_++; }
      a_(i_,j_) = v;
      return *this;
    }

  }; 

} // namespace sl 

#endif

