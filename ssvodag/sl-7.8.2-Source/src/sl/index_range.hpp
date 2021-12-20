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
#ifndef SL_INDEX_RANGE_HPP
#define SL_INDEX_RANGE_HPP

namespace sl {

  /// Ranges of integer indexes, and associated step-sizes for traversal 
  class index_range {
  protected:
    const int size_, start_, step_;
  public:

    /// Create the range [ start, end [, with a stepsize step
    inline index_range(int start, 
		       int end, 
		       int step=1) : 
      size_((end - start)/step), start_(start), step_(step) {
    }
    
    /// Init from other
    inline index_range(const index_range& other) : 
      size_(other.size_), start_(other.start_), step_(other.step_) {
    }

    /// The i-th index in the range
    inline int operator()(int i) const { 
      return start_ + step_*i ; 
    }

    /// The number of steps of this
    inline int size() const { 
      return size_; 
    }
    
    /// The start index of this
    inline int start() const { 
      return start_; 
    }
    
    /// The stepsize of this
    inline int step() const { 
      return step_; 
    }

    /// Translate the range by offset
    inline index_range operator+(int offset) const { 
      return index_range(start_+offset, start_+size_*step_+offset, step_); 
    }

    /// Translate the range by -offset
    inline index_range operator-(int offset) const { 
      return index_range(start_-offset, start_+size_*step_-offset, step_); 
    }

  };
}

#endif
