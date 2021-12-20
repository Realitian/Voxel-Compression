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
#ifndef SL_INDEX_HPP
#define SL_INDEX_HPP

#include <sl/assert.hpp>
#include <sl/operators.hpp>
#include <iostream>
#include <iomanip>

namespace sl {
  
  /// Multidimensional index
  template <size_t G_size>
  class index {
  public:
    typedef index<G_size> this_t;

    enum { max_rank = (G_size - 1) };

  public:

    // avoid problems with rank 0 indexes...
    size_t i_[G_size == 0 ? 1 : G_size]; 

  public:

    /// Default construct (0)
    inline index() {
      for (size_t i=0; i<G_size; i++) {
	i_[i] = 0;
      }
    }

    /// Explicit construction (dim >= 1)
    inline index(size_t i0) {
      SL_REQUIRE("Good dimension", count() >= 1);
      i_[0] = i0;
      for (size_t i=1; i<G_size; i++) i_[i] = 0;
    }

    /// Explicit construction (dim >= 2)
    inline index(size_t i0, size_t i1) {
      SL_REQUIRE("Good dimension", count() >= 2);
      i_[0] = i0;
      i_[1] = i1;
      for (size_t i=2; i<G_size; i++) i_[i] = 0;
    }

    /// Explicit construction (dim >= 3)
    inline index(size_t i0, size_t i1, size_t i2) {
      SL_REQUIRE("Good dimension", count() >= 3);
      i_[0] = i0;
      i_[1] = i1;
      i_[2] = i2;
      for (size_t i=3; i<G_size; i++) i_[i] = 0;
    }

    /// Explicit construction (dim >= 4)
    inline index(size_t i0, size_t i1, size_t i2, size_t i3) {
      SL_REQUIRE("Good dimension", count() >= 4);
      i_[0] = i0;
      i_[1] = i1;
      i_[2] = i2;
      i_[3] = i3;
      for (size_t i=4; i<G_size; i++) i_[i] = 0;
    }
    
    // etc...

  public: // Index slicing

    template <size_t G_size2>
    inline void from_lower_dimensional_extent(const index<G_size2>& other) {
      SL_REQUIRE("Good size", G_size2 <= G_size);
      for (size_t r = 0; r<G_size2; ++r) {
	i_[r] = other[r];
      }
      for (size_t r = G_size2; r<G_size; ++r) {
	i_[r] = 1;
      }      
    }

    template <size_t G_size2>
    inline void from_lower_dimensional_index(const index<G_size2>& other) {
      SL_REQUIRE("Good size",G_size2 <= G_size);
      for (size_t r = 0; r<G_size2; ++r) {
	i_[r] = other[r];
      }
      for (size_t r = G_size2; r<G_size; ++r) {
	i_[r] = 0;
      }
    }

  public: // Element access

    /// The i-th index (alias of operator())
    inline size_t operator[](size_t i) const {
      SL_REQUIRE("Good index", i < count());
      return i_[i];
    }

    /// Write access to the i-th index (alias of operator())
    inline size_t& operator[](size_t i) {
      SL_REQUIRE("Good index", i < count());
      return i_[i];
    }

    /// The i-th index
    inline size_t operator()(size_t i) const {
      SL_REQUIRE("Good index", i < count());
      return i_[i];
    }

    /// Write access to the i-th index
    inline size_t& operator()(size_t i) {
      SL_REQUIRE("Good index", i < count());
      return i_[i];
    }

    /// The number of indices
    inline size_t count() const {
      return G_size;
    }

    /// The number of elements from zero to this
    inline size_t element_count() const {
      size_t result=1;
      for (size_t r=0; r<count(); ++r) {
	result *= i_[r];
      }
      return result;
    }

  public: // Looping and increment

    /**
     *  The value of an index when exiting from a
     *  loop that has this as end bound.
     */
    inline this_t end_mark() const {
      this_t result;
      result[0] = (*this)[0];
      return result;
    }

    /**
     * The extent if this is the last element
     */
    inline this_t extent_from_last() const {
      this_t result = *this;
      for (size_t r=0; r<count(); ++r) {
        result[r] += 1;
      }
      return result;
    }
    
    /**
     *  Increment all elements of this, simulating
     *  nested loops that iterate on all indexes. 
     *  for (i[N] = 0; i[N] < extent[N]; i[N]++) {
     *    for (i[N-1] = 0; i[N-1] < extent[N-1]; i[N]++) {
     *      ...
     *      for (i[0] = 0; i[0] < extent[0]; i[0]++) {
     *      }
     *    }
     *  }
     *  On exit, this is set to the next index in the
     *  loop, or to the end_mark() on outer loop exit.
     */
    inline void increment(const this_t& extent) {
      SL_REQUIRE("Within bounds", all_lt(extent));
      size_t r=max_rank;
      ++i_[r];
      while ((r>0) && (i_[r] >= extent[r])) {
	i_[r] = 0; --r; ++i_[r];
      }
      SL_ENSURE("Within bounds or at end", all_lt(extent) || *this == extent.end_mark());
    }

    inline void increment(const this_t& l,
                          const this_t& h) {
      SL_REQUIRE("Within bounds", all_lt(l+l.box_extent(h)));
      size_t r=max_rank;
      ++i_[r];
      while ((r>0) && (i_[r] > h[r])) {
	i_[r] = l[r]; --r; ++i_[r];
      }
    }

  public: // Operations

    inline this_t box_extent(const this_t& h) const {
      this_t result;
      for (std::size_t r=0; r<count(); ++r) {
        result[r] = h[r]-(*this)[r]+1;
      }
      return result;
    }
      
    /// Increment this by other
    inline this_t operator += (const this_t& other) {
      for (std::size_t r=0; r<count(); ++r) {
        i_[r] += other[r];
      }
      return *this;
    }

    /// Decrement this by other.
    inline this_t operator -= (const this_t& other) {
      for (std::size_t r=0; r<count(); ++r) {
        i_[r] -= other[r];
      }
      return *this;
    }

    SL_OP_ADDITIVE_GROUP(this_t);
      
  public: // Comparison
    
    /// Is this[i] < other[i] for all i?
    inline bool all_lt(const this_t& other) const {
      bool result = true;
      for (size_t r=0; r<count() && result; ++r) {
	result = result && (i_[r] < other[r]);
      }
      return result;
    }

    /// Is this[i] <= other[i] for all i?
    inline bool all_le(const this_t& other) const {
      bool result = true;
      for (size_t r=0; r<count() && result; ++r) {
	result = result && (i_[r] <= other[r]);
      }
      return result;
    }

    /// Is l[i] <= this[i] <= h[i] for all i?
    inline bool all_in(const this_t& l, const this_t& h) const {
      bool result = true;
      for (size_t r=0; r<count() && result; ++r) {
	result = result && (l[r] <= i_[r]) && (i_[r] <= h[r]);
      }
      return result;
    }

    /// Is this < other in a sequence where the fastest varying index is the highest in rank?
    inline bool operator<(const this_t& other) const {
      for (size_t r=0; r<count(); ++r) {
	if (i_[r] < other(r)) {
	  return true;
	} else if (i_[r] > other(r)) {
	  return false;
	}
      }
      return false;
    }

    /// Is this == other in a sequence where the fastest varying index is the highest in rank?
    inline bool operator==(const this_t& other) const {
      bool result = true;
      for (size_t r=0; r<count() && result; ++r) {
	result = result && (i_[r] == other(r));
      }
      return result;
    }
      
    /// Is this != other in a sequence where the fastest varying index is the highest in rank?
    inline bool operator!=(const this_t& other) const {
      return !(*this == other);
    }
    
    /// Is this > other in a sequence where the fastest varying index is the highest in rank?
    inline bool operator>(const this_t& other) const {
      return other < *this;
    }

    /// Is this >= other in a sequence where the fastest varying index is the highest in rank?
    inline bool operator>=(const this_t& other) const {
      return !(*this < other);
    }

    /// Is this <= other in a sequence where the fastest varying index is the highest in rank?
    inline bool operator<=(const this_t& other) const {
      return !(*this > other);
    }
  };
}

//---------------------------------------------------------------------------
// Comparison

template <size_t G_size>
inline int compare(const sl::index<G_size>& arg1,
		   const sl::index<G_size>& arg2) {
  for (size_t r=0; r<G_size; ++r) {
    if (arg1[r] < arg2[r]) return -1;
    if (arg1[r] > arg2[r]) return  1;
  }
  return 0;
}


//---------------------------------------------------------------------------
// I/O

/// Write rhs to stream s
template <size_t G_rank>
std::ostream& operator<<(std::ostream& s, 
			 const sl::index<G_rank>& rhs) {
  s << "[ ";
  for (size_t i=0; i<rhs.count(); i++) {
    s << std::setw(3) << rhs(i) << " ";
  }
  s << "] ";
  return s;
}

//---------------------------------------------------------------------------
// typedef

namespace sl {
  typedef index<1> index1_t;
  typedef index<2> index2_t;
  typedef index<3> index3_t;
  typedef index<4> index4_t;
}

 
#endif
