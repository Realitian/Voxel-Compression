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
#ifndef SL_TIME_POINT_HPP
#define SL_TIME_POINT_HPP

#include <sl/time_duration.hpp>

namespace sl {

  /**
   *  A point in time.
   */
  class time_point {
  protected: 
    uint64_t value;
    
  public: 
    
    /// Initialize to zero.
    inline time_point() { value = 0; }

    /// Initialize to the value passed as parameter (usecs from zero).
    inline explicit time_point(uint64_t usec) { value = usec; }
    

  public: // Serialization
    
    inline void store_to(output_serializer& s) const {
      s << value;
    }
    
    inline void retrieve_from(input_serializer& s) {
      s >> value;
    }

  public: // Comparison
	
    /// Is this equal to other?
    inline bool operator == (const time_point& other) const {
      return value == other.value;
    }

    /// Is this less than other?
    inline bool operator < (const time_point& other) const {
      return value < other.value;
    }
    
    SL_OP_COMPARABLE1(time_point);          // provides >, <=, >=
    SL_OP_EQUALITY_COMPARABLE1(time_point); // provides !=
	
    public: // Algebra
	
    /// Increment this by other.
    inline time_point& operator += (const time_duration& other) {
      value += other.as_microseconds();
      return *this;
    }
    
    /// Decrement this by other.
    inline time_point& operator -= (const time_duration& other) {
      value -= other.as_microseconds();
      return *this;
    }
	
    /// Subtract other from this.
    inline time_duration operator - (const time_point& other) {
      return time_duration(value - other.value);
    }
	
    SL_OP_ADDABLE2(time_point, time_duration);
    SL_OP_SUBTRACTABLE2(time_point, time_duration);
    
  };
    
};

// Arithmetic operations overload 

SL_OP_ADDABLE2_OVERLOADS(SL_OP_NO_TEMPLATE, sl::time_point, sl::time_duration);
SL_OP_SUBTRACTABLE2_OVERLOADS(SL_OP_NO_TEMPLATE, sl::time_point, sl::time_duration);

#endif
