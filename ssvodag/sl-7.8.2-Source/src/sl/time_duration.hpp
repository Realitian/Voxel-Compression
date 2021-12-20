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
#ifndef SL_TIME_DURATION_HPP
#define SL_TIME_DURATION_HPP

#include <sl/operators.hpp>
#include <sl/cstdint.hpp>
#include <sl/serializer.hpp>
#include <iostream>

namespace sl {

  /**
   * A duration of time.
   */
  class time_duration {

  protected:
    int64_t value;

  public: // Init

    /// Default initialization.
    inline time_duration() { value = 0; }

    /// Initialization to usec.
    inline explicit time_duration(int64_t usec) { value = usec; }

  public: // Queries

    /// This converted to microseconds
    inline int64_t as_microseconds() const { return value; }

    /// This converted to milliseconds
    inline int64_t as_milliseconds() const { return value/1000; }

    /// This converted to seconds
    inline int64_t as_seconds() const { return value/1000/1000; }

  public: // Serialization
    
    inline void store_to(output_serializer& s) const {
      s << value;
    }
    
    inline void retrieve_from(input_serializer& s) {
      s >> value;
    }

  public: // Comparison

    /// Is this equal to other?
    inline bool operator == (const time_duration& other) const {
      return value == other.value;
    }

    /// Is this less than other?
    inline bool operator < (const time_duration& other) const {
      return value < other.value;
    }
      
    SL_OP_COMPARABLE1(time_duration);          // provides >, <=, >=
    SL_OP_EQUALITY_COMPARABLE1(time_duration); // provides !=
      
  public: // Algebra
      
    /// Increment this by other.
    inline time_duration operator += (const time_duration& other) {
      value += other.value;
      return *this;
    }

    /// Decrement this by other.
    inline time_duration operator -= (const time_duration& other) {
      value -= other.value;
      return *this;
    }
      
    /// Scale this by scalar.
    inline time_duration operator *= (double scalar) {
      value = static_cast<int64_t>(value * scalar);
      return *this;
    }
      
    /// Divide this by scalar.
    inline time_duration operator /= (double scalar) {
      value = static_cast<int64_t>(value * scalar);
      return *this;
    }
      
    SL_OP_LINEAR_SPACE(time_duration, double);
      
  public: // Useful constants
      
    /// A duration of length zero.
    inline static time_duration zero()            { return time_duration(0); }
    

    /// A duration of length one microsecond.
    inline static time_duration one_microsecond() { return time_duration(1); };

    /// A duration of length one millisecond.
    inline static time_duration one_millisecond() { return one_microsecond() * 1000.; }

    /// A duration of length one second.
    inline static time_duration one_second()      { return one_millisecond() * 1000.; }

    /// A duration of length one minute.
    inline static time_duration one_minute()      { return one_second() * 60.; }

    /// A duration of length one hour.
    inline static time_duration one_hour(void)    { return one_minute() * 60.; }

    /// A duration of length one day.
    inline static time_duration one_day(void)     { return one_hour() * 24.; }

  };

  
};

// Arithmetic operations overload

SL_OP_LINEAR_SPACE_OVERLOADS(SL_OP_NO_TEMPLATE,sl::time_duration, double);

// I/O


/// Write duration to to std::ostream s.
std::ostream& operator<<(std::ostream& s, const sl::time_duration& t);

#endif
