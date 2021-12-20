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
#ifndef SL_CLOCK_HPP
#define SL_CLOCK_HPP

#include <sl/time_point.hpp> // For time_point and time_duration

namespace sl {

  /// Abstract clock for measuring time
  class clock {

    sl::time_duration resolution_;


  public: // Creation

    inline clock() {
    }

    virtual inline ~clock() {
    }
    
  public: // Commands

    /// reset clock
    virtual void restart() = 0;

    /// wait until next tick
    inline void wait_for_tick() const { 
      sl::time_duration d0 = elapsed();
      do { } while (elapsed() == d0);
    }

  public: // Queries

    /// the amount of time elapsed since last restart
    virtual sl::time_duration elapsed() const = 0;

    /// the resolution of the clock (minimum measurable duration)
    virtual sl::time_duration resolution() const {
      if (resolution_ == sl::time_duration(0)) {
	// Default estimation method
	for (int i=0; i<10; i++) {
	  wait_for_tick();
	  sl::time_duration d0 = elapsed();
	  wait_for_tick();
	  sl::time_duration d1 = elapsed();
	  sl::time_duration delta = d1 - d0;
	  if (i==0 || delta < resolution_) {
	    ((clock*)this)->resolution_ = delta;
	  }
	}
      }
      return resolution_;
    }

  };

  /// Clock measuring real time (a.k.a wall clock)
  class real_time_clock: public clock {
  protected: // Implementation
    sl::time_point reference_value_;
    
  public: // Access to wall clock time

    /// Current time
    static sl::time_point now();

  public: // Creation

    /// default init
    inline real_time_clock() {
      reference_value_ = real_time_clock::now(); 
    }
    inline virtual ~real_time_clock() {
    }

  public: // Clock implementation

    inline void restart() { 
      reference_value_ = real_time_clock::now(); 
    }

    inline sl::time_duration elapsed() const {
      return real_time_clock::now() - reference_value_;
    }
  };

  /// Clock measuring CPU time
  class cpu_time_clock: public clock {
  protected: // Implementation
    sl::time_duration reference_value_;
    
  public: // Access to CPU clock time

    static sl::time_duration now();

  public: // Creation

    inline cpu_time_clock() {
      reference_value_ = cpu_time_clock::now(); 
    }

    virtual inline ~cpu_time_clock() {
    }

  public: // Clock implementation

    inline void restart() { 
      reference_value_ = cpu_time_clock::now(); 
    }

    inline sl::time_duration elapsed() const {
      return cpu_time_clock::now() - reference_value_;
    }
  };
        
};

#endif
