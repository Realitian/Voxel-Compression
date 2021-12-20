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
#ifndef SL_BENCHMARKER_HPP
#define SL_BENCHMARKER_HPP

#include <sl/clock.hpp>

namespace sl {

  /// An abstract operations
  class operation {
  public:
    virtual ~operation() {}
    
    /// Initialize - must be done once before the first call to operate()
    virtual void initialize() {}

    /// Operate - might be done many time after a single call to initialize()
    virtual void operate() = 0;

    /// Finalize - must be done once after 0-N calls to operate()
    virtual void finalize() {}
  };

  class null_operation: public operation {
  public:
    virtual ~null_operation() {}
    virtual void initialize() {}
    virtual void operate() {};
    virtual void finalize() {}
  };


  /// A class for measuring how much time is spent in an operation
  class benchmarker {
  protected:
    double last_cpu_time_rate_;
    double last_real_time_rate_;

    cpu_time_clock cpu_time_clock_;
    real_time_clock real_time_clock_;

  public:

    /// Default init
    benchmarker() { 
      last_cpu_time_rate_ = 0.;
      last_real_time_rate_ = 0.;
    }
    
    /// Execute a benchmark for operation cmd.
    void measure_rate(operation &cmd);
    
    /// The last measured number operations/CPU second
    inline double last_cpu_time_rate() {
      return last_cpu_time_rate_;
    }

    /// The last measured number operations/real second
    inline double last_real_time_rate() {
      return last_real_time_rate_;
    }
    
  }; //benchmarker

}; // namespace sl

#endif
