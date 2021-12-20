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
#include <sl/benchmarker.hpp>
#include <sl/utility.hpp>

using namespace sl;

void benchmarker::measure_rate(operation &cmd) {
  int reps;
  time_duration finalization_cpu_time_overhead;
  time_duration finalization_real_time_overhead;
  time_duration real_time_elapsed;
  time_duration cpu_time_elapsed;
	
  // Choose minimum run time
  time_duration min_run_time = 
    median(100.0 * max(real_time_clock_.resolution(),
		       cpu_time_clock_.resolution()),
	   0.5 * time_duration::one_second(),
	   5.0 * time_duration::one_second());

  // Calibrate, measuring overhead of finalization and timing
  cmd.initialize();
  reps = 0;
  cpu_time_clock_.restart();
  real_time_clock_.restart();
  do {
    cmd.finalize();
    reps++;
    real_time_elapsed = real_time_clock_.elapsed();
    cpu_time_elapsed = cpu_time_clock_.elapsed();
  } while (real_time_elapsed < min_run_time);
  finalization_real_time_overhead = (1.0/(double)reps) * real_time_elapsed;
  finalization_cpu_time_overhead = (1.0/(double)reps) * cpu_time_elapsed;
	
  // Measure
  reps = 0;
  int reps_factor = 2;

  time_duration min_expected_run_time = min_run_time + finalization_real_time_overhead;
  do {
    reps = max(reps * reps_factor, reps+1);

    cmd.initialize();
    cpu_time_clock_.restart();
    real_time_clock_.restart();
    {
      for (int i = reps; i > 0; --i) {
	cmd.operate();
      }
      cmd.finalize();
    }
    real_time_elapsed = real_time_clock_.elapsed();
    cpu_time_elapsed = cpu_time_clock_.elapsed();

    reps_factor = static_cast<int>(((double) min_expected_run_time.as_microseconds()) / 
				   max((double) real_time_elapsed.as_microseconds(), 1.0));
    reps_factor = median(reps_factor+1, 2, 10);
  } while (real_time_elapsed < min_expected_run_time);
	
  last_real_time_rate_ = (double)reps / (1.E-6 * (real_time_elapsed - finalization_real_time_overhead).as_microseconds());
  last_cpu_time_rate_  = (double)reps / (1.E-6 * (cpu_time_elapsed - finalization_cpu_time_overhead).as_microseconds());
}

