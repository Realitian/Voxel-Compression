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
/////// ALWAYS TEST IN DEBUG MODE, UNLESS IT BREAKS THE COMPILER
# if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
# endif 
///////
#include <sl/tester.hpp>
#include <sl/clock.hpp>
#include <iostream>
#ifdef _WIN32
// ??? sleep ???
#else
#include <unistd.h> // for sleep
#endif

static std::size_t failed_test_count = 0;
//////////////////////////////////////////////////////////////////////

class test_time {
public:
  static void do_it() {
    sl::tester tester("real_time_clock & cpu_time_clock");
    
    sl::real_time_clock rt_clock;
    sl::cpu_time_clock  cpu_clock;
    
    tester.test("Real-time clock resolution", int(rt_clock.resolution().as_microseconds()) <= int(10000));
    tester.test("CPU -time clock resolution", int(rt_clock.resolution().as_microseconds()) <= int(10000));
    
    rt_clock.restart(); cpu_clock.restart();

#ifndef _WIN32
    // Only Unix sleeps...
    sleep(1);
#endif
    sl::time_duration rt_elapsed = rt_clock.elapsed();
    sl::time_duration cpu_elapsed = cpu_clock.elapsed();
#ifndef _WIN32
    tester.test("Real-time elapsed after 1 s sleep", rt_elapsed >= 0.9 * sl::time_duration::one_second());
    tester.test("CPU -time elapsed after 1 s sleep", cpu_elapsed <= 0.1 * sl::time_duration::one_second());
#endif

    rt_clock.restart(); cpu_clock.restart();
    do { } while (cpu_clock.elapsed() < sl::time_duration::one_second());
    rt_elapsed = rt_clock.elapsed();
    cpu_elapsed = cpu_clock.elapsed();
    
    tester.test("Real-time elapsed after 1 s loop", rt_elapsed >= 0.9 * sl::time_duration::one_second());
    tester.test("CPU -time elapsed after 1 s loop", cpu_elapsed >= 0.9 * sl::time_duration::one_second());

    failed_test_count += tester.failed_test_count();
  }
};

int main() {

  test_time::do_it();

  return (int)failed_test_count;
}



