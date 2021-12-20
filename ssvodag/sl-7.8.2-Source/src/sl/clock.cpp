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
#include <sl/clock.hpp>

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/time.h>  // for struct timeval - CHECK PORTABILITY
#  include <sys/times.h> // for times()        - CHECK PORTABILITY 
#  include <unistd.h>    // for gettimeofday   - CHECK PORTABILITY
#endif
#include <time.h>        // for clock          - CHECK PORTABILITY

sl::time_point sl::real_time_clock::now() {
#ifdef _WIN32
  static bool inited = false;

  static LARGE_INTEGER  pc_freq;
  if (!inited) {
    QueryPerformanceFrequency(&pc_freq);
    inited = true;
  }

  LARGE_INTEGER  pc_now;
  QueryPerformanceCounter( &pc_now );
  return  sl::time_point(static_cast<uint64_t>(pc_now.QuadPart * 1000000 / pc_freq.QuadPart));
#else
#if 1
  struct timeval tv;  
  gettimeofday(&tv, 0L);
  return sl::time_point(static_cast<uint64_t>(tv.tv_usec) + 
			static_cast<uint64_t>(1000000) * 
			static_cast<uint64_t>(tv.tv_sec));
#else
  // <FIXME> Check real time clock.
  struct timespec tv;
  clock_gettime(CLOCK_REALTIME, &tv);
  return sl::time_point(static_cast<uint64_t>(tv.tv_nsec) / static_cast<uint64_t>(1000) + 
			static_cast<uint64_t>(1000000) * 
			static_cast<uint64_t>(tv.tv_sec));
#endif
#endif
}

sl::time_duration sl::cpu_time_clock::now() {
#ifdef _WIN32
  return sl::time_duration(static_cast<uint64_t>(::clock()) * 
			   (static_cast<uint64_t>(1000000) / 
			    static_cast<uint64_t>(CLOCKS_PER_SEC))); 
#else
  static uint64_t Microseconds_per_tick = 0;
  if (!Microseconds_per_tick) {
    Microseconds_per_tick = (static_cast<uint64_t>(1000000)/
                             static_cast<uint64_t>(::sysconf(_SC_CLK_TCK)));
  }
  struct tms timebuf;
  ::times(&timebuf);
  
  return sl::time_duration(static_cast<uint64_t>(timebuf.tms_utime + timebuf.tms_stime) * 
			   Microseconds_per_tick);
#endif
}
