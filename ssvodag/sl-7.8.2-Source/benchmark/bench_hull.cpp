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
#include <sl/convex_hull.hpp>
#include <vector>
#include <stdlib.h>

static inline void random_seed(unsigned int seed) {
  srand(seed);
}

static inline sl::point2f random_uniform_point() {
  float x = (rand() / (float)(RAND_MAX));
  float y = (rand() / (float)(RAND_MAX));
  return sl::point2f(x,y);
}


// -------------------------------------------

std::vector<sl::point2f> uniform_points(std::size_t sz) {
  std::vector<sl::point2f> result;

  random_seed(1);
  for (std::size_t i=0; i<sz; ++i) {
    result.push_back(random_uniform_point());
  }
  return result;
}

// -------------------------------------------

static void benchmark_hull(const char* id,
			   const std::vector<sl::point2f>& points) {
  sl::convex_hull_builder<2,float, std::vector<sl::point2f> > ch_builder;

    const size_t hull_method_count = 3;

    const char* hull_method_id[hull_method_count] = {
      "DEFAULT ALGORITHM",
      "CHAIN ALGORITHM",
      "INCREMENTAL ALGORITHM"
    };

    for (size_t hull_method = 0; hull_method<hull_method_count; ++hull_method) {
      
      int ch_create_ms = 0;
      {
	sl::cpu_time_clock clk;
	clk.restart();
	switch (hull_method) {
	case 0: ch_builder.build(points); break;
	case 1: ch_builder.build_chain(points); break;
	case 2: ch_builder.build_incremental(points); break;
	default: SL_FAIL("Unknown hull method");
	}
	ch_create_ms = (int)(clk.elapsed().as_microseconds()/1000);
      }

      std::cout << std::endl;
      std::cout << "== HULL - " << id << " - " << hull_method_id[hull_method] << std::endl;
      std::cout << "   count   : " << points.size()   << " points" << std::endl;
      std::cout << "   create  : " << ch_create_ms    << " ms" << std::endl;
    }
}

// -------------------------------------------

int main() {
  benchmark_hull("Uniform 100", uniform_points(100));
  benchmark_hull("Uniform 1'000", uniform_points(1000));
  benchmark_hull("Uniform 10'000", uniform_points(10000));
  benchmark_hull("Uniform 100'000", uniform_points(100000));
}
