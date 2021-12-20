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
#include <sl/oriented_box.hpp>
#include <vector>
#include <stdlib.h>

static inline void random_seed(unsigned int seed) {
  srand(seed);
}

static inline sl::point3f random_uniform_point() {
  float x = (rand() / (float)(RAND_MAX));
  float y = (rand() / (float)(RAND_MAX));
  float z = (rand() / (float)(RAND_MAX));
  return sl::point3f(x,y,z);
}


// -------------------------------------------

std::vector<sl::point3f> uniform_points(std::size_t sz) {
  std::vector<sl::point3f> result;

  random_seed(1);
  for (std::size_t i=0; i<sz; ++i) {
    result.push_back(random_uniform_point());
  }
  return result;
}

// -------------------------------------------

static void benchmark_obox(const char* id,
			   const std::vector<sl::point3f>& points) {
  sl::oriented_box_builder<3,float> obox_builder;

  sl::row_vector3f preferred_normal;

  int obox_create_no_normal_ms = 0;
  {
    sl::cpu_time_clock clk;
    clk.restart();
    obox_builder.begin_model();
    for (std::size_t i=0; i<points.size(); ++i) {
      obox_builder.put_point(points[i]);
    }
    obox_builder.end_model();
    preferred_normal = obox_builder.normal();
    obox_create_no_normal_ms = (int)(clk.elapsed().as_microseconds()/1000);
  }

  int obox_create_normal_ms = 0;
  {
    sl::cpu_time_clock clk;
    clk.restart();
    obox_builder.begin_model();
    obox_builder.set_normal(preferred_normal);
    for (std::size_t i=0; i<points.size(); ++i) {
      obox_builder.put_point(points[i]);
    }
    obox_builder.end_model();
    preferred_normal = obox_builder.normal();
    obox_create_normal_ms = (int)(clk.elapsed().as_microseconds()/1000);
  }

  std::cout << std::endl;
  std::cout << "== " << id << std::endl;
  std::cout << "   count                     : " << points.size() << " points" << std::endl;
  std::cout << "   create with cloud         : " << obox_create_no_normal_ms << " ms" << std::endl;
  std::cout << "   create with cloud+normal  : " << obox_create_normal_ms << " ms" << std::endl;
}

// -------------------------------------------

int main() {
  benchmark_obox("OBOX - Uniform 100", uniform_points(100));
  benchmark_obox("OBOX - Uniform 1'000", uniform_points(1000));
  benchmark_obox("OBOX - Uniform 10'000", uniform_points(10000));
  benchmark_obox("OBOX - Uniform 100'000", uniform_points(100000));
}
