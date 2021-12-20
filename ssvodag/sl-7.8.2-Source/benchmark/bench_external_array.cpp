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
#include "sl/sort.hpp"
#include "sl/external_array.hpp"
#include <sl/utility.hpp>
#include <sl/clock.hpp>
#include <sl/random.hpp>
#include <iostream>
#include <iterator>
#include <cstdlib>

void xsort(sl::external_array1<float>* xarray, int sorttp) {
  sl::real_time_clock ck;
  ck.restart();
  switch (sorttp) {
  case 0:
    std::cerr << "### METHOD= sl::quicksort" << std::endl;
    sl::quicksort(xarray->begin(), xarray->end());
    break;
  default:
    std::cerr << "### METHOD= std::sort" << std::endl;
    std::sort(xarray->begin(), xarray->end());
  }
  std::cerr << "==== Done: Time=" << ck.elapsed().as_seconds() << "s" << std::endl;
  std::cerr << "Checking order..." << std::endl;
  std::cerr << (sl::is_sorted(xarray->begin(), xarray->end()) ? "OK" : "ERROR!!") << std::endl;
}

int main(int argc, const char* argv[]) {
  if (argc!=3) {
    std::cerr << "Usage: " << argv[0] << " <sorttp> <thousands of elements>" << std::endl;
    return 1;
  }

  std::size_t N = atoi(argv[2])*1000;
  std::size_t sorttp = atoi(argv[1]);

  sl::external_array1<float>* xarray = new sl::external_array1<float>("test.arr","w");
  xarray->clear();

  sl::real_time_clock ck;
  std::cerr << "==== Inserting " << sl::human_readable_quantity(N) << " random elements..." << std::endl;
  ck.restart();
  sl::random::uniform<float> rng;
  for (std::size_t i=0; i<=N; ++i) {
    xarray->push_back(rng.value());
  }
  std::cerr << "Done: Time=" << ck.elapsed().as_seconds() << "s" << std::endl;
  
  std::cerr << "==== Sorting " << sl::human_readable_quantity(N) << " random elements..." << std::endl;
  xsort(xarray, sorttp);
  
  std::cerr << "==== Sorting " << sl::human_readable_quantity(N) << " ordered elements..." << std::endl;
  xsort(xarray, sorttp);

  delete xarray; xarray=0;

  return 0;
}
