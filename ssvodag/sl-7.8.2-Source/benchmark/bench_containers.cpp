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
#include <sl/utility.hpp>
#include <sl/sorted_vector_set.hpp>
#include <sl/adaptive_packed_memory_array.hpp>
#include <sl/utility.hpp>
#include <sl/clock.hpp>
#include <iostream>
#include <set>
#include <algorithm>
#include <iterator>
#include <cstdlib>

template <class CONTAINER_T>
void do_bench(const std::string& s,
	      std::size_t N) {
  std::cerr << "====================================================" << std::endl;
  std::cerr << "==== BENCHMARK: " << s << std::endl;
  std::cerr << "====================================================" << std::endl;

  typedef typename CONTAINER_T::value_type value_t;

  CONTAINER_T    myset;

  std::vector<value_t> data;
  for (std::size_t i=0; i<N; ++i) {
    data.push_back(i);
    if (i%4==0) data.push_back(i);
  }
  sl::real_time_clock ck;
  
  std::cerr << "---- Inserting " << sl::human_readable_quantity(data.size()) << " forward ordered elements... ";
  ck.restart();
  myset.clear();
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    myset.insert(x);
  }
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(data.size())/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;

  std::cerr << "---- Inserting " << sl::human_readable_quantity(data.size()) << " reverse ordered elements... ";
  myset.clear();
  ck.restart();
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[data.size()-1-i];
    myset.insert(x);
  }
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(data.size())/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;

  std::random_shuffle(data.begin(), data.end());
  std::cerr << "---- Inserting " << sl::human_readable_quantity(data.size()) << " random  ordered elements... ";
  myset.clear();
  ck.restart();
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    myset.insert(x);
  }
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(data.size())/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;

  std::cerr << "---- Traversing " << sl::human_readable_quantity(myset.size()) << " elements forward...        ";
  ck.restart();
  std::size_t traversed_count=0;
  for (typename CONTAINER_T::iterator it = myset.begin(); it!= myset.end(); ++it) {
    // value_t x = *it;
    ++traversed_count;
  }

  if (traversed_count !=myset.size()) std::cerr << "ERROR -- traversed " << traversed_count << " elts instead of " << myset.size() << "elts." << std::endl;
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(data.size())/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;
  std::cerr << "---- Traversing " << sl::human_readable_quantity(myset.size()) << " elements backward...       ";
  ck.restart();
  traversed_count=0;
  for (typename CONTAINER_T::reverse_iterator it = myset.rbegin(); it!= myset.rend(); ++it) {
    // value_t x = *it;
    ++traversed_count;
  }
  if (traversed_count !=myset.size()) std::cerr << "ERROR -- traversed " << traversed_count << " elts instead of " << myset.size() << "elts." << std::endl;
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(data.size())/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;

  std::random_shuffle(data.begin(), data.end());
  std::cerr << "---- Finding   " << sl::human_readable_quantity(data.size()) << " random  ordered elements... ";
  ck.restart();
  std::size_t found_count=0;
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    found_count += (myset.find(x) == myset.end()) ? 0 : 1;
  }
  if (found_count !=data.size()) std::cerr << "ERROR -- found " << found_count << " elts instead of " << data.size() << "elts." << std::endl;
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(data.size())/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;

  std::random_shuffle(data.begin(), data.end());
  const std::size_t M= N/3;
  std::cerr << "---- Erasing   " << sl::human_readable_quantity(M) << " random  ordered elements... ";
  ck.restart();
  for (std::size_t i=0; i<M; ++i) {
    value_t x = data[i];
    myset.erase(x);
  }
  std::cerr << "     SPEED: " << sl::human_readable_quantity(1000.0*double(M)/double(ck.elapsed().as_milliseconds())) << "items/s" << std::endl;
}

int main(int argc, const char* argv[]) {
  std::size_t N = 1000*1000;
  if (argc == 1) {
    // Use defaults
  } else if (argc==2) {
    N = atoi(argv[1])*1000;
  } else {
    std::cerr << "Usage: " << argv[0] << " <thousands of elements>" << std::endl;
    return 1;
  }
  
  do_bench<std::set<int> >("std::set<int>", N);
  do_bench<std::multiset<int> >("std::multiset<int>", N);
  do_bench<sl::adaptive_packed_memory_array_set<int> >("sl::adaptive_packed_memory_array_set<int>", N);
  do_bench<sl::adaptive_packed_memory_array_multiset<int> >("sl::adaptive_packed_memory_array_multiset<int>", N);
  // -- Too slow do_bench<sl::sorted_vector_set<int> >("sl::sorted_vector_set<int>", N);
  // -- Too slow do_bench<sl::sorted_vector_multiset<int> >("sl::sorted_vector_multiset<int>", N);
  
  return 0;
}
