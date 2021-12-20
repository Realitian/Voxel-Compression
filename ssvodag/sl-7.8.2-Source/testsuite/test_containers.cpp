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
#include <sl/utility.hpp>
#include <sl/sorted_vector_set.hpp>
#include <sl/adaptive_packed_memory_array.hpp>

#include <sl/circular_buffer.hpp>

#include <iostream>
#include <set>
#include <map>
#include <algorithm>

static std::size_t failed_test_count = 0;

std::ostream& operator<< (std::ostream& os, const std::pair<int,int>& p) {
  return os << "(" << p.first << " " << p.second << ")"; 
}

template <class CONTAINER_T>
std::string as_string(CONTAINER_T& container) {
  std::ostringstream os;
  os << "{";
  for (typename CONTAINER_T::iterator it=container.begin(); it!= container.end(); ++it) {
    os << " " << *it;
  }
  os << " }";
  return os.str();
}

void do_test_circular_buffer() {
  sl::tester tester("sl::circular_buffer<int>");

  sl::circular_buffer<int> cb;

  tester.test("Creation", as_string(cb), std::string("{ }"));

  cb.set_capacity(4);

  tester.test("Set capacity", as_string(cb), std::string("{ }"));

  for (std::size_t i=0; i<5; ++i) {
    cb.push_back(i);
  }

  tester.test("Insert with wrap", as_string(cb), std::string("{ 1 2 3 4 }"));

  cb.push_front(0);
  tester.test("Push front with wrap", as_string(cb), std::string("{ 0 1 2 3 }"));

  cb.pop_front();
  tester.test("Pop front", as_string(cb), std::string("{ 1 2 3 }"));

  cb.pop_back();
  tester.test("Pop back", as_string(cb), std::string("{ 1 2 }"));

  tester.test("Size", cb.size(), std::size_t(2));

  cb.set_capacity(5);
  tester.test("Reset capacity -- grow/1", as_string(cb), std::string("{ 1 2 }"));
  cb.push_back(0);
  tester.test("Reset capacity -- grow/2", as_string(cb), std::string("{ 1 2 0 }"));

  for (std::size_t i=1; i<5; ++i) {
    cb.push_back(i);
  }
  tester.test("Reset capacity -- grow/2", as_string(cb), std::string("{ 0 1 2 3 4 }"));

  cb.set_capacity(4);
  tester.test("Reset capacity -- shrink", as_string(cb), std::string("{ 0 1 2 3 }"));

  cb.resize(8, 11);
  tester.test("Resize", as_string(cb), std::string("{ 0 1 2 3 11 11 11 11 }"));
  
  cb.rotate(cb.begin()+2);
  tester.test("Rotate full", as_string(cb), std::string("{ 2 3 11 11 11 11 0 1 }"));

  cb.pop_back();
  cb.pop_front();
  tester.test("Pop back+front", as_string(cb), std::string("{ 3 11 11 11 11 0 }"));

  cb.rotate(cb.begin()+2);
  tester.test("Rotate non full", as_string(cb), std::string("{ 11 11 11 0 3 11 }"));

  failed_test_count += tester.failed_test_count();
}
  
template <class CONTAINER_T, class REFCONTAINER_T>
void do_test_set(const std::string& s) {
  sl::tester tester(s);

  typedef typename CONTAINER_T::value_type value_t;
  
  REFCONTAINER_T refset;
  CONTAINER_T    myset;

  std::vector<value_t> data;
  for (std::size_t i=0; i<100; ++i) {
    data.push_back(i);
    if (i%4==0) data.push_back(i);
  }

  tester.test("Creation", as_string(myset), as_string(refset));

  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    myset.insert(x);
    refset.insert(x);
  }
  tester.test("Insert ordered size", myset.size(), refset.size());
  tester.test("Insert ordered", as_string(myset), as_string(refset));

  myset.clear();
  refset.clear();
  tester.test("Clear", as_string(myset), as_string(refset));

  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[data.size()-1-i];
    myset.insert(x);
    refset.insert(x);
  }
  tester.test("Insert back ordered size", myset.size(), refset.size());
  tester.test("Insert back ordered", as_string(myset), as_string(refset));

  myset.clear();
  refset.clear();
  tester.test("Clear", as_string(myset), as_string(refset));

  std::random_shuffle(data.begin(), data.end());
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    myset.insert(x);
    refset.insert(x);
  }
  tester.test("Insert random size", myset.size(), refset.size());
  tester.test("Insert random", as_string(myset), as_string(refset));

  std::random_shuffle(data.begin(), data.end());
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    if (i%5==0) {
      myset.erase(x);
      refset.erase(x);
    }
  }
  tester.test("Erase random size", myset.size(), refset.size());
  tester.test("Erase random", as_string(myset), as_string(refset));

  myset.clear();
  refset.clear();
  tester.test("Clear", as_string(myset), as_string(refset));

  failed_test_count += tester.failed_test_count();
}

template <class CONTAINER_T, class REFCONTAINER_T>
void do_test_map(const std::string& s) {
  sl::tester tester(s);

  typedef typename CONTAINER_T::key_type key_t;
  typedef typename CONTAINER_T::mapped_type mapped_t;
  typedef typename CONTAINER_T::value_type value_t;
  
  REFCONTAINER_T refmap;
  CONTAINER_T    mymap;

  std::vector< std::pair<key_t,mapped_t> > data;
  for (int i=0; i<100; ++i) {
    data.push_back(std::make_pair(key_t(i),mapped_t(2*i)));
    if (i%4==0) data.push_back(std::make_pair(key_t(i),mapped_t(2*i)));
  }

  tester.test("Creation", as_string(mymap), as_string(refmap));

  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    mymap.insert(x);
    refmap.insert(x);
  }
  tester.test("Insert ordered size", mymap.size(), refmap.size());
  tester.test("Insert ordered", as_string(mymap), as_string(refmap));

  mymap.clear();
  refmap.clear();
  tester.test("Clear", as_string(mymap), as_string(refmap));

  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[data.size()-1-i];
    mymap.insert(x);
    refmap.insert(x);
  }
  tester.test("Insert back ordered size", mymap.size(), refmap.size());
  tester.test("Insert back ordered", as_string(mymap), as_string(refmap));

  mymap.clear();
  refmap.clear();
  tester.test("Clear", as_string(mymap), as_string(refmap));

  std::random_shuffle(data.begin(), data.end());
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    mymap.insert(x);
    refmap.insert(x);
  }
  tester.test("Insert random size", mymap.size(), refmap.size());
  tester.test("Insert random", as_string(mymap), as_string(refmap));

  mymap.clear();
  refmap.clear();
   std::random_shuffle(data.begin(), data.end());
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    mymap[x.first] = x.second;
    refmap[x.first] = x.second; // FIXME [x.first] = x.second;
  }
  tester.test("Insert random with [] size", mymap.size(), refmap.size());
  tester.test("Insert random with []", as_string(mymap), as_string(refmap));

  std::random_shuffle(data.begin(), data.end());
  for (std::size_t i=0; i<data.size(); ++i) {
    value_t x = data[i];
    if (i%5==0) {
      mymap.erase(x.first);
      refmap.erase(x.first);
    }
  }
  tester.test("Erase random size", mymap.size(), refmap.size());
  tester.test("Erase random", as_string(mymap), as_string(refmap));

  mymap.clear();
  refmap.clear();
  tester.test("Clear", as_string(mymap), as_string(refmap));

  failed_test_count += tester.failed_test_count();
}

int main() {
  do_test_set<sl::sorted_vector_set<int>, std::set<int> >("sl::sorted_vector_set<int>");
  do_test_set<sl::sorted_vector_multiset<int>, std::multiset<int> >("sl::sorted_vector_multiset<int>");
  do_test_set<sl::adaptive_packed_memory_array_set<int>, std::set<int> >("sl::adaptive_packed_memory_array_set<int>");
  do_test_set<sl::adaptive_packed_memory_array_multiset<int>, std::multiset<int> >("sl::adaptive_packed_memory_array_multiset<int>");

  do_test_map<sl::adaptive_packed_memory_array_map<int,int>, std::map<int,int> >("sl::adaptive_packed_memory_array_map<int>");
  // do_test_map<sl::adaptive_packed_memory_array_multimap<int,int>, std::multimap<int,int> >("sl::adaptive_packed_memory_array_multimap<int>");

  
  do_test_circular_buffer();
  
  return failed_test_count;
}

  
  
  

  
  
