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
/////// ALWAYS TEST IN DEBUG MODE
#if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////

#include <sl/utility.hpp>
#include <sl/tester.hpp>
#include <sl/floating_cone.hpp>

//--------------------

static std::size_t failed_test_count = 0;

//--------------------

typedef double value_t;
typedef sl::vector3d vector_t;
typedef sl::point3d point_t;
typedef sl::intervald interval_t;

void test_minifloating_cone() {
  const std::size_t d = 3;
  const std::size_t n = 10000;
  typedef sl::floating_cone<d,value_t> cone_t;

  sl::floating_cone_builder<d,value_t> mb;
  sl::random::uniform<value_t> rng;

  sl::tester tester(std::string() + "floating_cone_builder < " + sl::to_string(d) + ", " + sl::numeric_traits<value_t>::what() + " >");
  
  rng.set_seed(1999);

  vector_t n0;  n0.to_unit(0);
  value_t a0 = value_t(3.14f/4.0f);

  std::vector<vector_t> testset;
    
  for (std::size_t i=0; i<n; ++i) {
    sl::floating_cone_builder<d,value_t>::vector_t n;
    for (std::size_t j=0; j<d; ++j) {
      n[j] = rng.value();
    }
    n.ok_normalized();
    if (n.angle(n0) < a0) {
      testset.push_back(n0);
    }
  }

  mb.begin_model();
  for (std::vector<vector_t>::iterator it = testset.begin();
       it != testset.end();
       ++it) {
    mb.put_vector(*it);
  }
  mb.end_model();
  cone_t bcone = mb.last_bounding_cone(); 

  tester.test("mbp angle", bcone.angle() <= a0);
  for (std::vector<vector_t>::iterator it = testset.begin();
       it != testset.end();
       ++it) {
    tester.test("bounding cone containment", bcone.contains(*it));
  }
  failed_test_count += tester.failed_test_count();
}


int main() {
  test_minifloating_cone();
  return (int)failed_test_count;
}
