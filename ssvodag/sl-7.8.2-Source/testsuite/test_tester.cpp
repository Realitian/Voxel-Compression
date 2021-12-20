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

//--------------------

static std::size_t failed_test_count = 0;

////////////////////////////////////////////////////////////////////

SL_COMPILE_TIME_CHECK("true", true);

////////////////////////////////////////////////////////////////////

class test_tester {
public:
  static const char * abc() {
    return "abc";
  }
  
  static const char * bcd() {
    return "bcd";
   }
  
  static void do_it() {
    sl::tester t("Tested Tester");
    t.set_silent(true);

    sl::tester tester("sl::tester");
    
    t.reset();
    t.test("String-OK",   abc(), abc());
    tester.test("String-OK", t.failed_test_count(), std::size_t(0));

    t.reset();
    t.test("String-Fail", abc(), bcd());
    tester.test("String-Fail", t.failed_test_count(), std::size_t(1));

    t.reset();
    t.test("int-OK", 10, 10);
    tester.test("int-OK", t.failed_test_count(), std::size_t(0));

    t.reset();
    t.test("int-Fail", 10, 11);
    tester.test("int-Fail", t.failed_test_count(), std::size_t(1));
    
    t.reset();
    t.test("Interval-Ok", 10.0, sl::interval<double>("10.0"));
    tester.test("interval-OK", t.failed_test_count(), std::size_t(0));

    t.reset();
    t.test("Interval-Fail", 10.0, sl::interval<double>("9.0"));
    tester.test("interval-Fail", t.failed_test_count(), std::size_t(1));

    t.reset();

    failed_test_count += tester.failed_test_count();
  }
}; // test_tester

int main() {
  test_tester::do_it();

  return (int)failed_test_count;
}







