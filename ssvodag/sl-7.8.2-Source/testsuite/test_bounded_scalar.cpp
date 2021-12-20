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
#include <sl/tester.hpp>
#include <sl/bounded_scalar.hpp>
#include <iostream>
//////////////////////////////////////////////////////////////////////

static std::size_t failed_test_count = 0;

template <typename T_SCALAR>
class test_bounded_scalar {
public:
  typedef sl::interval<T_SCALAR> interval_t;
  typedef sl::bounded_scalar<T_SCALAR> bounded_scalar_t;

  static void do_it(const std::string& scalar_id) {
    sl::tester tester("bounded_scalar<" + scalar_id + ">");
  
    tester.test("10 * 20", (bounded_scalar_t(10) *  bounded_scalar_t(20)).value(), interval_t("200.0"));
    tester.test("abs([-15])",  (sl::abs(bounded_scalar_t(-15))).value(), interval_t("15.0"));
    
    bounded_scalar_t a = bounded_scalar_t(2);
    bounded_scalar_t b = bounded_scalar_t(4);
    
    tester.test("a + b",  (a+b).value(), interval_t("6.0"));
    tester.test("a ** 3", std::pow(a,3).value(), interval_t("8.0"));
    tester.test("a ** b", std::pow(a,b).value(), interval_t("16.0"));
    tester.test("cos(a)", std::cos(a).value(), interval_t("-0.416"));

    failed_test_count += tester.failed_test_count();
  }
};

int main() {

  test_bounded_scalar<float>::do_it("float");
  test_bounded_scalar<double>::do_it("double");

  return (int)failed_test_count;
}

