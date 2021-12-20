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
#include <sl/interval.hpp>
#include <iostream>

//--------------------

static std::size_t failed_test_count = 0;

//////////////////////////////////////////////////////////////////////

template <class T_SCALAR>
class test_interval {
public:
  typedef sl::interval<T_SCALAR> interval_t;
 
  static void dump_math(const std::string& scalar_id) {
    std::cerr << "sl::scalar_math<" << scalar_id << ">" << std::endl;
    std::cerr << "  size                          = " << sizeof(T_SCALAR)                                          << std::endl;
    std::cerr << "  epsilon                       = " << sl::scalar_math<T_SCALAR>::epsilon()                      << std::endl;
    std::cerr << "  eta                           = " << sl::scalar_math<T_SCALAR>::eta()                          << std::endl;
    std::cerr << "  basic_op_max_error_ulp        = " << sl::scalar_math<T_SCALAR>::basic_op_max_error_ulp()       << std::endl;
    std::cerr << "  basic_op_max_error            = " << sl::scalar_math<T_SCALAR>::basic_op_max_error()           << std::endl;
    std::cerr << "  basic_op_zero_test            = " << sl::scalar_math<T_SCALAR>::basic_op_zero_test()           << std::endl;
    std::cerr << "  one_minus_basic_op_max_error  = " << sl::scalar_math<T_SCALAR>::one_minus_basic_op_max_error() << std::endl;
    std::cerr << "  one_plus_basic_op_max_error   = " << sl::scalar_math<T_SCALAR>::one_plus_basic_op_max_error()  << std::endl;
  }

  static void do_it(const std::string& scalar_id) {
    dump_math(scalar_id);

    sl::tester tester("interval<" + scalar_id + ">");
    
    tester.test("[10..10] * [20..20] | inf",  (interval_t("10.000") *  interval_t("20.000"))[0], interval_t("200.0"));
    tester.test("[10..10] * [20..20] | sup",  (interval_t("10.000") *  interval_t("20.000"))[1], interval_t("200.0"));
    
    tester.test("abs([-15..10]) | inf",  sl::abs(interval_t(-15,10))[0], interval_t("0.0"));
    tester.test("abs([-15..10]) | sup",  sl::abs(interval_t(-15,10))[1], interval_t("15.0"));

    interval_t a = interval_t(1.,2.);
    interval_t b = interval_t(3.,4.);
    
    tester.test("a + b | inf",  (a+b)[0], interval_t("4.0000"));
    tester.test("a + b | sup",  (a+b)[1], interval_t("6.0000"));
    
    tester.test("a ** 3 | inf", std::pow(a,3)[0], interval_t("1.0000"));
    tester.test("a ** 3 | sup", std::pow(a,3)[1], interval_t("8.0000"));

    tester.test("a ** b | inf", std::pow(a,b)[0], interval_t("1.0000"));
    tester.test("a ** b | sup", std::pow(a,b)[1], interval_t("16.0000"));
    
    tester.test("cos(a) | inf", std::cos(a)[0], interval_t("-0.416"));
    tester.test("cos(a) | sup", std::cos(a)[1], interval_t("+0.540"));
    
    sl::interval<double> a_double = sl::interval<double>(a);
    tester.test("a_double | inf",  a_double[0], sl::interval<double>("1.0000"));
    tester.test("a_double | sup",  a_double[1], sl::interval<double>("2.0000"));

    failed_test_count += tester.failed_test_count();
  };

};

//////////////////////////////////////////////////////////////////////

int main() {
  test_interval<float>::do_it("float");
  test_interval<double>::do_it("double");

  return (int)failed_test_count;
}







