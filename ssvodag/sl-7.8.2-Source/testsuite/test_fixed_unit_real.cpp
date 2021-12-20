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
#include <sl/fixed_unit_real.hpp>
#include <iostream>

//--------------------

static std::size_t failed_test_count = 0;

//////////////////////////////////////////////////////////////////////

template <std::size_t N_Bits, bool B_Signed>
class test_signed_fixed_unit_real {
public:
  typedef sl::fixed_unit_real<N_Bits,B_Signed> real_t;
  typedef sl::interval<double>             interval_t;
 
  static void dump_math() {
    std::cerr << "  size                          = " << sizeof(real_t)                                          << std::endl;
    std::cerr << "  epsilon                       = " << sl::scalar_math<real_t>::epsilon()                      << std::endl;
    std::cerr << "  eta                           = " << sl::scalar_math<real_t>::eta()                          << std::endl;
    std::cerr << "  basic_op_max_error_ulp        = " << sl::scalar_math<real_t>::basic_op_max_error_ulp()       << std::endl;
    std::cerr << "  basic_op_max_error            = " << sl::scalar_math<real_t>::basic_op_max_error()           << std::endl;
    std::cerr << "  basic_op_zero_test            = " << sl::scalar_math<real_t>::basic_op_zero_test()           << std::endl;
    std::cerr << "  one_minus_basic_op_max_error  = " << sl::scalar_math<real_t>::one_minus_basic_op_max_error() << std::endl;
    std::cerr << "  one_plus_basic_op_max_error   = " << sl::scalar_math<real_t>::one_plus_basic_op_max_error()  << std::endl;
  }

  static void do_it() {
    sl::tester tester(sl::numeric_traits<real_t>::what());
    
    tester.test("add: 0.1 + 0.1", (real_t(0.1) + real_t(0.1)).value(), interval_t("0.2"));
    tester.test("sub: 0.2 - 0.1", (real_t(0.2) - real_t(0.1)).value(), interval_t("0.1"));
    tester.test("mul: 0.9 * 0.9", (real_t(0.9) * real_t(0.9)).value(), interval_t("0.8"));
    tester.test("div: 0.2 / 0.5", (real_t(0.2) / real_t(0.5)).value(), interval_t("0.4"));

    if (real_t::is_signed) {
      tester.test("add: 0.1 + -0.2", (real_t(0.1) + real_t(-0.2)).value(), interval_t("-0.1"));
      tester.test("sub: 0.2 -  0.3", (real_t(0.2) - real_t( 0.3)).value(), interval_t("-0.1"));
      tester.test("mul: 0.9 * -0.9", (real_t(0.9) * real_t(-0.9)).value(), interval_t("-0.8"));
      tester.test("div: 0.2 / -0.5", (real_t(0.2) / real_t(-0.5)).value(), interval_t("-0.4"));

      tester.test("ofl: 0.8 + 0.9", (real_t( 0.8) + real_t(0.9)).value(), interval_t("1.0"));
      tester.test("ofl: -0.8 - 0.9", (real_t(-0.8) - real_t(0.9)).value(), interval_t("-1.0"));
      tester.test("ofl: 0.8 / 0.1", (real_t( 0.8) / real_t(0.1)).value(), interval_t("1.0"));
    } else {
      tester.test("ofl: 0.8 + 0.9", (real_t( 0.8) + real_t(0.9)).value(), interval_t("1.0"));
      tester.test("ofl: 0.8 - 0.9", (real_t( 0.8) - real_t(0.9)).value(), interval_t("0.0"));
      tester.test("ofl: 0.8 / 0.1", (real_t( 0.8) / real_t(0.1)).value(), interval_t("1.0"));
    }

    failed_test_count += tester.failed_test_count();
  };

};

//////////////////////////////////////////////////////////////////////

int main() {
  test_signed_fixed_unit_real<8,true>::do_it();
  test_signed_fixed_unit_real<16,true>::do_it();
  test_signed_fixed_unit_real<32,true>::do_it();

  test_signed_fixed_unit_real<8,false>::do_it();
  test_signed_fixed_unit_real<16,false>::do_it();
  test_signed_fixed_unit_real<32,false>::do_it();

  return (int)failed_test_count;
}







