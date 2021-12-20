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
#include <sl/math.hpp>
#include <sl/cie.hpp>

static std::size_t failed_test_count = 0;

template <class T_SCALAR>
class test_cie {
  typedef T_SCALAR               value_t;
  typedef sl::interval<value_t> interval_t;
public:
  static void do_it() {
    sl::tester tester("cie < " + sl::numeric_traits<value_t>::what() + " >");
    
    sl::cie<value_t> cie;

    value_t rgb_in[3];
    value_t rgb_out[3];
    value_t xyz_out[3];

    rgb_in[0] = value_t(0.5);
    rgb_in[1] = value_t(0.6);
    rgb_in[2] = value_t(0.7);

    cie.rgb_to_xyz(rgb_in[0], rgb_in[1], rgb_in[2],
		   xyz_out[0], xyz_out[1], xyz_out[2]);
    cie.xyz_to_rgb(xyz_out[0], xyz_out[1], xyz_out[2],
		   rgb_out[0], rgb_out[1], rgb_out[2]);

    tester.test("xyz->rgb R", rgb_in[0] - rgb_out[0], interval_t("0.000"));
    tester.test("xyz->rgb G", rgb_in[1] - rgb_out[1], interval_t("0.000"));
    tester.test("xyz->rgb B", rgb_in[2] - rgb_out[2], interval_t("0.000"));

    failed_test_count += tester.failed_test_count();
  } 
};

template <class T_SCALAR>
class test_suite {
public:

  static void do_it() {
    test_cie<T_SCALAR>::do_it();
  }    

};

int main() {
  test_suite<float>::do_it();
  test_suite<double>::do_it();

  return (int)failed_test_count;
}


