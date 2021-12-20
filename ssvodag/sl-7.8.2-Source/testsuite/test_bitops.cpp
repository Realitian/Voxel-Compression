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
#include <sl/bitops.hpp>

//--------------------

static std::size_t failed_test_count = 0;

template <class G_int>
class test_bitops {
public:
  typedef G_int                                                value_t;
  typedef sl::bitops<value_t>                                  bitops;
  typedef sl::morton_bitops<value_t,2>                         morton2d;
  typedef sl::morton_bitops<value_t,3>                         morton3d;

public:

  static void do_it() {
    sl::tester tester("bitops <" + 
		      sl::numeric_traits<value_t>::what() + 
		      " >");
    
    tester.test("one_count(0xf0)",  bitops::one_count(value_t(0xf0)), value_t(4));
    tester.test("zero_count(0xf0)", bitops::zero_count(value_t(0xf0)), value_t(bitops::bits - 4));
    tester.test("reverse",          (bitops::reverse(value_t(0x10)) != value_t(0x10)) &&
		(bitops::reverse(bitops::reverse(value_t(0x10))) == 0x10));
    tester.test("low_order_zero_count(0xf0)", bitops::one_count(value_t(0xf0)), value_t(4));
    tester.test("log2(16)", bitops::log2(16), value_t(4));
    tester.test("is_power2(16)", (bitops::is_power2(16)));
    tester.test("is_power2(15)", !(bitops::is_power2(15)));
    tester.test("next_power2(16)", bitops::next_power2(15), value_t(16));
    tester.test("next_power2(16)", bitops::next_power2(16), value_t(16));
    tester.test("parity(16)", bitops::parity(16), true);
    tester.test("parity(15)", bitops::parity(15), false);

    tester.test("gray(19)", bitops::binary_from_gray(bitops::gray_from_binary(19)), value_t(19));
    tester.test("gray(254)", bitops::binary_from_gray(bitops::gray_from_binary(254)), value_t(254));
    tester.test("gray(255)", bitops::binary_from_gray(bitops::gray_from_binary(255)), value_t(255));

    if (sizeof(value_t)>1) {
      tester.test("morton2d(11,12)_x", morton2d::decoded((morton2d::encoded(11, 0) |
							morton2d::encoded(12, 1)), 0), value_t(11));
      tester.test("morton2d(11,12)_y", morton2d::decoded((morton2d::encoded(11, 0) |
							  morton2d::encoded(12, 1)), 1), value_t(12));
      tester.test("morton3d(11,12,13)_x", morton3d::decoded((morton3d::encoded(11, 0) |
							     morton3d::encoded(12, 1) |
							     morton3d::encoded(13, 2)), 0), value_t(11));
      tester.test("morton3d(11,12,13)_y", morton3d::decoded((morton3d::encoded(11, 0) |
							     morton3d::encoded(12, 1) |
							     morton3d::encoded(13, 2)), 1), value_t(12));
      tester.test("morton3d(11,12,13)_z", morton3d::decoded((morton3d::encoded(11, 0) |
							     morton3d::encoded(12, 1) |
							     morton3d::encoded(13, 2)), 2), value_t(13));
      for (std::size_t i=0; i<5; ++i) {
	sl::uint32_t Full_bits   = 8*sizeof(value_t);
	sl::uint32_t Morton_bits = Full_bits/3;
	
	value_t x = (value_t(1)<<Morton_bits)-1-i;
	value_t y = (value_t(1)<<Morton_bits)-1-i-1;
	value_t z = (value_t(1)<<Morton_bits)-1-5*i-2;
	
	value_t mxyz = (morton3d::encoded(x,0) |
			morton3d::encoded(y,1) |
			morton3d::encoded(z,2));
	tester.test(std::string()+"morton3d("+sl::to_string(x)+","+sl::to_string(y)+","+sl::to_string(z)+").x",
		    morton3d::decoded(mxyz,0), x);
	tester.test(std::string()+"morton3d("+sl::to_string(x)+","+sl::to_string(y)+","+sl::to_string(z)+").y",
		    morton3d::decoded(mxyz,1), y);
	tester.test(std::string()+"morton3d("+sl::to_string(x)+","+sl::to_string(y)+","+sl::to_string(z)+").z",
		    morton3d::decoded(mxyz,2), z);
      }
    }
    failed_test_count += tester.failed_test_count();
  }
};

int main() {
  test_bitops<sl::uint8_t>::do_it();
  test_bitops<sl::uint16_t>::do_it();
  test_bitops<sl::uint32_t>::do_it();
  test_bitops<sl::uint64_t>::do_it();

  return (int)failed_test_count;

}


