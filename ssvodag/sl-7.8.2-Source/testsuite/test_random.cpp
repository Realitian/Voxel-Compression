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
# if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////
#include <sl/tester.hpp>
#include <sl/random.hpp>

static std::size_t failed_test_count = 0;

template <class T_SCALAR>
class test_random_uniform {
  typedef T_SCALAR               value_t;
  typedef sl::interval<value_t> interval_t;
public:
  static void do_it() {
    sl::tester tester("sl::random::uniform_open < " + sl::numeric_traits<value_t>::what() + " >");
    
    sl::random::uniform_open<value_t> rng;

    static const std::size_t N = 100;
    std::size_t d_hist[N];
    for (std::size_t i=0; i<N; ++i) {
      d_hist[i] = 0;
    }
    static const std::size_t M = 100000;
    for (std::size_t i=0; i<N*M; ++i) {
      ++d_hist[sl::uint32_t(rng.value()*N)];
    }
    for (std::size_t i=0; i<N; ++i) {
      tester.test("random bucket", value_t(double(d_hist[i])/double(M)), interval_t(value_t(0.98), value_t(1.02)));
    }

    failed_test_count += tester.failed_test_count();
  } 
};

template <class T_SCALAR>
class test_suite {
public:

  static void test_group_pick() {
    sl::random::irng_marsaglia rng;

    const sl::uint32_t N = 10;
    sl::uint32_t h[N];
    for (sl::uint32_t k=1; k<=N; ++k) {
      sl::uint32_t a[N];
      // Clear histogram
      for (sl::uint32_t i=0; i<N; ++i) {
	h[i] = 0;
      }
      sl::uint32_t picked_count = 0;
      for (sl::uint32_t m=0; m<100000; ++m) {
	rng.pick_k_out_of_n_in(k, a, N);
	for (sl::uint32_t j=0; j<k; ++j) {
	  ++h[a[j]];
	  ++picked_count;
	}
      }
      
      double expected = double(picked_count) / double(N);

      sl::tester tester(std::string()+"random pick " + sl::to_string(k) + " out of " + sl::to_string(N));
      for (sl::uint32_t i=0; i<N; ++i) {
	tester.test(std::string()+"histo bucket " + sl::to_string(i),
		    double(h[i])/double(expected),
		    sl::intervald(0.95,1.05));
      }
      failed_test_count += tester.failed_test_count();
    }
  }
  
  static void do_it() {
    test_random_uniform<T_SCALAR>::do_it();
    test_group_pick();
  }    

};

int main() {
  test_suite<float>::do_it();
  test_suite<double>::do_it();

  return (int)failed_test_count;
}


