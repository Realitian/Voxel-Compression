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
#include <sl/cubature_rule.hpp>
#include <sl/linear_map_factory.hpp>
#include <sl/random.hpp>
#include <iostream>
#include <iomanip>

//--------------------

static std::size_t failed_test_count = 0;

template <class G_scalar, std::size_t G_dimension>
class monomial {
  enum { dimension = G_dimension };

  typedef G_scalar                                  value_t;
  typedef sl::interval<value_t>                     interval_t;
  typedef sl::fixed_size_point<dimension, value_t>  point_t;
  typedef sl::fixed_size_array<dimension, int>      exponent_t;

protected:
  
  exponent_t n_;

public:

  inline value_t brute_force_pow(value_t x, int n) const {
    value_t pow = value_t(1.);
    while (n>0) {
      pow *= x;
      n--;
    }
    return pow;
  }

  inline monomial(const exponent_t& n): n_(n) {};

  inline value_t operator()(const point_t& x) const {
    SL_REQUIRE("Check power", sl::abs( sl::ipow(x[0],n_[0]) -  brute_force_pow(x[0],n_[0]) ) < 0.0001);
    value_t result = sl::ipow(x[0],n_[0]);
    for (std::size_t i=1; i<dimension; ++i) {
      result *= sl::ipow(x[i],n_[i]);
    }
    return result;
  }

}; // monomial

template <class G_scalar>
class test_cubature {
public:
  typedef G_scalar                                             value_t;
  typedef sl::interval<G_scalar>                               interval_t;

  typedef sl::fixed_size_point<2, value_t>                     point_t;
  typedef sl::fixed_size_array<2, int>                         exponent_t;
  typedef monomial<value_t, 2>                                 monomial_t;

  typedef sl::cubature_rule<value_t, 2>                        cubarule_t;

public:

  static value_t brute_force_integral(const monomial_t& f) {
    const std::size_t N = 100;
    value_t result = value_t(0.0);
    for (std::size_t i=0; i<N; ++i) {
      value_t u = value_t(i+0.5)/value_t(N);
      for (std::size_t j=0; j<N; ++j) {
	value_t v = value_t(j+0.5)/value_t(N);

	result += f(point_t(u,v));
      }
    }
    result /= value_t(N*N);
    return result;
  }   

  static void test_monomial(const cubarule_t& r) {
    sl::tester tester("Quad cubature rule <" + 
		      sl::numeric_traits<value_t>::what() + 
		      " > of degree " + sl::to_string(r.degree()) + 
		      " with " + sl::to_string(r.node_count()) + " nodes");
    
    std::size_t max_degree = r.degree();
    for (std::size_t i=0; i<=max_degree; ++i) {
      for (std::size_t j=0; j<=max_degree; ++j) {
	if (i+j <= max_degree) {
	  exponent_t e; e = i, j;

	  tester.test("Monomial u^" + sl::to_string(i) + "*v^" + sl::to_string(j),
		      sl::integral<value_t>(r,monomial_t(e)),
		      interval_t(brute_force_integral(monomial_t(e))-value_t(1.0e-3),
				 brute_force_integral(monomial_t(e))+value_t(1.0e-3)));
	}
      }
    }

    failed_test_count += tester.failed_test_count();
  }

  static void do_it() {
    sl::quad01_cubature_rule_factory<value_t> cf;
    for (std::size_t i=cf.minimum_degree(); i<=cf.maximum_degree(); ++i) {
      test_monomial(cf.rule_from_degree(i));
    }
  }

};

int main() {

  test_cubature<float>::do_it();

  return (int)failed_test_count;

}


