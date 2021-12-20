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
#include <sl/interval.hpp>
#include <sl/bounded_scalar.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/fixed_size_vector.hpp>
#include <sl/fixed_size_packed_matrix.hpp>


#include <iostream>

static std::size_t failed_test_count = 0;

static void test_array() {
  
  sl::fixed_size_array<5, int> a1; a1 = 1, 2, 3, 4, 5;
  sl::fixed_size_array<5, int> a2; a2 = 1;

  std::cerr << "a1 = " << a1 << std::endl;
  std::cerr << "a2 = " << a2 << std::endl;
  
}


template <typename T_SCALAR>
class test_matrix {
public:
  typedef sl::interval<T_SCALAR> interval_t;

  static void do_it(const std::string& t_scalar_id) {
    {
      sl::tester tester("Matrix2 <" + t_scalar_id + ">");
      
      sl::fixed_size_square_matrix<2,T_SCALAR> a(T_SCALAR(3.),T_SCALAR(0.),
					    T_SCALAR(0.),T_SCALAR(4.));
      
      tester.test("a.two_norm()", a.two_norm(), interval_t("5.0000"));
      tester.test("a * ~a", ((a * ~a) - sl::fixed_size_square_matrix_factory<2,T_SCALAR>::identity()).amax(), interval_t("0.0"));

      failed_test_count += tester.failed_test_count();
    }
    
    {
      sl::tester tester("Matrix4 <" + t_scalar_id + ">");
      
      sl::fixed_size_square_matrix<4,T_SCALAR> a(2.0f,0.0f,2.0f,1.0f,
					    1.0f,2.0f,2.0f,4.0f,
					    2.0f,2.0f,3.0f,2.0f,
					    5.0f,2.0f,2.0f,3.0f);
      
      tester.test("a.two_norm()", a.two_norm(), interval_t("9.85"));
      tester.test("a * ~a", ((a * ~a) - sl::fixed_size_square_matrix_factory<4,T_SCALAR>::identity()).amax(), interval_t("0.00"));

      sl::fixed_size_square_matrix<4,T_SCALAR> a_aff(1.1f,2.3f,3.1f,4.1f,
						5.2f,6.5f,7.3f,8.2f,
						9.3f,10.7f,11.2f,12.1f,
						0.0f,0.0f,0.0f,1.0f);

      tester.test("a_aff * ~a_aff", ((a_aff * ~a_aff) - sl::fixed_size_square_matrix_factory<4,T_SCALAR>::identity()).amax(), interval_t("0.00"));

      failed_test_count += tester.failed_test_count();
    }
    
    {
      sl::tester tester("packed_symmetric_matrix4 <" + t_scalar_id + ">");
      
      sl::fixed_size_packed_symmetric_matrix<4,T_SCALAR> a;
      a(0,0) = 1.0f; a(0,1)=2.0f; a(0,2)=3.0f; a(0,3)=4.0f;
      /*          */ a(1,1)=5.0f; a(1,2)=6.0f; a(1,3)=7.0f;
      /*                       */ a(2,2)=8.0f; a(2,3)=9.0f;
      /*                                    */ a(3,3)=10.0f;

      tester.test("a(3,1) == a(1,3)", a(3,1) == a(1,3));
      tester.test("a(2,3) == a(3,2)", a(2,3) == a(3,2));
      tester.test("a * ~a", ((a * ~a) - sl::fixed_size_packed_symmetric_matrix_factory<4,T_SCALAR>::identity()).amax(), interval_t("0.00"));
      failed_test_count += tester.failed_test_count();
    }

  }
};


int main() {
  test_array();

  test_matrix<float>::do_it("float");
  test_matrix<double>::do_it("double");
  test_matrix< sl::bounded_scalar<float> >::do_it("sl::bounded_scalar<float>");
  test_matrix< sl::bounded_scalar<double> >::do_it("sl::bounded_scalar<double>");

  return (int)failed_test_count;
}



