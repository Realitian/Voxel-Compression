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
#include <sl/benchmarker.hpp>
#include <sl/bounded_scalar.hpp>
#include <sl/affine_map.hpp>
#include <stdlib.h> 
//#include <sstream>

//---------------------------------------------------------
// op_matrix_multiply
//---------------------------------------------------------

template <typename U_SCALAR>
class op_matrix_multiply: public sl::operation {
public:
  typedef U_SCALAR value_t;
  typedef sl::fixed_size_square_matrix<4,value_t> matrix_t;
protected:
  matrix_t m1;
  matrix_t m2;
  matrix_t m3;
  bool is_affine;

  op_matrix_multiply() {
  }

public:

  op_matrix_multiply(bool affine) {
    make_matrices(affine);
  }

  virtual std::string method_id() {
    std::string s;
    // OLD s << m1.row_count() << "x" << m1.column_count() << " CPP   operator* matrix multiply";
    return s;
  }

  virtual void make_matrices(bool affine) {
    is_affine = affine;
    srand(1234);
    for (size_t i=0; i< m1.row_count(); i++) {
      for (size_t j=0; j< m1.column_count(); j++) {
	m1(i,j) = value_t(rand()/(double)(RAND_MAX));
	m2(i,j) = value_t(rand()/(double)(RAND_MAX));
	m3(i,j) = value_t(rand()/(double)(RAND_MAX));
      }
    }
    if (affine) {
      size_t n = m1.column_count();
      for (size_t j=0; j< n; j++) {
	m1(n-1,j) = (j == n-1 ? value_t(1.0) : value_t(0.0));
	m2(n-1,j) = (j == n-1 ? value_t(1.0) : value_t(0.0));
      }
    }
  }
  
  virtual void operate() {
    m3 = m1 * m2;
  }
  
  virtual matrix_t result() const {
    return m3;
  }
};

//---------------------------------------------------------
// op_c_matrix_multiply
//---------------------------------------------------------

template <class U_SCALAR>
class op_c_matrix_multiply: public op_matrix_multiply<U_SCALAR> {
protected:
  typedef op_matrix_multiply<U_SCALAR> super_t;
  typedef typename super_t::value_t value_t;
  typedef value_t c_matrix_4x4_t[4][4];
  c_matrix_4x4_t *m1_ptr;
  c_matrix_4x4_t *m2_ptr;
  c_matrix_4x4_t *m3_ptr;
public:
  
  op_c_matrix_multiply(bool affine) {
    super_t::make_matrices(affine);
    m1_ptr = (c_matrix_4x4_t *)(this->m1).to_pointer();
    m2_ptr = (c_matrix_4x4_t *)(this->m2).to_pointer();
    m3_ptr = (c_matrix_4x4_t *)(this->m3).to_pointer();
  }
  
  virtual std::string method_id() {
    std::string s;
    // OLD s << m1.row_count() << "x" << m1.column_count() << " RAW C matrix multiply";
    return s;
  }

  virtual void operate() {
    // GENERAL * GENERAL 
    (*m3_ptr)[0][0] = (*m2_ptr)[0][0]*(*m1_ptr)[0][0]+(*m2_ptr)[0][1]*(*m1_ptr)[1][0]+(*m2_ptr)[0][2]*(*m1_ptr)[2][0]+(*m2_ptr)[0][3]*(*m1_ptr)[3][0];
    (*m3_ptr)[0][1] = (*m2_ptr)[0][0]*(*m1_ptr)[0][1]+(*m2_ptr)[0][1]*(*m1_ptr)[1][1]+(*m2_ptr)[0][2]*(*m1_ptr)[2][1]+(*m2_ptr)[0][3]*(*m1_ptr)[3][1];
    (*m3_ptr)[0][2] = (*m2_ptr)[0][0]*(*m1_ptr)[0][2]+(*m2_ptr)[0][1]*(*m1_ptr)[1][2]+(*m2_ptr)[0][2]*(*m1_ptr)[2][2]+(*m2_ptr)[0][3]*(*m1_ptr)[3][2];
    (*m3_ptr)[0][3] = (*m2_ptr)[0][0]*(*m1_ptr)[0][3]+(*m2_ptr)[0][1]*(*m1_ptr)[1][3]+(*m2_ptr)[0][2]*(*m1_ptr)[2][3]+(*m2_ptr)[0][3]*(*m1_ptr)[3][3];
    
    (*m3_ptr)[1][0] = (*m2_ptr)[1][0]*(*m1_ptr)[0][0]+(*m2_ptr)[1][1]*(*m1_ptr)[1][0]+(*m2_ptr)[1][2]*(*m1_ptr)[2][0]+(*m2_ptr)[1][3]*(*m1_ptr)[3][0];
    (*m3_ptr)[1][1] = (*m2_ptr)[1][0]*(*m1_ptr)[0][1]+(*m2_ptr)[1][1]*(*m1_ptr)[1][1]+(*m2_ptr)[1][2]*(*m1_ptr)[2][1]+(*m2_ptr)[1][3]*(*m1_ptr)[3][1];
    (*m3_ptr)[1][2] = (*m2_ptr)[1][0]*(*m1_ptr)[0][2]+(*m2_ptr)[1][1]*(*m1_ptr)[1][2]+(*m2_ptr)[1][2]*(*m1_ptr)[2][2]+(*m2_ptr)[1][3]*(*m1_ptr)[3][2];
    (*m3_ptr)[1][3] = (*m2_ptr)[1][0]*(*m1_ptr)[0][3]+(*m2_ptr)[1][1]*(*m1_ptr)[1][3]+(*m2_ptr)[1][2]*(*m1_ptr)[2][3]+(*m2_ptr)[1][3]*(*m1_ptr)[3][3];
    
    (*m3_ptr)[2][0] = (*m2_ptr)[2][0]*(*m1_ptr)[0][0]+(*m2_ptr)[2][1]*(*m1_ptr)[1][0]+(*m2_ptr)[2][2]*(*m1_ptr)[2][0]+(*m2_ptr)[2][3]*(*m1_ptr)[3][0];
    (*m3_ptr)[2][1] = (*m2_ptr)[2][0]*(*m1_ptr)[0][1]+(*m2_ptr)[2][1]*(*m1_ptr)[1][1]+(*m2_ptr)[2][2]*(*m1_ptr)[2][1]+(*m2_ptr)[2][3]*(*m1_ptr)[3][1];
    (*m3_ptr)[2][2] = (*m2_ptr)[2][0]*(*m1_ptr)[0][2]+(*m2_ptr)[2][1]*(*m1_ptr)[1][2]+(*m2_ptr)[2][2]*(*m1_ptr)[2][2]+(*m2_ptr)[2][3]*(*m1_ptr)[3][2];
    (*m3_ptr)[2][3] = (*m2_ptr)[2][0]*(*m1_ptr)[0][3]+(*m2_ptr)[2][1]*(*m1_ptr)[1][3]+(*m2_ptr)[2][2]*(*m1_ptr)[2][3]+(*m2_ptr)[2][3]*(*m1_ptr)[3][3];
    
    (*m3_ptr)[3][0] = (*m2_ptr)[3][0]*(*m1_ptr)[0][0]+(*m2_ptr)[3][1]*(*m1_ptr)[1][0]+(*m2_ptr)[3][2]*(*m1_ptr)[2][0]+(*m2_ptr)[3][3]*(*m1_ptr)[3][0];
    (*m3_ptr)[3][1] = (*m2_ptr)[3][0]*(*m1_ptr)[0][1]+(*m2_ptr)[3][1]*(*m1_ptr)[1][1]+(*m2_ptr)[3][2]*(*m1_ptr)[2][1]+(*m2_ptr)[3][3]*(*m1_ptr)[3][1];
    (*m3_ptr)[3][2] = (*m2_ptr)[3][0]*(*m1_ptr)[0][2]+(*m2_ptr)[3][1]*(*m1_ptr)[1][2]+(*m2_ptr)[3][2]*(*m1_ptr)[2][2]+(*m2_ptr)[3][3]*(*m1_ptr)[3][2];
    (*m3_ptr)[3][3] = (*m2_ptr)[3][0]*(*m1_ptr)[0][3]+(*m2_ptr)[3][1]*(*m1_ptr)[1][3]+(*m2_ptr)[3][2]*(*m1_ptr)[2][3]+(*m2_ptr)[3][3]*(*m1_ptr)[3][3]; 
  }
};

//---------------------------------------------------------
// op_cpp_explicit_matrix_multiply
//---------------------------------------------------------

template <class U_SCALAR>
class op_cpp_explicit_matrix_multiply: public op_matrix_multiply<U_SCALAR> {
protected:
  typedef op_matrix_multiply<U_SCALAR> super_t;
  typedef typename super_t::value_t value_t;

public:
  
  op_cpp_explicit_matrix_multiply(bool affine) {
    super_t::make_matrices(affine);
  }
  
  virtual std::string method_id() {
    std::string s;
    // OLD s << "" << m1.row_count() << "x" << m1.column_count() << " CPP item(i,j) matrix multiply";
    return s;
  }

  virtual void operate() {
    const size_t M = (this->m1).column_count();
    const size_t N = (this->m1).row_count();
    const size_t K = (this->m2).column_count();
    for(size_t i=0; i < K; i++) {
      for(size_t j=0; j < N; j++) {
	value_t sum = (this->m1)(j,0) * (this->m2)(0,i);
	for(size_t k=1; k < M; k++) {
	  sum += (this->m1)(j,k) * (this->m2)(k,i);
	}
	(this->m3)(j,i) = sum;
      }           
    }
  }
};

//---------------------------------------------------------
// benchmark_geometry
//---------------------------------------------------------

template <typename T_SCALAR>
class benchmark_geometry {
public:

  static double bench_and_report(op_matrix_multiply<T_SCALAR>* op,
				 const std::string& t_scalar_id,
				 double reference_rate) {
    sl::benchmarker b;

    b.measure_rate(*op);
    double result = b.last_cpu_time_rate();
    std::cerr << op->method_id() << " (" << t_scalar_id << "): ";
    std::cerr << result << " ops/s ";
    if (reference_rate > 0) {
      std::cerr << (result / reference_rate) << " wrt ref";
    } 
    std::cerr << std::endl;

    return result;
  }

  static void do_it(const std::string& t_scalar_id) {
    for (bool affine = false; affine != true; affine = !affine) {
      std::cerr << "GENERAL ------------------------------------------" << std::endl;
      {
	double raw_c_cpu_time   = 
	  bench_and_report(new op_c_matrix_multiply<T_SCALAR>(affine), 
			   t_scalar_id, 0.0);
	(void)
	  bench_and_report(new op_matrix_multiply<T_SCALAR>(affine), 
			   t_scalar_id, raw_c_cpu_time);
	(void)
	  bench_and_report(new op_cpp_explicit_matrix_multiply<T_SCALAR>(affine), 
			   t_scalar_id, raw_c_cpu_time);
      }
    }
  }

}; // 

int main() {
  benchmark_geometry<float>::do_it("float");
  benchmark_geometry<double>::do_it("double");
  //  benchmark_geometry< sl::bounded_scalar<float> >::do_it("sl::bounded_scalar<float>");
  //  benchmark_geometry< sl::bounded_scalar<double> >::do_it("sl::bounded_scalar<double>");

}
