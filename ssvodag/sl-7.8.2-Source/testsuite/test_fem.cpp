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
#include <sl/fem_basis.hpp>
#include <sl/linear_map_factory.hpp>
#include <sl/random.hpp>

//--------------------

static std::size_t failed_test_count = 0;

template <class G_scalar>
class test_fem_base {
public:
  typedef G_scalar                                             value_t;
  typedef sl::interval<G_scalar>                               interval_t;

  typedef sl::linear_map_factory<2,value_t>                    linear_map_factory_t;
  typedef typename linear_map_factory_t::affine_map_t          affine_map_t;
  typedef typename linear_map_factory_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename linear_map_factory_t::projective_map_t      projective_map_t;

  typedef typename affine_map_t::point_t                       point_t;
  typedef typename affine_map_t::vector_t                      vector_t;
  typedef typename affine_map_t::dual_vector_t                 dual_vector_t;
  typedef typename affine_map_t::plane_t                       plane_t;

  typedef sl::fem_basis<value_t,2>                             basis_t;
  typedef sl::dense_array<value_t,1,void>                      coeff_t;

public:

  static sl::random::uniform<value_t> urandom;
  static sl::quad01_fem_basis_factory<G_scalar>* fem_basis_factory;

public:

  static void random_reset() {
    urandom.set_seed(0);
  }

  static value_t random_01_value() {
    return urandom.value();
  }

  static point_t random_01_point() {
    return point_t(urandom.value(),
		   urandom.value(),
		   urandom.value());

  }

  static coeff_t random_01_coefficients(std::size_t n) {
    coeff_t result = coeff_t(n);
    for (std::size_t i=0; i<n; ++i) {
      result(i) = urandom.value();
    }
    return result;
  }

  static void do_it() {
    
  }
};


template <class G_scalar>
sl::random::uniform<G_scalar>  test_fem_base<G_scalar>::urandom;

template <class G_scalar>
sl::quad01_fem_basis_factory<G_scalar>* test_fem_base<G_scalar>::fem_basis_factory = new sl::quad01_fem_basis_factory<G_scalar>();

template <class G_scalar>
class test_fem_quad01: public test_fem_base<G_scalar> {
protected:
  typedef test_fem_base<G_scalar> super_t;
  typedef G_scalar                                             value_t;
  typedef sl::interval<G_scalar>                               interval_t;

  typedef sl::linear_map_factory<2,value_t>                    linear_map_factory_t;
  typedef typename linear_map_factory_t::affine_map_t          affine_map_t;
  typedef typename linear_map_factory_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename linear_map_factory_t::projective_map_t      projective_map_t;

  typedef typename affine_map_t::point_t                       point_t;
  typedef typename affine_map_t::vector_t                      vector_t;
  typedef typename affine_map_t::dual_vector_t                 dual_vector_t;
  typedef typename affine_map_t::plane_t                       plane_t;

  typedef sl::fem_basis<value_t,2>                             basis_t;
  typedef sl::dense_array<value_t,1,void>                      coeff_t;
public:

  static void test_pushpull(sl::tester& tester,
			    const basis_t& basis) {
    const std::size_t N  = basis.count();
    
    /// Test consistency of push-pull operations for a number of cases
    for (std::size_t i=0; i<10; ++i) {
      // Define parent-child relative geometry
      // Child 0: 0..u_sub, 0..1
      // Child 1: u_sub..1, 0..1
      const value_t u_sub = value_t(1.0)/value_t(2+i);

      const affine_map_t parent_to_parent;
      const affine_map_t child0_to_parent = 
	linear_map_factory_t::scaling(vector_t(u_sub, value_t(1.0)));
      const affine_map_t child1_to_parent = 
	linear_map_factory_t::translation(vector_t(u_sub, value_t(0.0))) *
	linear_map_factory_t::scaling(vector_t(value_t(1.0)-u_sub, value_t(1.0)));
      const value_t parent_area = value_t(1.0);
      const value_t child0_area = u_sub;
      const value_t child1_area = parent_area - child0_area;

      // Compute push-pull matrix
      sl::dense_array<value_t,2,void> self_h(N,N);
      sl::fem_push_pull_matrix_in(self_h,
				  basis,
				  basis,
				  parent_to_parent);
      sl::dense_array<value_t,2,void> child0_h(N,N);
      sl::fem_push_pull_matrix_in(child0_h,
				  basis,
				  basis,
				  child0_to_parent);
      sl::dense_array<value_t,2,void> child1_h(N,N);
      sl::fem_push_pull_matrix_in(child1_h,
				  basis,
				  basis,
				  child1_to_parent);
      
      // Self Push-Pull
      coeff_t parent_coefficients = super_t::random_01_coefficients(N);
      coeff_t parent_coefficients_prime = coeff_t(N);
      coeff_t child0_coefficients = coeff_t(N);
      coeff_t child1_coefficients = coeff_t(N);
      coeff_t parent2_coefficients = coeff_t(N);

      sl::fem_push_in(parent_coefficients_prime.as_sized_raw_array_pointer(),
		      parent_coefficients.as_const_sized_raw_array_pointer(),
		      self_h);
      tester.test("Self Push consistency", compare(parent_coefficients, parent_coefficients_prime, value_t(0.01)), 0);
      parent2_coefficients.clear();
      sl::fem_accumulate_pull_in(parent2_coefficients.as_sized_raw_array_pointer(),
				 self_h,
				 parent_coefficients_prime.as_const_sized_raw_array_pointer(),
				 parent_area / parent_area);
      tester.test("Self Pull consistency", compare(parent_coefficients, parent2_coefficients, value_t(0.01)), 0);
      
      SL_TRACE_OUT(2) << "SELF-H  = " << std::endl << self_h << std::endl;
      SL_TRACE_OUT(2) << "PARENT-COEFF  = " << std::endl << parent_coefficients << std::endl;
      SL_TRACE_OUT(2) << "PUSHED-COEFF  = " << std::endl << parent_coefficients_prime << std::endl;
      SL_TRACE_OUT(2) << "PULLED-COEFF = " << std::endl << parent2_coefficients << std::endl;
      
      // Push
      sl::fem_push_in(child0_coefficients.as_sized_raw_array_pointer(),
		      parent_coefficients.as_const_sized_raw_array_pointer(),
		      child0_h);
      sl::fem_push_in(child1_coefficients.as_sized_raw_array_pointer(),
		      parent_coefficients.as_const_sized_raw_array_pointer(),
		      child1_h);


      // Pull what we pushed - should get the original values
      parent2_coefficients.clear();
      sl::fem_accumulate_pull_in(parent2_coefficients.as_sized_raw_array_pointer(),
				 child0_h,
				 child0_coefficients.as_const_sized_raw_array_pointer(),
				 child0_area / parent_area);
      sl::fem_accumulate_pull_in(parent2_coefficients.as_sized_raw_array_pointer(),
				 child1_h,
				 child1_coefficients.as_const_sized_raw_array_pointer(),
				 child1_area / parent_area);
      
      tester.test("Push-pull consistency", compare(parent_coefficients, parent2_coefficients, value_t(0.01)), 0);

      SL_TRACE_OUT(2) << "CHILD1-H  = " << std::endl << child0_h << std::endl;
      SL_TRACE_OUT(2) << "CHILD2-H  = " << std::endl << child1_h << std::endl;
      SL_TRACE_OUT(2) << "PARENT-COEFF  = " << std::endl << parent_coefficients << std::endl;
      SL_TRACE_OUT(2) << "CHILD1-COEFF  = " << std::endl << child0_coefficients << std::endl;
      SL_TRACE_OUT(2) << "CHILD2-COEFF  = " << std::endl << child1_coefficients << std::endl;
      SL_TRACE_OUT(2) << "PARENT2-COEFF = " << std::endl << parent2_coefficients << std::endl;
      
      // TODO: Eval consistency
      tester.test("Eval consistency - child0", basis.value_at(point_t(0.5f,0.5f), child0_coefficients.as_const_sized_raw_array_pointer()), 
		  interval_t(basis.value_at(transformation(child0_to_parent, point_t(0.5f,0.5f)), parent_coefficients.as_const_sized_raw_array_pointer())-value_t(0.001),
			     basis.value_at(transformation(child0_to_parent, point_t(0.5f,0.5f)), parent_coefficients.as_const_sized_raw_array_pointer())+value_t(0.001)));
      tester.test("Eval consistency - child1", basis.value_at(point_t(0.5f,0.5f), child1_coefficients.as_const_sized_raw_array_pointer()), 
		  interval_t(basis.value_at(transformation(child1_to_parent, point_t(0.5f,0.5f)), parent_coefficients.as_const_sized_raw_array_pointer())-value_t(0.001),
			     basis.value_at(transformation(child1_to_parent, point_t(0.5f,0.5f)), parent_coefficients.as_const_sized_raw_array_pointer())+value_t(0.001)));

      // TODO: Generalized push pull
    }
  }

  static void test_overlap(sl::tester& tester,
			    const basis_t& basis) {
    const std::size_t N  = basis.count();
    
    for (std::size_t i=0; i<N; ++i) {
      for (std::size_t j=0; j<N; ++j) {
	const sl::cubature_rule<value_t,2>& cr = basis.get_cubature_rule_factory().rule_from_degree(basis(i).maximum_degree() + basis(j).maximum_degree());
	value_t O_i_j = sl::scalar_math<value_t>::zero();
	for (std::size_t k=0; k<cr.node_count(); ++k) {
	  O_i_j += cr.weight(k) * basis(i)(cr.position(k)) * basis(j)(cr.position(k));
	}
	interval_t expected = interval_t(value_t(-0.001), value_t(0.001));
	if (i==j) expected = interval_t(value_t(1.0-0.001), value_t(1.0+0.001));

	tester.test("Overlap " + sl::to_string(i) + "-" + sl::to_string(j),
		    O_i_j,
		    expected);
      }
    }
  }
  
  static void do_it(sl::fem_basis_type tp) {
    sl::tester tester("fem_basis < " + sl::numeric_traits<value_t>::what() + " > for unit quad of type " + sl::to_string(tp));
    const basis_t& basis = super_t::fem_basis_factory->basis(tp);

    test_overlap(tester, basis);
    test_pushpull(tester, basis);

    failed_test_count += tester.failed_test_count();
  }
  

  static void do_it() {
    union { 
      sl::fem_basis_type tp;
      int val;
    } tp;

    for (tp.tp = sl::FEM_BASIS_FIRST_TYPE; 
	 tp.tp != sl::FEM_BASIS_LAST_TYPE;
	 ++tp.val) {
      do_it(tp.tp);
    }
  }

};

template <class G_scalar>
class test_suite {
public:

  static void do_it() {
    test_fem_quad01<G_scalar>::do_it();
  }    

};

int main() {

  test_suite<float>::do_it();

  return (int)failed_test_count;

}


