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

# ifdef _MSC_VER
#   pragma warning(disable:4244)
# endif //_MSC_VER

/////// ALWAYS TEST IN DEBUG MODE
#if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////

#include <sl/tester.hpp>
#include <sl/interval.hpp>
#include <sl/bounded_scalar.hpp>
#include <sl/affine_map.hpp>
#include <sl/projective_map.hpp>
#include <sl/cartesian_frame.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/oriented_box.hpp>
#include <sl/polygon.hpp>
#include <sl/convex_hull.hpp>
#include <sl/shaft.hpp>
#include <sl/rigid_body_map.hpp>
#include <sl/linear_map_factory.hpp>
#include <sl/fixed_size_quadric_matrix.hpp>
#include <sl/minimal_area_triangulator.hpp>
#include <sl/random.hpp>

//--------------------

static std::size_t failed_test_count = 0;

template <class T_SCALAR>
class test_geometry_base {
public:
  typedef T_SCALAR                           value_t;
  typedef sl::interval<T_SCALAR>                      interval_t;

  typedef sl::linear_map_factory<3,T_SCALAR>                   linear_map_factory_t;
  typedef typename linear_map_factory_t::affine_map_t          affine_map_t;
  typedef typename linear_map_factory_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename linear_map_factory_t::projective_map_t      projective_map_t;

  typedef typename affine_map_t::point_t              point_t;
  typedef typename affine_map_t::vector_t             vector_t;
  typedef typename affine_map_t::dual_vector_t        dual_vector_t;
  typedef typename affine_map_t::quaternion_t         quaternion_t;
  typedef typename affine_map_t::plane_t              plane_t;

  static sl::random::uniform<value_t> urandom;

  static void random_reset() {
    urandom.set_seed(0);
  }

  static value_t random_unit_value() {
    return urandom.value();
  }

  static point_t random_unit_point() {
    return point_t(-1.0f + 2.0f * urandom.value(),
		   -1.0f + 2.0f * urandom.value(),
		   -1.0f + 2.0f * urandom.value());

  }

  static vector_t random_unit_vector() {
    return vector_t(-1.0f + 2.0f * urandom.value(),
		    -1.0f + 2.0f * urandom.value(),
		    -1.0f + 2.0f * urandom.value()).ok_normalized();
		    
  }

  static dual_vector_t random_unit_dual_vector() {
    return dual_vector_t(-1.0f + 2.0f * urandom.value(),
			 -1.0f + 2.0f * urandom.value(),
			 -1.0f + 2.0f * urandom.value()).ok_normalized();
		    
  }

  static quaternion_t random_quaternion() {
    return quaternion_t(random_unit_vector(), urandom.value() * sl::scalar_math<value_t>::Pi());
  }

  static plane_t random_plane() {
    return plane_t(random_unit_dual_vector(), random_unit_point());
  }

  static void do_it() {
    
  }
};


template <class T_SCALAR>
sl::random::uniform<T_SCALAR>  test_geometry_base<T_SCALAR>::urandom;

template <class T_SCALAR>
class test_point: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:
  static void do_it() {
    sl::tester tester("point < " + sl::numeric_traits<value_t>::what() + " >");

    {
      point_t P0; P0 = value_t(0.0f), value_t(0.0f), value_t(0.0f);
      point_t P1; P1 = value_t(0.0f), value_t(0.0f), value_t(1.0f);
      point_t P2; P2 = value_t(0.0f), value_t(1.0f), value_t(1.0f);
      point_t P;  P  = value_t(0.0f), value_t(1.0f), value_t(0.0f);
      point_t Q;  Q  = value_t(1.0f), value_t(1.0f), value_t(0.0f);
    
      tester.test("(0,0,0).distance_squared_to((0,0,0),(0,0,1),(0,1,1)", P0.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(0,0,1).distance_squared_to((0,0,0),(0,0,1),(0,1,1)", P1.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(0,1,1).distance_squared_to((0,0,0),(0,0,1),(0,1,1)", P2.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(0,1,0).distance_squared_to((0,0,0),(0,0,1),(0,1,1)",  P.distance_squared_to(P0,P1,P2), interval_t("0.500"));
      tester.test("(1,1,0).distance_squared_to((0,0,0),(0,0,1),(0,1,1)",  Q.distance_squared_to(P0,P1,P2), interval_t("1.500"));
    }

    {
      point_t P0; P0 = value_t(23.0f), value_t(3.0f), value_t(1.0f);
      point_t P1; P1 = value_t(25.0f), value_t(5.0f), value_t(1.0f);
      point_t P2; P2 = value_t(22.0f), value_t(5.0f), value_t(1.0f);
    
      tester.test("(23,3,1).distance_squared_to((23,3,1),(25,5,1),(22,5,1)", P0.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(25,5,1).distance_squared_to((23,3,1),(25,5,1),(22,5,1)", P1.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(22,5,1).distance_squared_to((23,3,1),(25,5,1),(22,5,1)", P2.distance_squared_to(P0,P1,P2), interval_t("0.000"));
    }
    {
      point_t P0; P0 = value_t(0.0f), value_t(3.0f), value_t(1.0f);
      point_t P1; P1 = value_t(2.0f), value_t(4.0f), value_t(1.0f);
      point_t P2; P2 = value_t(4.0f), value_t(5.0f), value_t(1.0f);
    
      tester.test("(0,3,1).distance_squared_to((0,3,1),(2,4,1),(4,5,1)", P0.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(2,4,1).distance_squared_to((0,3,1),(2,4,1),(4,5,1)", P1.distance_squared_to(P0,P1,P2), interval_t("0.000"));
      tester.test("(4,5,1).distance_squared_to((0,3,1),(2,4,1),(4,5,1)", P2.distance_squared_to(P0,P1,P2), interval_t("0.000"));
    }

    
    failed_test_count += tester.failed_test_count();
  } 
};

template <class T_SCALAR>
class test_plane: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:
  static void do_it() {
    sl::tester tester("plane < " + sl::numeric_traits<value_t>::what() + " >");
    
    plane_t plane(dual_vector_t(1.0,0.0,0.0), value_t(-2.0)); // "x == 2"
    
    tester.test("plane.normal()", (plane.normal() - dual_vector_t(1.0f,0.0f,0.0f)).two_norm(), interval_t("0.0"));
    tester.test("plane.value(2, 0, 0)", plane.value(point_t(2.0f,0.0f,0.0f)), interval_t("0.0"));
    tester.test("plane.value(3, 0, 0)", plane.value(point_t(3.0f,0.0f,0.0f)), interval_t("1.0"));
    tester.test("plane.value(1, 0, 0)", plane.value(point_t(1.0f,0.0f,0.0f)), interval_t("-1.0"));

    failed_test_count += tester.failed_test_count();
  } 
};

template <class T_SCALAR>
class test_quadric: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;

  typedef sl::fixed_size_quadric_matrix<3,value_t>      quadric_matrix_t;
  
public:
  
  static void do_it() {
    sl::tester tester("quadric matrix < " + sl::numeric_traits<value_t>::what() + " >");
    
    quadric_matrix_t Qz(point_t(0.0f,0.0f,0.0f),
                        point_t(0.0f,1.0f,0.0f),
                        point_t(1.0f,1.0f,0.0f));
    
    quadric_matrix_t Qx(point_t(0.0f,0.0f,0.0f),
                        point_t(0.0f,1.0f,0.0f),
                        point_t(0.0f,1.0f,1.0f));
    
    quadric_matrix_t Qy(point_t(0.0f,0.0f,0.0f),
                        point_t(0.0f,0.0f,1.0f),
                        point_t(1.0f,0.0f,1.0f));

    tester.test("Qz.evaluate(0.0,0.0,0.0)", Qz.evaluate(point_t(0.0f,0.0f,0.0f)), interval_t("0.000"));
    tester.test("Qz.evaluate(1.0,0.0,0.0)", Qz.evaluate(point_t(1.0f,0.0f,0.0f)), interval_t("0.000"));
    tester.test("Qz.evaluate(0.0,1.0,0.0)", Qz.evaluate(point_t(0.0f,1.0f,0.0f)), interval_t("0.000"));
    tester.test("Qz.evaluate(0.0,0.0,1.0)", Qz.evaluate(point_t(0.0f,0.0f,1.0f)), interval_t("1.000"));

    quadric_matrix_t Qz2 = quadric_matrix_t(plane_t(dual_vector_t(0.0f,0.0f,1.0f), value_t(1.0f)));

    tester.test("Qz2.evaluate(0.0,0.0,0.0)", Qz2.evaluate(point_t(0.0f,0.0f,0.0f)), interval_t("1.000"));
    tester.test("Qz2.evaluate(1.0,0.0,0.0)", Qz2.evaluate(point_t(1.0f,0.0f,0.0f)), interval_t("1.000"));
    tester.test("Qz2.evaluate(0.0,1.0,0.0)", Qz2.evaluate(point_t(0.0f,1.0f,0.0f)), interval_t("1.000"));
    tester.test("Qz2.evaluate(0.0,0.0,1.0)", Qz2.evaluate(point_t(0.0f,0.0f,1.0f)), interval_t("0.000"));
    
    quadric_matrix_t Qxyz = Qx+Qy+Qz;

    tester.test("Qxyz.evaluate(0.0,0.0,0.0)", Qz.evaluate(point_t(0.0f,0.0f,0.0f)), interval_t("0.000"));
    tester.test("Qxyz.evaluate(1.0,0.0,0.0)", Qz.evaluate(point_t(1.0f,0.0f,0.0f)), interval_t("0.000"));
    tester.test("Qxyz.evaluate(0.0,1.0,0.0)", Qxyz.evaluate(point_t(0.0f,1.0f,0.0f)), interval_t("1.000"));
    tester.test("Qxyz.evaluate(0.0,0.0,1.0)", Qxyz.evaluate(point_t(0.0f,0.0f,1.0f)), interval_t("1.000"));

    tester.test("Qxyz.evaluate(0.0,0.0,0.0)", Qxyz.evaluate(point_t(0.0f,0.0f,0.0f)), interval_t("0.000"));

    tester.test("Qxyz.optimized()", (Qxyz.optimized() - point_t(0.0f,0.0f,0.0f)).two_norm(), interval_t("0.000"));

    point_t P       = point_t(0.0f, 0.0f, 0.0f);
    point_t Q       = point_t(0.0f, 1.0f, 0.0f);
    point_t R       = point_t(1.0f, 1.0f, 0.0f);
    point_t P_prime = point_t(0.0f, 0.0f, 1.0f);

    value_t V_tetra = (1.0f/6.0f)*(P_prime-P).dot(sl::cross(P,Q,R));
    value_t V_tetra_squared = V_tetra*V_tetra;

    value_t A_tri_squared = sl::triangle_area_squared(P,Q,R);
    value_t Q_weight = A_tri_squared/9.0f;
    quadric_matrix_t Q_tri = quadric_matrix_t(P,Q,R);
    Q_tri *= Q_weight;

    tester.test("Q_tri", sl::abs(Q_tri.evaluate(P_prime)-V_tetra_squared), interval_t("0.0000"));

    typedef sl::fixed_size_point<6,value_t> point6_t;
    typedef sl::fixed_size_quadric_matrix<6,value_t> quadric_matrix6_t;

    point6_t P6; P6 = 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
    point6_t Q6; Q6 = 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f;
    point6_t R6; R6 = 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f;
    point6_t P6_prime; P6_prime = 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f;

    quadric_matrix6_t Q6_tri = quadric_matrix6_t(P6,Q6,R6);
    tester.test("Q6", Q6_tri.evaluate(P6_prime), interval_t("1.000"));
    
    failed_test_count += tester.failed_test_count();
  } 
};

template <class T_SCALAR>
class test_quaternion: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:
  static void do_it() {
    sl::tester tester("quaternion < " + sl::numeric_traits<value_t>::what() + " >");
    
    quaternion_t a(vector_t(0,0,1), value_t(0.2)*sl::Pi(value_t()));
    
    tester.test("a.two_norm()", a.two_norm(), interval_t("1.0"));
    tester.test("a.axis()", (a.axis() - vector_t(0,0,1)).two_norm(), interval_t("0.0"));
    tester.test("a.angle()", sl::abs(a.angle()-value_t(0.2)*sl::Pi(value_t())), interval_t("0.0"));
    // EGO REMOVED    tester.test("a * ~a", ((a * ~a) - a.identity()).two_norm(), interval_t("0.0"));

    failed_test_count += tester.failed_test_count();
  } 
};


template <class T_SCALAR>
class test_affine_map: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:
  static void do_it() {
    sl::tester tester("affine_map < 3, " + sl::numeric_traits<value_t>::what() + " >");

    affine_map_t map = linear_map_factory_t::translation(10.0f,20.0f,30.0f);
    
    plane_t hp = plane_t(dual_vector_t(1.0f, 0.0f, 0.0f), point_t(2.0f, 3.0f, 4.0f));

    tester.test("map * point",  (map * point_t(1.0f,2.0f,3.0f) - point_t(11.0f,22.0f,33.0f)).two_norm(), interval_t("0.000"));
    tester.test("map * vector", (map * vector_t(1.0f,2.0f,3.0f) - vector_t(1.0f,2.0f,3.0f)).two_norm(), interval_t("0.000"));
    tester.test("map * plane",  ((hp.value(point_t(1.0f, 2.0f, 3.0f)) -
				  (map * hp).value(map * point_t(1.0f, 2.0f, 3.0f)))), interval_t("0.000"));

    tester.test("map * ~map",   ((map * ~map).as_matrix() - map.identity().as_matrix()).two_norm(), interval_t("0.00"));

    {
      const size_t N_TEST = 20;

      for (size_t i=0; i<N_TEST; ++i) {
	vector_t      t_0 = vector_t(10.1*i,20.2*i,30.3*i);
	vector_t      sh_0 = vector_t(0.0,0.0,0.0);
	vector_t      sc_0 = vector_t(1.1*i+1.0f,1.1*i+1.0f,1.1*i+1.0f);
	quaternion_t  r_0 = quaternion_t(vector_t(0,0,1),value_t(0.17*i));

	plane_t hp = plane_t(dual_vector_t::unit(i%3), point_t(2.*i, -3.*i, 4.*i));
	point_t p0 = point_t(1.0*i, -2.0*i, 3.0*i);
	point_t p1 = point_t(-1.0*i, -1.1*i, 6.1*i);
     
	affine_map_t map = 
	  linear_map_factory_t::translation(t_0) *
	  linear_map_factory_t::rotation(r_0) *
	  linear_map_factory_t::scaling(sc_0);
        
	tester.test("affine map * ~map",   ((map * ~map).as_matrix() - map.identity().as_matrix()).infinite_norm(), interval_t("0.00"));
	
	tester.test("affine map * vector vs. map * point",
		    (map * (p1 - p0) - (map * p1 - map * p0)).infinite_norm(), 
		    interval_t("0.000"));
      	tester.test("affine inverse xform point",
		    (inverse_transformation(map, p0) - transformation(~map, p0)).infinite_norm(), 
		    interval_t("0.000"));
      	tester.test("affine inverse xform vector",
		    (inverse_transformation(map, p0-p1) - transformation(~map, p0-p1)).infinite_norm(), 
		    interval_t("0.000"));
	tester.test("affine inverse xform plane",
		    (inverse_transformation(map, hp).as_vector() - transformation(~map, hp).as_vector()).infinite_norm(), 
		    interval_t("0.000"));

	vector_t      t;
	quaternion_t  r;
	vector_t      sh;
	vector_t      sc;
      
	map.factorize_to(sc,sh,r,t);

	tester.test("affine map -> translation",   (t - t_0).infinite_norm(), interval_t("0.000"));
	tester.test("affine map -> rotation",   (r - r_0).infinite_norm(), interval_t("0.000"));
	tester.test("affine map -> shear",   (sh - sh_0).infinite_norm(), interval_t("0.000"));
	tester.test("affine map -> scale",   (sc - sc_0).infinite_norm(), interval_t("0.000"));
	
      }
    }

    failed_test_count += tester.failed_test_count();
  } 
}; 

template <class T_SCALAR>
class test_rigid_body_map: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:
  static void do_it() {
    sl::tester tester("rigid_body_map < 3, " + sl::numeric_traits<value_t>::what() + " >");

    rigid_body_map_t map = linear_map_factory_t::translation(10.0f,20.0f,30.0f);
    
    plane_t hp = plane_t(dual_vector_t(1.0f, 0.0f, 0.0f), point_t(2.0f, 3.0f, 4.0f));

    tester.test("rigid map * point",  (map * point_t(1.0f,2.0f,3.0f) - point_t(11.0f,22.0f,33.0f)).two_norm(), interval_t("0.000"));
    tester.test("rigid map * vector", (map * vector_t(1.0f,2.0f,3.0f) - vector_t(1.0f,2.0f,3.0f)).two_norm(), interval_t("0.000"));
    tester.test("rigid map * plane",  ((hp.value(point_t(1.0f, 2.0f, 3.0f)) -
				  (map * hp).value(map * point_t(1.0f, 2.0f, 3.0f)))), interval_t("0.000"));
    {
      const size_t N_TEST = 20;

      super_t::random_reset();
      for (size_t i=0; i<N_TEST; ++i) {
	vector_t      t_0 = super_t::random_unit_vector();
	quaternion_t  r_0 = super_t::random_quaternion();

	plane_t hp = super_t::random_plane();
	point_t p0 = super_t::random_unit_point();
	point_t p1 = super_t::random_unit_point();
     
	rigid_body_map_t map = 
	  linear_map_factory_t::translation(t_0) *
	  linear_map_factory_t::rotation(r_0);

	tester.test("rigid map * ~map",   ((map * ~map).as_matrix() - map.identity().as_matrix()).infinite_norm(), interval_t("0.00"));
	
	tester.test("rigid map * vector vs. map * point",
		    ((map * (p1 - p0)) - (map * p1 - map * p0)).infinite_norm(), 
		    interval_t("0.000"));
	tester.test("rigid map * point preserves distances", 
		    (map * p1 - map* p0).two_norm() - (p1 - p0).two_norm(), 
		    interval_t("0.000"));
	tester.test("rigid map * plane preserves distances",  
		    hp.value(p0) - (map * hp).value(map * p0), 
		    interval_t("0.000"));
      	tester.test("rigid inverse xform point",
		    (inverse_transformation(map, p0) - transformation(~map, p0)).infinite_norm(), 
		    interval_t("0.000"));
      	tester.test("rigid inverse xform vector",
		    (inverse_transformation(map, p0-p1) - transformation(~map, p0-p1)).infinite_norm(), 
		    interval_t("0.000"));
	tester.test("rigid inverse xform plane",
		    (inverse_transformation(map, hp).as_vector() - transformation(~map, hp).as_vector()).infinite_norm(), 
		    interval_t("0.000"));

	vector_t      t;
	quaternion_t  r;
      
	map.factorize_to(r,t);

	tester.test("rigid map -> translation",   (t - t_0).infinite_norm(), interval_t("0.000"));
	tester.test("rigid map -> rotation",   
		    interval_t("0.00").contains((r * ~r_0).angle()) ||
		    interval_t("6.28").contains((r * ~r_0).angle()));

	tester.test("vector to vector mapper general", 
		    (linear_map_factory_t::rotation(t_0, vector_t::unit(2)) * t_0 -
		     vector_t::unit(2)).infinite_norm(), 
		    interval_t("0.000"));
	tester.test("vector to vector mapper parallel", 
		    (linear_map_factory_t::rotation(vector_t::unit(i%3), vector_t::unit(i%3)) * vector_t::unit(i%3) -
		     vector_t::unit(i%3)).infinite_norm(), 
		    interval_t("0.000"));
	tester.test("vector to vector mapper antiparallel", 
		    (linear_map_factory_t::rotation(-vector_t::unit(i%3), vector_t::unit(i%3)) * (-vector_t::unit(i%3)) -
		     vector_t::unit(i%3)).infinite_norm(), 
		    interval_t("0.000"));
      }
    } 

    failed_test_count += tester.failed_test_count();
  }
}; 

template <class T_SCALAR>
class test_projective_map: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static void do_it() {  
    sl::tester tester("projective_map < 3, " + sl::numeric_traits<value_t>::what() + " >");
    {
      projective_map_t map = linear_map_factory_t::perspective(value_t(30)*sl::Pi(value_t())/value_t(180.0),
							       value_t(1.0),
							       value_t(1.0),
							       value_t(10.0));
      
      tester.test("map * point",  (map * point_t(0.0f,0.0f,-1.0f) - point_t(0.0f,0.0f,-1.0)).two_norm(), interval_t("0.000"));
      tester.test("map * point",  (map * point_t(0.0f,0.0f,-10.0f) - point_t(0.0f,0.0f,1.0)).two_norm(), interval_t("0.000"));
      // EGO-REMOVED      tester.test("map * ~map",   ((map * ~map).as_matrix() - map.identity().as_matrix()).two_norm(), interval_t("0.00"));

      failed_test_count += tester.failed_test_count();
    }
    
  } // projective_map test
};


template <class T_SCALAR>
class test_cartesian_frame: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static void do_it() {  
    sl::tester tester("cartesian_frame < 3, " + sl::numeric_traits<value_t>::what() + " >");

    sl::cartesian_frame<3,value_t> frame = sl::cartesian_frame<3,value_t>(linear_map_factory_t::translation(10.0f,20.0f,30.0f));
    
    tester.test("frame * point",  (frame * point_t(1.0f,2.0f,3.0f) - point_t(11.0f,22.0f,33.0f)).two_norm(), interval_t("0.000"));
    tester.test("frame * vector", (frame * vector_t(1.0f,2.0f,3.0f) - vector_t(1.0f,2.0f,3.0f)).two_norm(), interval_t("0.000"));
    tester.test("frame * normal", (frame * vector_t(1.0f,2.0f,3.0f) - vector_t(1.0f,2.0f,3.0f)).two_norm(), interval_t("0.000"));

    failed_test_count += tester.failed_test_count();
  } 
};

template <class T_SCALAR>
class test_axis_aligned_box: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static void do_it() {  
    sl::tester tester("axis_aligned_box < 3, " + sl::numeric_traits<value_t>::what() + " >");

    sl::axis_aligned_box<3,value_t> box;
    point_t p0 = point_t(-1.0f,-2.0f,-3.0f);
    point_t p1 = point_t(1.0f,2.0f,3.0f);
    
    box.to(p0);
    box.merge(p1);
    
    tester.test("box diagonal", (box.diagonal() - (p1-p0)).infinite_norm(), interval_t("0.000"));
    tester.test("box contains", box.contains(p0.lerp(p1, value_t( 0.01))), true);
    tester.test("box contains", box.contains(p0.lerp(p1, value_t( 0.99f))), true);
    tester.test("box contains", box.contains(p0.lerp(p1, value_t( 1.01))), false);
    tester.test("box contains", box.contains(p0.lerp(p1, value_t(-0.01))), false);

    failed_test_count += tester.failed_test_count();
  } 
}; 

template <class T_SCALAR>
class test_minimal_area_triangulator: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static void check(sl::tester& tester,
		    const std::vector<point_t>& vertices,
		    const interval_t& expected_area) {
    sl::minimal_area_triangulator<3,value_t> tr;

    std::vector<sl::triangle_connectivity> triangles;

    tr.triangulation_in(triangles, vertices);
    float A = 0.0;
    for (std::size_t i=0; i<triangles.size(); ++i) {
      A += sl::triangle_area(vertices[triangles[i][0]],
			     vertices[triangles[i][1]],
			     vertices[triangles[i][2]]);
    }
    tester.test("triangle count", triangles.size(), vertices.size()-2);
    tester.test("area", A, expected_area);

  }

  static void do_it() {  
    sl::tester tester("minimal_area_triangulator < 3, " + sl::numeric_traits<value_t>::what() + " >");


    std::vector<point_t> vertices;

    {
      vertices.clear();
      vertices.push_back(point_t(0.0f, 0.0f, 0.0f));
      vertices.push_back(point_t(1.0f, 0.0f, 0.0f));
      vertices.push_back(point_t(0.0f, 2.0f, 0.0f));

      check(tester, vertices, interval_t("1.000"));
    }

    {
      vertices.clear();
      vertices.push_back(point_t(0.0f, 0.0f, 0.0f));
      vertices.push_back(point_t(1.0f, 0.0f, 0.0f));
      vertices.push_back(point_t(1.0f, 2.0f, 0.0f));
      vertices.push_back(point_t(0.0f, 2.0f, 0.0f));

      check(tester, vertices, interval_t("2.000"));
    }

    {
      vertices.clear();
      vertices.push_back(point_t(0.0f, 0.0f, 0.0f));
      vertices.push_back(point_t(3.0f, 0.0f, 0.0f));
      vertices.push_back(point_t(3.0f, 2.0f, 0.0f));
      vertices.push_back(point_t(2.0f, 2.0f, 0.0f));
      vertices.push_back(point_t(2.0f, 1.0f, 0.0f));
      vertices.push_back(point_t(1.0f, 1.0f, 0.0f));
      vertices.push_back(point_t(1.0f, 2.0f, 0.0f));
      vertices.push_back(point_t(0.0f, 2.0f, 0.0f));

      check(tester, vertices, interval_t("5.000"));
    }

    failed_test_count += tester.failed_test_count();
  } 
}; 

template <class T_SCALAR>
class test_shaft: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static void do_it() {  
    sl::tester tester("shaft < " + sl::numeric_traits<value_t>::what() + " >");

    {
      sl::axis_aligned_box<3,value_t> box0(point_t(-1.0f,-2.0f,-3.0f), point_t(1.0f,2.0f,3.0f));
      sl::axis_aligned_box<3,value_t> box1(point_t( 1.0f, 2.0f, 3.0f), point_t(3.0f,6.0f,9.0f));
      sl::axis_aligned_box<3,value_t> boxin(point_t( 0.9, 1.9, 2.9), point_t(1.1,2.1,3.1));
      sl::axis_aligned_box<3,value_t> boxout(point_t( 10.9, 11.9, 12.9), point_t(11.1,12.1,13.1));
      sl::oriented_box<3, value_t>    oboxin(boxin);
      sl::oriented_box<3, value_t>    oboxout(boxout);

      sl::shaft<value_t> shaft(box0, box1);
      
      tester.test("shaft contains", shaft.contains(point_t(-0.5f, -1.5f, -2.5f)),   true);
      tester.test("shaft contains", shaft.contains(point_t( 2.5f,  5.5f,  8.5f)),   true);
      tester.test("shaft contains", shaft.contains(point_t( 12.5f,  15.5f,  18.5f)),false);
      tester.test("shaft aabox test", shaft.contains(boxin),                     true);
      tester.test("shaft aabox test", shaft.contains(boxout),                    false);
      tester.test("shaft obox test", shaft.contains(oboxin),                     true);
      tester.test("shaft obox test", shaft.contains(oboxout),                    false);
      tester.test("shaft culls_out", shaft.culls_out(point_t(-0.5f, -1.5f, -2.5f)),   false);
      tester.test("shaft culls_out", shaft.culls_out(point_t( 2.5f,  5.5f,  8.5f)),   false);
      tester.test("shaft culls_out", shaft.culls_out(point_t( 12.5f,  15.5f,  18.5f)),true);
      tester.test("shaft aabox test", shaft.culls_out(boxin),                     false);
      tester.test("shaft aabox test", shaft.culls_out(boxout),                    true);
      tester.test("shaft obox test", shaft.culls_out(oboxin),                     false);
      tester.test("shaft obox test", shaft.culls_out(oboxout),                    true);
    }

    {
      sl::axis_aligned_box<3,value_t> box0(point_t(0.0f,0.0f,0.0f), point_t(1.0f,1.0f,1.0f));
      sl::axis_aligned_box<3,value_t> box1(point_t(-1.0f,-1.0f,0.5f), point_t(-0.5f,1.5f,1.5f));

      sl::shaft<value_t> shaft(box0, box1);
      sl::axis_aligned_box<3,value_t> box2(point_t(0.5f,0.5f,0.5f), point_t(3.5f,3.5f,3.5f));
 
      tester.test("shaft aabox test", shaft.contains(box2),                     false);
      tester.test("shaft aabox cull", shaft.culls_out(box2),                    false);

    }

    {
      const size_t N_TEST = 20;

      for (size_t i=0; i<N_TEST; ++i) {
	vector_t      t_0 = vector_t(10.1f*i,20.2f*i,30.3f*i);
	quaternion_t  r_0 = quaternion_t(vector_t(0.0f,0.0f,1.0f),value_t(0.17f*i));
	vector_t      t_1 = vector_t(10.3f*i,-20.3f*i,-30.3f*i);
	quaternion_t  r_1 = quaternion_t(vector_t(0.0f,1.0f,0.0f),value_t(0.17f*i));
	vector_t      h_l = vector_t(0.5f*i, 0.1f*i, 0.2f*i);

	rigid_body_map_t map1 = 
	  linear_map_factory_t::translation(t_0) *
	  linear_map_factory_t::rotation(r_0);
	
	rigid_body_map_t map2 = 
	  linear_map_factory_t::translation(t_1) *
	  linear_map_factory_t::rotation(r_1);

	sl::oriented_box<3,value_t> obox1(map1, h_l);
	sl::oriented_box<3,value_t> obox2(map2, h_l);

	sl::shaft<value_t> shaft(obox1, obox2);
 
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(1), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(1), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(2), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(2), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(3), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(3), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(5), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(5), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(6), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(6), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(7), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.corner(7), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(1), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(1), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(2), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(2), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(3), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(3), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(5), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(5), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(6), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(6), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(7), value_t(0.01f))));
	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox1.corner(7), value_t(0.99f))));

	tester.test("shaft contain", shaft.contains(sl::lerp(obox1.center(), obox2.center(), value_t(0.5f))));
      }
    }

    failed_test_count += tester.failed_test_count();
  } 
}; 


template <class T_SCALAR>
class test_oriented_box: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static void do_it3() {  
    sl::tester tester("oriented_box < 3, " + sl::numeric_traits<value_t>::what() + " >");
    
    sl::oriented_box<3,value_t> box;
    sl::oriented_box_builder<3, value_t> builder;
    
    point_t p0 = point_t(-1.0f,-2.0f, 0.0f);
    point_t p1 = point_t( 1.0f, 2.0f, 0.0f);
    
    builder.begin_model();
    builder.put_point(p0);
    builder.put_point(p1);
    builder.end_model();
    box = builder.last_bounding_volume();
    
    tester.test("box contains - 0", box.contains(p0.lerp(p1, value_t( 0.01f))), true);
    tester.test("box contains - 1", box.contains(p0.lerp(p1, value_t( 0.99f))), true);

#if 1
    //tester.test("box contains", box.contains(p0.lerp(p1, value_t( 1.01))), false);
    //tester.test("box contains", box.contains(p0.lerp(p1, value_t(-0.01f))), false);
    tester.test("box size - 0", box.half_side_lengths()[(box.half_side_lengths().iamax()+0)%3] - value_t(0.5f) * (p1-p0).two_norm(), interval_t("0.00")); 
      if (tester.has_last_test_failed()) {
	std::cerr << "PO = " << as_dual(p0.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "P1 = " << as_dual(p1.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "BOX.hsl = " << as_dual(box.half_side_lengths()) << std::endl;
	std::cerr << "BOX.map = " << std::endl << box.from_box_space_map() << std::endl;
      }
    tester.test("box size - 1", box.half_side_lengths()[(box.half_side_lengths().iamax()+1)%3] - value_t(0)    , interval_t("0.00")); 
      if (tester.has_last_test_failed()) {
	std::cerr << "PO = " << as_dual(p0.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "P1 = " << as_dual(p1.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "BOX.hsl = " << as_dual(box.half_side_lengths()) << std::endl;
	std::cerr << "BOX.map = " << std::endl << box.from_box_space_map() << std::endl;
      }
    tester.test("box size - 2", box.half_side_lengths()[(box.half_side_lengths().iamax()+2)%3] - value_t(0)    , interval_t("0.00")); 
      if (tester.has_last_test_failed()) {
	std::cerr << "PO = " << as_dual(p0.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "P1 = " << as_dual(p1.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "BOX.hsl = " << as_dual(box.half_side_lengths()) << std::endl;
	std::cerr << "BOX.map = " << std::endl << box.from_box_space_map() << std::endl;
      }
#endif

    {
      sl::oriented_box<3,value_t> obox;
      sl::axis_aligned_box<3,value_t> abox;
      sl::oriented_box_builder<3, value_t> obuilder;
      sl::axis_aligned_box_builder<3, value_t> abuilder;
      
      const size_t N_TEST = 40;
      std::vector< sl::fixed_size_point<3,value_t> > pts;
      for (size_t i=0; i<N_TEST; ++i) {
	sl::fixed_size_point<3,value_t> p = super_t::random_unit_point();
	p[2] = p[0]+p[1]+0.0001f*super_t::random_unit_value(); // Give a preferential directionm to test oriented box
	pts.push_back(p);
      }
      
      abuilder.begin_model();
      obuilder.begin_model();
      for (size_t i=0; i<N_TEST; ++i) {
	abuilder.put_point(pts[i]);
	obuilder.put_point(pts[i]);
      }
      obuilder.end_model();
      abuilder.end_model();
      
      obox = obuilder.last_bounding_volume();
      abox = abuilder.last_bounding_volume();
      
      tester.test("obox 3 volume <= abox volume", obox.volume() <= abox.volume());
      if (tester.has_last_test_failed()) {
	std::cerr << "O-BOX VOUME = " << obox.volume() << std::endl;
	std::cerr << "A-BOX VOUME = " << abox.volume() << std::endl;
      }

      for (size_t i=0; i<N_TEST; ++i) {
	point_t p = obox.center().lerp(pts[i], value_t(0.99f));

	tester.test("obox contains", obox.contains(p), true);
	if (tester.has_last_test_failed()) {
	  std::cerr << "LOCAL PT = " << as_dual(obox.to_box_space(p).as_vector()) << std::endl;
	  std::cerr << "HSL      = " << as_dual(obox.half_side_lengths()) << std::endl;
	}
	  
	tester.test("abox contains", abox.contains(p), true);
      }
    }

    {
      const size_t N_TEST = 20;

      for (size_t i=0; i<N_TEST; ++i) {
	vector_t      t_0 = vector_t(10.1*i,20.2*i,30.3*i);
	quaternion_t  r_0 = quaternion_t(vector_t(0,0,1),value_t(0.17*i));
	vector_t      h_l = vector_t(0.5f*i, 0.1*i, 0.2*i);

	rigid_body_map_t map = 
	  linear_map_factory_t::translation(t_0) *
	  linear_map_factory_t::rotation(r_0);

	sl::oriented_box<3,value_t> obox(map, h_l);
	sl::axis_aligned_box<3,value_t> abox = obox.bounding_axis_aligned_box();
 
	for (size_t j=0; j<8; ++j) {
	  tester.test("box contain", abox.contains(sl::lerp(obox.center(), obox.corner(j), value_t(0.01f))));
	  tester.test("box contain", abox.contains(sl::lerp(obox.center(), obox.corner(j), value_t(0.99f))));
	}

	tester.test("box dist", obox.signed_distance_range(plane_t::perpendicular_to(0))[0] - abox[0][0], interval_t("0.000"));
	tester.test("box dist", obox.signed_distance_range(plane_t::perpendicular_to(0))[1] - abox[1][0], interval_t("0.000"));
	tester.test("box dist", obox.signed_distance_range(plane_t::perpendicular_to(1))[0] - abox[0][1], interval_t("0.000"));
	tester.test("box dist", obox.signed_distance_range(plane_t::perpendicular_to(1))[1] - abox[1][1], interval_t("0.000"));
	tester.test("box dist", obox.signed_distance_range(plane_t::perpendicular_to(2))[0] - abox[0][2], interval_t("0.000"));
	tester.test("box dist", obox.signed_distance_range(plane_t::perpendicular_to(2))[1] - abox[1][2], interval_t("0.000"));
      }

      failed_test_count += tester.failed_test_count();
    }    
  } 

  static void do_it2() {  
    {
      sl::tester tester("oriented_box < 2, " + sl::numeric_traits<value_t>::what() + " >");
      
      sl::oriented_box<2,value_t> box;
      sl::oriented_box_builder<2, value_t> builder;
      
      sl::fixed_size_point<2,value_t> p0(-2.0f,-2.0f);
      sl::fixed_size_point<2,value_t> p1( 1.0f, 1.0f);
    
      builder.begin_model();
      builder.put_point(p0);
      builder.put_point(p1);
      builder.end_model();

      box = builder.last_bounding_volume();
      
#if 1
      tester.test("obox size - 0", box.half_side_lengths()[(box.half_side_lengths().iamax()+0)%2] - value_t(0.5f) * (p1-p0).two_norm(), interval_t("0.00")); 
      if (tester.has_last_test_failed()) {
	std::cerr << "PO = " << as_dual(p0.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "P1 = " << as_dual(p1.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "BOX.hsl = " << as_dual(box.half_side_lengths()) << std::endl;
	std::cerr << "BOX.map = " << std::endl << box.from_box_space_map() << std::endl;
      }
      tester.test("obox size - 1", box.half_side_lengths()[(box.half_side_lengths().iamax()+1)%2] - value_t(0.0f)    , interval_t("0.00")); 
      if (tester.has_last_test_failed()) {
	std::cerr << "PO = " << as_dual(p0.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "P1 = " << as_dual(p1.as_vector()) << "  - local = " << as_dual((box.from_box_space_map().inverse() * p0).as_vector()) << std::endl;
	std::cerr << "BOX.hsl = " << as_dual(box.half_side_lengths()) << std::endl;
	std::cerr << "BOX.map = " << std::endl << box.from_box_space_map() << std::endl;
      }
#endif
      
      {
	sl::oriented_box<2,value_t> obox;
	sl::axis_aligned_box<2,value_t> abox;
	sl::oriented_box_builder<2, value_t> obuilder;
	sl::axis_aligned_box_builder<2, value_t> abuilder;

	const size_t N_TEST = 20;
	std::vector< sl::fixed_size_point<2,value_t> > pts;
	for (size_t i=0; i<N_TEST; ++i) {
	  pts.push_back(sl::fixed_size_point<2,value_t>(value_t(-1.0f * (i * 13 % 17)),
							value_t(1.0f * (i * 15 % 23))));
	}

	abuilder.begin_model();
	obuilder.begin_model();
	for (size_t i=0; i<N_TEST; ++i) {
	  abuilder.put_point(pts[i]);
	  obuilder.put_point(pts[i]);
	}
	obuilder.end_model();
	abuilder.end_model();
	
	obox = obuilder.last_bounding_volume();
	abox = abuilder.last_bounding_volume();

	SL_TRACE_OUT(1) << "###############O-BOX VOUME = " << obox.volume() << std::endl;
	SL_TRACE_OUT(1) << "###############A-BOX VOUME = " << abox.volume() << std::endl;

	tester.test("obox cvolume <= abox volume", obox.volume() <= abox.volume());
	for (size_t i=0; i<N_TEST; ++i) {
	  sl::fixed_size_point<2,value_t> p = obox.center().lerp(pts[i], value_t(0.99f));

	  tester.test("obox contains", obox.contains(p), true);
	  if (tester.has_last_test_failed()) {
	    std::cerr << "LOCAL PT = " << as_dual(obox.to_box_space(p).as_vector()) << std::endl;
	    std::cerr << "HSL      = " << as_dual(obox.half_side_lengths()) << std::endl;
	  }
	  tester.test("abox contains", abox.contains(p), true);
	}
      }

      failed_test_count += tester.failed_test_count();
    }
  }


  static void do_it() {
    do_it3();
    do_it2();
  }

};

template <class T_SCALAR>
class test_convex_hull: public test_geometry_base<T_SCALAR> {
public:
  typedef test_geometry_base<T_SCALAR>       super_t;

  typedef typename super_t::value_t          value_t;
  typedef typename super_t::interval_t       interval_t;

  typedef typename super_t::linear_map_factory_t  linear_map_factory_t;
  typedef typename super_t::affine_map_t          affine_map_t;
  typedef typename super_t::rigid_body_map_t      rigid_body_map_t;
  typedef typename super_t::projective_map_t      projective_map_t;

  typedef typename super_t::point_t              point_t;
  typedef typename super_t::vector_t             vector_t;
  typedef typename super_t::dual_vector_t        dual_vector_t;
  typedef typename super_t::quaternion_t         quaternion_t;
  typedef typename super_t::plane_t              plane_t;
public:

  static bool has_vertex(const std::vector< sl::simplex_connectivity<2> >& h, size_t idx) {
    bool result = false;
    for (size_t i=0; i<h.size() && !result; ++i) {
      result = (idx == h[i][0]);
    }
    return result;
  }

  static void do_it() {  
    const size_t hull_method_count = 3;

    const char* hull_method_id[hull_method_count] = {
      "DEFAULT ALGORITHM",
      "CHAIN ALGORITHM",
      "INCREMENTAL ALGORITHM"
    };

    for (size_t hull_method = 0; hull_method<hull_method_count; ++hull_method) {
      sl::tester tester("convex_hull_builder <2," + sl::numeric_traits<value_t>::what() + " >" + 
			" -- " + hull_method_id[hull_method]);
	
      typedef sl::fixed_size_point<2,value_t> point_t;
      
      { 
	std::vector<point_t> pts;
	pts.push_back(point_t(value_t(10.0f),value_t(20.0f)));
	pts.push_back(point_t(value_t(30.0f),value_t(40.0f)));
	pts.push_back(point_t(value_t(10.0f),value_t(10.0f)));
	pts.push_back(point_t(value_t( 5.0f),value_t( 5.1f)));
	pts.push_back(point_t(value_t( 1.0f),value_t( 1.0f)));
	
	sl::convex_hull_builder< 2,value_t,std::vector<point_t> > chull;
	switch (hull_method) {
	case 0: chull.build(pts); break;
	case 1: chull.build_chain(pts); break;
	case 2: chull.build_incremental(pts); break;
	default: SL_FAIL("Unknown hull method");
	}
	
	typedef typename sl::convex_hull_builder< 2,value_t, std::vector<point_t> >::hull_t hull_t;
        typedef typename sl::convex_hull_builder< 2,value_t, std::vector<point_t> >::simplex_t simplex_t;
        hull_t known_hull;
	known_hull.push_back(simplex_t(4,2));
	known_hull.push_back(simplex_t(2,1));
	known_hull.push_back(simplex_t(1,0));
	known_hull.push_back(simplex_t(0,4));
 
	tester.test("Hull size", chull.last_hull().size(), known_hull.size());
	if (chull.last_hull().size() == known_hull.size()) {
	  const size_t n =  known_hull.size();
	  size_t i0 = 0;
	  for (i0 = 0; i0<n; ++i0) {
	    if (chull.last_hull()[i0][0] == known_hull[0][0]) break;
	  }
	  for (size_t i=0; i<n; ++i) {
	    tester.test("Hull vertex org", chull.last_hull()[(i+i0)%n][0], known_hull[i][0]);
	    tester.test("Hull vertex ext", chull.last_hull()[(i+i0)%n][1], known_hull[i][1]);
	  }
	}
      }

      {
	const size_t N_TEST = 20;
	std::vector< point_t > pts;
	for (size_t i=0; i<N_TEST; ++i) {
	  pts.push_back(sl::fixed_size_point<2,value_t>(value_t(-1.0f * (i * 13 % 17)),
							value_t(1.0f * (i * 15 % 23))));
	}

	sl::convex_hull_builder< 2,value_t,std::vector<point_t> > chull;
	switch (hull_method) {
	case 0: chull.build(pts); break;
	case 1: chull.build_chain(pts); break;
	case 2: chull.build_incremental(pts); break;
	default: SL_FAIL("Unknown hull method");
	}
	
	std::vector< sl::simplex_connectivity<2> > h = chull.last_hull();
	sl::simplex_to_poly_wrapper<value_t, std::vector< point_t> > poly( h, pts );
	
	tester.test("Hull is convex", sl::is_ccw_convex_polygon(poly));
	
	for (size_t i=0; i<N_TEST; ++i) {
	  tester.test("Hull containment", 
		      has_vertex(h,i) || 
		      sl::point_in_ccw_convex_polygon(poly,pts[i]));
	}
	
      }

      failed_test_count += tester.failed_test_count();
    }
  }
};

class test_conv_to {
public:

  static void do_it() {
    sl::tester tester("conv_to<...>");

    tester.test("point3f -> point3d",
		(sl::point3f(0.1f,0.2f,-10.1f)-
		 sl::conv_to<sl::point3f>::from(sl::point3d(0.1,0.2,-10.1))).infinite_norm(),
		sl::intervalf("0.000000"));
    tester.test("point3d -> point3f",
		(sl::point3d(0.1f,0.2f,-10.1f)-
		 sl::conv_to<sl::point3d>::from(sl::point3f(0.1,0.2,-10.1))).infinite_norm(),
		sl::intervald("0.000000"));
   }
};


template <class T_SCALAR>
class test_suite {
public:

  static void do_it() {
    test_quadric<T_SCALAR>::do_it();
    test_point<T_SCALAR>::do_it();
    test_plane<T_SCALAR>::do_it();
    test_quaternion<T_SCALAR>::do_it();
    test_affine_map<T_SCALAR>::do_it();
    test_rigid_body_map<T_SCALAR>::do_it();
    test_projective_map<T_SCALAR>::do_it();
    test_cartesian_frame<T_SCALAR>::do_it();
    test_axis_aligned_box<T_SCALAR>::do_it();
    test_convex_hull<T_SCALAR>::do_it();
    test_oriented_box<T_SCALAR>::do_it();
    test_shaft<T_SCALAR>::do_it();
    test_minimal_area_triangulator<T_SCALAR>::do_it();
  }    

};

int main() {

  test_suite<float>::do_it();
  test_conv_to::do_it();
  
#if 0
  test_suite<double>::do_it();
  test_suite<sl::bounded_scalar<float> >::do_it();
  test_suite<sl::bounded_scalar<double> >::do_it();
#endif

  return (int)failed_test_count;

}


