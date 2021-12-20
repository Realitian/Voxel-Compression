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
#include <sl/linear_map_factory.hpp>
#include <sl/lsq_absolute_orientation.hpp>
#include <sl/ransac_absolute_orientation.hpp>
#include <vector>
#include <sl/tester.hpp>

static std::size_t failed_test_count = 0;

void test_lsq_absolute_orientation_no_scaling() {
  sl::tester tester("Least squares absolute orientation - no scaling");

  std::size_t N = 12;

  std::vector<sl::point3f> a;
  std::vector<sl::point3f> b;

  sl::rigid_body_map3f tR = 
    sl::linear_map_factory3f::translation(-10,-20,30)*
    sl::linear_map_factory3f::rotation(0.3,0.5,0.7);

  for (std::size_t i=0; i<N; ++i) {
    sl::point3f a_i = sl::point3f((i*13)%17,
				  -(i*17+15)%13,
				  -(i*23)%7);
    sl::point3f b_i = tR * a_i;
    a.push_back(a_i);
    b.push_back(b_i);
  }

  sl::rigid_body_map3f tR_star;
  sl::lsq_absolute_orientation_in(tR_star,
				  N,
				  (float*)0,
				  &(a[0]),
				  &(b[0]));
  float error = 0.0;
  for (std::size_t i=0; i<N; ++i) {
    sl::point3f a_i = a[i];
    sl::point3f b_i = b[i];
    sl::point3f b_i_star = tR_star * a_i;

    error += b_i.distance_squared_to(b_i_star);
  }

  tester.test("Distance squared", error, sl::intervalf("0.0000"));

  failed_test_count += tester.failed_test_count();
}

void test_lsq_absolute_orientation_with_scaling() {
  sl::tester tester("Least squares absolute orientation - with scaling");

  std::size_t N = 12;

  std::vector<sl::point3f> a;
  std::vector<sl::point3f> b;

  sl::affine_map3f tRs = 
    sl::linear_map_factory3f::translation(-10,-20,30) *
    sl::linear_map_factory3f::rotation(0.3,0.5,0.7)*
    sl::linear_map_factory3f::scaling(0.3);

  for (std::size_t i=0; i<N; ++i) {
    sl::point3f a_i = sl::point3f((i*13)%17,
				  -(i*17+15)%13,
				  -(i*23)%7);
    sl::point3f b_i = tRs * a_i;
    a.push_back(a_i);
    b.push_back(b_i);
  }

  sl::affine_map3f tRs_star;
  sl::lsq_absolute_orientation_in(tRs_star,
				  N,
				  (float*)0,
				  &(a[0]),
				  &(b[0]));
  float error = 0.0;
  for (std::size_t i=0; i<N; ++i) {
    sl::point3f a_i = a[i];
    sl::point3f b_i = b[i];
    sl::point3f b_i_star = tRs_star * a_i;

    error += b_i.distance_squared_to(b_i_star);
  }

  tester.test("Distance squared", error, sl::intervalf("0.0000"));

  failed_test_count += tester.failed_test_count();
}

void test_ransac_absolute_orientation() {
  sl::tester tester("RANSAC Absolute orientation - with scaling");

  std::size_t N = 1000;
    
  std::vector<sl::point3f> a_inlier;
  std::vector<sl::point3f> b_inlier;

  std::vector<sl::point3f> a_contaminated;
  std::vector<sl::point3f> b_contaminated;

  sl::affine_map3f tRs = 
    sl::linear_map_factory3f::translation(-10,-20, 30) *
    sl::linear_map_factory3f::rotation(0.3,0.5,0.7)*
    sl::linear_map_factory3f::scaling(0.3);

  sl::random::uniform<float> rng;

  float noise_amount = 1.0f;
  for (std::size_t i=0; i<N; ++i) {
    sl::point3f a_i;
    sl::point3f b_i;
    
    a_i = sl::point3f(rng.value(), rng.value(), rng.value()).scaled_by(100.0f);
    if (i%3==0) {
      // Good
      b_i = tRs * a_i + noise_amount / std::sqrt(3.0f) * sl::vector3f(-0.5f+rng.value(),
								      -0.5f+rng.value(),
								      -0.5f+rng.value());
      a_inlier.push_back(a_i);
      b_inlier.push_back(b_i);
    } else {
      // Bad
      b_i = sl::point3f(rng.value(), rng.value(), rng.value()).scaled_by(100.0f);
    }
    a_contaminated.push_back(a_i);
    b_contaminated.push_back(b_i);
  }

  // Reconstruct using contaminated set
  sl::affine_map3f tRs_star;
  sl::ransac_absolute_orientation_in(tRs_star,
				     a_contaminated.size(),
				     &(a_contaminated[0]),
				     &(b_contaminated[0]),
				     2*noise_amount);

  // Check inlier fitting
  for (std::size_t i=0; i<a_inlier.size(); ++i) {
    sl::point3f a_i = a_inlier[i];
    sl::point3f b_i = b_inlier[i];
    sl::point3f b_i_star = tRs_star * a_i;

    float err_i = b_i.distance_squared_to(b_i_star);

    // std::cerr << err_i << std::endl;
    tester.test("Fitting", err_i, sl::intervalf(-noise_amount,noise_amount));
  }
  
  failed_test_count += tester.failed_test_count();
}

int main() {
  test_lsq_absolute_orientation_no_scaling();
  test_lsq_absolute_orientation_with_scaling();
  test_ransac_absolute_orientation();
  return failed_test_count;
}
