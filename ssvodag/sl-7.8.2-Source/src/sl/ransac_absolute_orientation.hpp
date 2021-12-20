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
#ifndef SL_RANSAC_ABSOLUTE_ORIENTATION_HPP
#define SL_RANSAC_ABSOLUTE_ORIENTATION_HPP

#include <sl/lsq_absolute_orientation.hpp>
#include <sl/random.hpp>
#include <sl/utility.hpp>

namespace sl {

  /**
   * Robustly find affine map such that b = s R(a) + t
   *
   * is_forcing_unit_scaling is true iff ignoring scaling in the optimization
   * confidence_threshold is the confidence with which we don't find a solution
   * with all inliers
   * min_inlier_fraction is the estimated minimum number of inliers (used only
   * as a cutoff on the number of iterations)
   *
   * Based on:
   *   D. Capel, An effective bail-out test for RANSAC consensus scoring,
   *   Proc. BMVC, 629-638, 2005.
   *   O. Chum, J. Matas, J. Kittler, Locally optimized RANSAC,
   *   Pattern Recognition, 2003: 236-242.
   *   P. H. S. Torr and A. Zisserman. MLESAC: A new robust estimator with
   *   application to estimating image geometry. Computer Vision and Image
   *   Understanding, 78: 138–156, 2000.
   *   M. A. Fischler and R. C. Bolles. Random sample consensus: A paradigm
   *   for model fitting with applications to image analysis and
   *   automated cartography. Comm. Assoc. Comp. Mach., 24(6):381–395, 1981.
   */   
  template <class T>
  void ransac_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				      sl::quaternion<T>&                                 q_star,
				      T&                                                 s_star,
				      const std::size_t N,
				      const sl::fixed_size_point<3,T>* a,
				      const sl::fixed_size_point<3,T>* b,
				      bool  is_forcing_unit_scaling,
				      const T& inlier_threshold,
				      double confidence_threshold = 0.01,
				      double min_inlier_fraction = 0.1) {
    assert(N>0);

    const std::size_t dimension = 3;
    typedef T                                                value_t;
    typedef sl::fixed_size_point<dimension, value_t>         point_t;
    typedef typename point_t::vector_t                       vector_t;
    typedef sl::quaternion<value_t>                          quaternion_t;
    typedef sl::rigid_body_map<3,value_t>                    rigid_body_map_t;
    typedef sl::affine_map<3,value_t>                        affine_map_t;
    typedef sl::linear_map_factory<3,value_t>                linear_map_factory_t;
    
    random::std_irng_t rng; // FIXME what about seed?

    t_star.to_zero();
    q_star.to_identity();
    s_star = value_t(1.0);
    double e_star = 1e30;

    const std::size_t M = 3; // Size of basis set

    if (N<=M) {
      // Simple problem, use least squares solution
      if (is_forcing_unit_scaling) {
	lsq_absolute_orientation_in(t_star, q_star, 
				    N, a, b);
      } else {
	lsq_absolute_orientation_in(t_star, q_star, s_star,
				    N, a, b);
      }	
    } else {
      // More than 3 correspondences, use RANSAC approach

      // Estimate max iterations from number of inliers
      const double      log_eta = std::log(confidence_threshold);
      const double      x       = 1.0-std::pow(min_inlier_fraction,int(M));
      const std::size_t n_max0  = sl::median(std::size_t(10),
					     std::size_t(10000),
					     std::size_t(log_eta/std::log(sl::median(x,0.001,0.9999))+0.5f));
      
      std::size_t N_star = 0;
      std::size_t n = 0;
      std::size_t n_max = n_max0;
      while (n<n_max) {
	// Sample a basis set of size m from correspondances a,b
	uint32_t idx_i[M];
	rng.pick_k_out_of_n_in(M, idx_i,uint32_t(N));

	// std::cerr << "ITER: " << n << " BASIS = ";
	std::vector<point_t> a_i;
	std::vector<point_t> b_i;
	for (std::size_t k=0; k<M; ++k) {
	  a_i.push_back(a[idx_i[k]]);
	  b_i.push_back(b[idx_i[k]]);

	  // std::cerr << idx_i[k] << " ";
	}
	// std::cerr << std::endl;
	
	// Generate current alignment hypothesis
	vector_t     t_i; 
	quaternion_t q_i; 
	value_t      s_i;
	double       e_i = e_star;

	// Estimate error and refine if better than current best (LO-RANSAC)
	for (std::size_t lo_i=0; lo_i<2; ++lo_i) {
	  const value_t thr2 = inlier_threshold*inlier_threshold;
	  t_i.to_zero();
	  q_i.to_zero();
	  s_i = value_t(1.0);
	  if (is_forcing_unit_scaling) {
	    lsq_absolute_orientation_in(t_i, q_i, a_i.size(), &(a_i[0]), &(b_i[0]));
	  } else {
	    lsq_absolute_orientation_in(t_i, q_i, s_i, a_i.size(), &(a_i[0]), &(b_i[0]));
	  }
	  affine_map_t tRs_i =
	    rigid_body_map_t(t_i) *
	    rigid_body_map_t(q_i) *
	    linear_map_factory_t::scaling(s_i);

	  // Rebuild inlier set and update error;
	  a_i.clear();
	  b_i.clear();
	  e_i = 0.0;
	  
	  for (std::size_t k=0; k<N; ++k) {
	    // MLESAC robust estimator
	    point_t b_prime = tRs_i * a[k];
	    value_t d2 = b[k].distance_squared_to(b_prime);
	    if (d2<thr2) {
	      a_i.push_back(a[k]);
	      b_i.push_back(b[k]);
	      e_i += d2;
	    } else {
	      e_i += thr2;
	    }
	    
	    // Bail out if worse than current best
	    if (e_i > e_star) break;
	  }

	  // std::cerr << "LO ITER: " << lo_i << " => e = " << e_i << " vs. " << e_star << std::endl;
	  // Bail out if worse than current best
	  if (e_i > e_star) break;

	  // Bail out if not enough inliers
	  if (a_i.size()<M) break;
	} // local optimization loop

	// Update best alignment and error
	if (e_i<e_star) {
	  // std::cerr << "UPDATE BEST!: INNER = " << a_i.size() << std::endl;
	  t_star = t_i;
	  q_star = q_i;
	  s_star = s_i;
	  e_star = e_i;
	  N_star = a_i.size();
	}

	// Recompute estimate of number of inliers over
	// number of elements
	double eps = double(N_star) / double(N);

	// Given the true fraction of inlying correspondences eps,
	// the probability of selecting a basis set of size M that
	// consists entirely of inliers is eps^M. Hence the probability
	// of sampling n basis sets all of which are polluted by at
	// least one outlier is eta = (1−eps^M)^n
	// Therefore the minimum number of samples n_max that must be
	// taken in order that this probability falls below a given
	// confidence threshold eta_star is given by:
	double x = 1.0-std::pow(eps,int(M));
	if (x<=0.0) {
	  // All inliers
	  n_max = n; 
	} else if (x>=1.0) {
	  // No inliers at all -- should not occurr
	  n_max = 100000;
	} else {
	  // Estimate iterations from number of inliers
	  n_max = std::size_t(log_eta/std::log(x)+0.5f);
	}
	n_max = std::min(n_max0, n_max);
	
	++n;

	// std::cerr << "RANSAC: Iter " << n << "/" << n_max << ": ERR=" << e_star << " IN: " << N_star << "/" << N << " = " << sl::human_readable_percent(100.0*eps) << std::endl;
      }
    }
  }

  /**
   * Robustly find affine map such that b = s R(a) + t
   *
   * is_forcing_unit_scaling is true iff ignoring scaling in the optimization
   * confidence_threshold is the confidence with which we don't find a solution
   * with all inliers
   * min_inlier_fraction is the estimated minimum number of inliers (used only
   * as a cutoff on the number of iterations)
   *
   * Based on:
   *   D. Capel, An effective bail-out test for RANSAC consensus scoring,
   *   Proc. BMVC, 629-638, 2005.
   *   O. Chum, J. Matas, J. Kittler, Locally optimized RANSAC,
   *   Pattern Recognition, 2003: 236-242.
   *   P. H. S. Torr and A. Zisserman. MLESAC: A new robust estimator with
   *   application to estimating image geometry. Computer Vision and Image
   *   Understanding, 78: 138–156, 2000.
   *   M. A. Fischler and R. C. Bolles. Random sample consensus: A paradigm
   *   for model fitting with applications to image analysis and
   *   automated cartography. Comm. Assoc. Comp. Mach., 24(6):381–395, 1981.
   */   
  template <class T>
  void ransac_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				      sl::quaternion<T>&                                 q_star,
				      T&                                                 s_star,
				      const std::size_t N,
				      const sl::fixed_size_point<3,T>* a,
				      const sl::fixed_size_point<3,T>* b,
				      bool  is_forcing_unit_scaling,
				      const T* inlier_threshold,
				      double confidence_threshold = 0.01,
				      double min_inlier_fraction = 0.1) {
    assert(N>0);

    const std::size_t dimension = 3;
    typedef T                                                value_t;
    typedef sl::fixed_size_point<dimension, value_t>         point_t;
    typedef typename point_t::vector_t                       vector_t;
    typedef sl::quaternion<value_t>                          quaternion_t;
    typedef sl::rigid_body_map<3,value_t>                    rigid_body_map_t;
    typedef sl::affine_map<3,value_t>                        affine_map_t;
    typedef sl::linear_map_factory<3,value_t>                linear_map_factory_t;
    
    random::std_irng_t rng; // FIXME what about seed?

    t_star.to_zero();
    q_star.to_identity();
    s_star = value_t(1.0);
    double e_star = 1e30;

    const std::size_t M = 3; // Size of basis set

    if (N<=M) {
      // Simple problem, use least squares solution
      if (is_forcing_unit_scaling) {
	lsq_absolute_orientation_in(t_star, q_star, 
				    N, a, b);
      } else {
	lsq_absolute_orientation_in(t_star, q_star, s_star,
				    N, a, b);
      }	
    } else {
      // More than 3 correspondences, use RANSAC approach

      // Estimate max iterations from number of inliers
      const double      log_eta = std::log(confidence_threshold);
      const double      x       = 1.0-std::pow(min_inlier_fraction,int(M));
      const std::size_t n_max0  = sl::median(std::size_t(10),
					     std::size_t(10000),
					     std::size_t(log_eta/std::log(sl::median(x,0.001,0.9999))+0.5f));
      
      std::size_t N_star = 0;
      std::size_t n = 0;
      std::size_t n_max = n_max0;
      while (n<n_max) {
	// Sample a basis set of size m from correspondances a,b
	uint32_t idx_i[M];
	rng.pick_k_out_of_n_in(M, idx_i,uint32_t(N));

	// std::cerr << "ITER: " << n << " BASIS = ";
	std::vector<point_t> a_i;
	std::vector<point_t> b_i;
	for (std::size_t k=0; k<M; ++k) {
	  a_i.push_back(a[idx_i[k]]);
	  b_i.push_back(b[idx_i[k]]);

	  // std::cerr << idx_i[k] << " ";
	}
	// std::cerr << std::endl;
	
	// Generate current alignment hypothesis
	vector_t     t_i; 
	quaternion_t q_i; 
	value_t      s_i;
	double       e_i = e_star;

	// Estimate error and refine if better than current best (LO-RANSAC)
	for (std::size_t lo_i=0; lo_i<2; ++lo_i) {
	  t_i.to_zero();
	  q_i.to_zero();
	  s_i = value_t(1.0);
	  if (is_forcing_unit_scaling) {
	    lsq_absolute_orientation_in(t_i, q_i, a_i.size(), &(a_i[0]), &(b_i[0]));
	  } else {
	    lsq_absolute_orientation_in(t_i, q_i, s_i, a_i.size(), &(a_i[0]), &(b_i[0]));
	  }
	  affine_map_t tRs_i =
	    rigid_body_map_t(t_i) *
	    rigid_body_map_t(q_i) *
	    linear_map_factory_t::scaling(s_i);

	  // Rebuild inlier set and update error;
	  a_i.clear();
	  b_i.clear();
	  e_i = 0.0;
	  
	  for (std::size_t k=0; k<N; ++k) {
	    // MLESAC robust estimator
	    point_t b_prime = tRs_i * a[k];
	    value_t d2 = b[k].distance_squared_to(b_prime);
	    const value_t thr = inlier_threshold[k];
	    const value_t thr2 = thr*thr;
	    if (d2<thr2) {
	      a_i.push_back(a[k]);
	      b_i.push_back(b[k]);
	      e_i += d2;
	    } else {
	      e_i += thr2;
	    }
	    
	    // Bail out if worse than current best
	    if (e_i > e_star) break;
	  }

	  // std::cerr << "LO ITER: " << lo_i << " => e = " << e_i << " vs. " << e_star << std::endl;
	  // Bail out if worse than current best
	  if (e_i > e_star) break;
	} // local optimization loop

	// Update best alignment and error
	if (e_i<e_star) {
	  // std::cerr << "UPDATE BEST!: INNER = " << a_i.size() << std::endl;
	  t_star = t_i;
	  q_star = q_i;
	  s_star = s_i;
	  e_star = e_i;
	  N_star = a_i.size();
	}

	// Recompute estimate of number of inliers over
	// number of elements
	double eps = double(N_star) / double(N);

	// Given the true fraction of inlying correspondences eps,
	// the probability of selecting a basis set of size M that
	// consists entirely of inliers is eps^M. Hence the probability
	// of sampling n basis sets all of which are polluted by at
	// least one outlier is eta = (1−eps^M)^n
	// Therefore the minimum number of samples n_max that must be
	// taken in order that this probability falls below a given
	// confidence threshold eta_star is given by:
	double x = 1.0-std::pow(eps,int(M));
	if (x<=0.0) {
	  // All inliers
	  n_max = n; 
	} else if (x>=1.0) {
	  // No inliers at all -- should not occurr
	  n_max = 100000;
	} else {
	  // Estimate iterations from number of inliers
	  n_max = std::size_t(log_eta/std::log(x)+0.5f);
	}
	n_max = std::min(n_max0, n_max);
	
	++n;

	// std::cerr << "RANSAC: Iter " << n << "/" << n_max << ": ERR=" << e_star << " IN: " << N_star << "/" << N << " = " << sl::human_readable_percent(100.0*eps) << std::endl;
      }
    }
  }

  // Robustly find affine map such that b = s R(a) + t
  template <class T>
  void ransac_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				      sl::quaternion<T>&                                 q_star,
				      T&                                                 s_star,
				      const std::size_t N,
				      const sl::fixed_size_point<3,T>* a,
				      const sl::fixed_size_point<3,T>* b,
				      const T& inlier_threshold,
				      double confidence_threshold = 0.01,
				      double min_inlier_fraction = 0.1) {
    ransac_absolute_orientation_in(t_star, q_star, s_star,
				   N, a, b,
				   false, inlier_threshold,
				   confidence_threshold, min_inlier_fraction);
  }

  // Robustly find affine map such that b = R(a) + t
  template <class T>
  void ransac_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				      sl::quaternion<T>&                                 q_star,
				      const std::size_t N,
				      const sl::fixed_size_point<3,T>* a,
				      const sl::fixed_size_point<3,T>* b,
				      const T& inlier_threshold,
				      double confidence_threshold = 0.01,
				      double min_inlier_fraction = 0.1) {
    T s_star,
    ransac_absolute_orientation_in(t_star, q_star, s_star,
				   N, a, b,
				   true, inlier_threshold,
				   confidence_threshold, min_inlier_fraction);
  }

  // Robustly find affine map such that b = s R(a) + t
  template <class T>
  void ransac_absolute_orientation_in(sl::affine_map<3,T>& tRs,
				      const std::size_t N,
				      const sl::fixed_size_point<3,T>* a,
				      const sl::fixed_size_point<3,T>* b,
				      const T& inlier_threshold,
				      double confidence_threshold = 0.01,
				      double min_inlier_fraction = 0.1) {
    sl::fixed_size_vector<sl::column_orientation,3,T>  t_star;
    sl::quaternion<T>                                  q_star;
    T                                                  s_star;
    ransac_absolute_orientation_in(t_star, q_star, s_star,
				   N, a, b,
				   false, inlier_threshold,
				   confidence_threshold, min_inlier_fraction);
    tRs =
      sl::rigid_body_map<3,T>(t_star) *
      sl::rigid_body_map<3, T>(q_star) *
      sl::linear_map_factory<3,T>::scaling(s_star);
  }

  // Robustly find affine map such that b = R(a) + t
  template <class T>
  void ransac_absolute_orientation_in(sl::rigid_body_map<3,T>& tR,
				      const std::size_t N,
				      const sl::fixed_size_point<3,T>* a,
				      const sl::fixed_size_point<3,T>* b,
				      const T& inlier_threshold,
				      double confidence_threshold = 0.01,
				      double min_inlier_fraction = 0.1) {
    sl::fixed_size_vector<sl::column_orientation,3,T> t_star;
    sl::quaternion<T>                                 q_star;
    T                                                 s_star;
    ransac_absolute_orientation_in(t_star, q_star, s_star,
				   N, a, b,
				   true, inlier_threshold,
				   confidence_threshold, min_inlier_fraction);
    tR =
      sl::rigid_body_map<3,T>(t_star) *
      sl::rigid_body_map<3, T>(q_star);
   }

} // namespace sl

#endif
