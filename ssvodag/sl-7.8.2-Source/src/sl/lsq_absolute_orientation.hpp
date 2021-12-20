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
#ifndef SL_LSQ_ABSOLUTE_ORIENTATION_HPP
#define SL_LSQ_ABSOLUTE_ORIENTATION_HPP

#include <iostream>
#include <sl/rigid_body_map.hpp>
#include <sl/linear_map_factory.hpp>
#include <sl/quaternion.hpp>
#include <sl/fixed_size_square_matrix.hpp>
#include <sl/fixed_size_vector.hpp>
#include <sl/fixed_size_point.hpp>
#include <cassert>

namespace sl {

  /** 
   *  Find weighted least square optimal affine map such that b=s R(a)+t
   *  The map minimizes the error sum_i ||m_i * (b_i-sR(a_i)-t) ||
   *  Based on:
   *  HORN B.: Closed-form solution of absolute orientation
   *  using unit quaternions. Journal of Optical Society
   *  of America 4 (April 1987), 629-642.
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				   sl::quaternion<T>&                                 q_star,
				   T&                                                 s_star,
				   const std::size_t N,
				   const T* mass,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b,
				   bool  with_scaling = true) {
    assert(N>0);

    const std::size_t dimension = 3;
    typedef T                                                value_t;
    typedef sl::fixed_size_point<dimension, value_t>         point_t;
    typedef typename point_t::vector_t                       vector_t;
    typedef sl::fixed_size_square_matrix<4, value_t>         matrix4_t;
    typedef typename matrix4_t::vector_t                     vector4_t;
    typedef sl::quaternion<value_t>                          quaternion_t;
    typedef sl::rigid_body_map<3,value_t>                    rigid_body_map_t;

    // Default values
    t_star.to_zero();
    q_star.to_identity();
    s_star = value_t(1.0);
    
    // Compute centers of mass
    point_t a_cm;
    point_t b_cm;
    value_t sum_m_i = value_t(0.0);
    for (std::size_t i=0; i<N; ++i) {
      value_t m_i = (mass) ? (mass[i]) : value_t(1.0);
      a_cm+= m_i * a[i].as_vector();
      b_cm+= m_i * b[i].as_vector();
      sum_m_i+= m_i;
    }
    if (sum_m_i == value_t(0.0)) {
      // Zero mass?? return identity xform
      return;
    } else {
      // Find centers of mass
      const value_t inv_sum_m_i = sl::reciprocal(sum_m_i);
      for (std::size_t d=0; d<dimension; ++d) {
	a_cm[d] *= inv_sum_m_i;
	b_cm[d] *= inv_sum_m_i;
      }
      
      // Find optimal rotations via eigenvalue problem
      // Find optimal scaling via weighted averaging of distances to centers of mass
      value_t sum_d2a_i=value_t(0.0);
      value_t sum_d2b_i=value_t(0.0);
      matrix4_t sumAtB;
      for (std::size_t i=0; i<N; ++i) {
	value_t m_i = (mass) ? (mass[i]) : value_t(1.0);
	value_t w_i = m_i * inv_sum_m_i; // Not necessary, but keeps weight in nice range
	  
	vector_t alpha_i = a[i]-a_cm;
	vector_t beta_i  = b[i]-b_cm;
	matrix4_t Ai; Ai =
			0.0,        -alpha_i[0], -alpha_i[1], -alpha_i[2],
			alpha_i[0],         0.0,  alpha_i[2], -alpha_i[1],
			alpha_i[1], -alpha_i[2],         0.0,  alpha_i[0],
			alpha_i[2],  alpha_i[1], -alpha_i[0],        0.0;
	matrix4_t Bi; Bi = 
			0.0,       -beta_i[0], -beta_i[1], -beta_i[2],
			beta_i[0],        0.0, -beta_i[2],  beta_i[1],
			beta_i[1],  beta_i[2],        0.0, -beta_i[0],
			beta_i[2], -beta_i[1],  beta_i[0],       0.0;
      
	sumAtB += w_i * Ai.transposed()*Bi;

	sum_d2a_i += w_i * alpha_i.dot(alpha_i);
	sum_d2b_i += w_i * beta_i.dot(beta_i);
      }

      // Extract scaling
      if (with_scaling) {
	s_star = sum_d2a_i ? value_t(std::sqrt(sum_d2b_i/sum_d2a_i)) : value_t(1.0);
      }
      
      // Solve eigenproblem
      vector4_t sumAtB_eigenvals;
      matrix4_t sumAtB_eigenvecs;
      bool      ok;
      sumAtB.symmetric_sorted_eigen_in(sumAtB_eigenvals, sumAtB_eigenvecs, &ok);
      if (!ok) {
	// FIXME unable to solve eigenvalue problem - Assume identity rotation ???
	t_star = b_cm-a_cm;
      } else {
	q_star = quaternion_t(sumAtB_eigenvecs(1,0), 
			      sumAtB_eigenvecs(2,0), 
			      sumAtB_eigenvecs(3,0), 
			      sumAtB_eigenvecs(0,0));
	rigid_body_map_t R_star = rigid_body_map_t(q_star);
	t_star = b_cm - (R_star * a_cm).scaled_by(s_star); 
      }
    }
  }

  /** 
   *  Find weighted least square optimal affine map such that b=s R(a)+t
   *  The map minimizes the error sum_i ||m_i * (b_i-R(a_i)-t) ||
   *  Based on:
   *  HORN B.: Closed-form solution of absolute orientation
   *  using unit quaternions. Journal of Optical Society
   *  of America 4 (April 1987), 629-642.
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				   sl::quaternion<T>&                                 q_star,
				   const std::size_t N,
				   const T* mass,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);

    T s_star;
    lsq_absolute_orientation_in(t_star, q_star, s_star, N, mass, a, b, false);
  }
  
  /** 
   *  Find least square optimal affine map such that b=s R(a)+t
   *  The map minimizes the error sum_i || (b_i-sR(a_i)-t) ||
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				   sl::quaternion<T>&                                 q_star,
				   T&                                                 s_star,
				   const std::size_t N,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);
    lsq_absolute_orientation_in(t_star, q_star, s_star, N, (const T*)0, a, b);
  }

  /** 
   *  Find weighted least square optimal affine map such that b= s R(a)+t
   *  The map minimizes the error sum_i ||m_i * (b_i-s R(a_i)-t) ||
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::affine_map<3,T>& tRs,
				   const std::size_t N,
				   const T* mass,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);

    sl::fixed_size_vector<sl::column_orientation,3,T> t_star;
    sl::quaternion<T>                                 q_star;
    T                                                 s_star;

    lsq_absolute_orientation_in(t_star, q_star, s_star, N, mass, a, b);

    tRs =
      sl::rigid_body_map<3,T>(t_star) *
      sl::rigid_body_map<3, T>(q_star) *
      sl::linear_map_factory<3,T>::scaling(s_star);
  }
  
  /** 
   *  Find least square optimal affine map such that b= s R(a)+t
   *  The map minimizes the error sum_i ||(b_i-s R(a_i)-t) ||
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::affine_map<3,T>& tRs,
				   const std::size_t N,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);

    lsq_absolute_orientation_in(tRs, N, (const T*)0, a, b);
  }

  /** 
   *  Find least square optimal rigid body map such that b= R(a)+t
   *  The map minimizes the error sum_i || (b_i-R(a_i)-t) ||
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::fixed_size_vector<sl::column_orientation,3,T>& t_star,
				   sl::quaternion<T>&                                 q_star,
				   const std::size_t N,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);

    lsq_absolute_orientation_in(t_star, q_star, N, (const T*)0, a, b);
  }

  /** 
   *  Find least square optimal rigid body map such that b= R(a)+t
   *  The map minimizes the error sum_i ||m_i * (b_i- R(a_i)-t) ||
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::rigid_body_map<3,T>& tR,
				   const std::size_t N,
				   const T* mass,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);

    sl::fixed_size_vector<sl::column_orientation,3,T> t_star;
    sl::quaternion<T>                                 q_star;

    lsq_absolute_orientation_in(t_star, q_star, N, mass, a, b);

    tR =
      sl::rigid_body_map<3,T>(t_star) *
      sl::rigid_body_map<3, T>(q_star);
  }

  /** 
   *  Find least square optimal rigid body map such that b= R(a)+t
   *  The map minimizes the error sum_i || (b_i- R(a_i)-t) ||
   */
  template <class T>
  void lsq_absolute_orientation_in(sl::rigid_body_map<3,T>& tR,
				   const std::size_t N,
				   const sl::fixed_size_point<3,T>* a,
				   const sl::fixed_size_point<3,T>* b) {
    assert(N>0);

    lsq_absolute_orientation_in(tR, N, (const T*)0, a, b);
  }
  
} // namespace sl

////////////////////////// OBSOLETE STUFF - BACKWARD COMPATIBILITY
namespace sl {

  /**
   *  ** OBSOLETE ** use lsq_absolute_orientation_in
   */
  template <class T>
  sl::rigid_body_map<3,T> absolute_orientation(const std::size_t N,
					       const T* mass,
					       const sl::fixed_size_point<3,T>* a,
					       const sl::fixed_size_point<3,T>* b) {
    //
    sl::rigid_body_map<3,T> Rt;
    lsq_absolute_orientation_in(Rt, N, mass, a, b);
    return Rt;
  }
  
  /** 
   *   ** OBSOLETE ** use lsq_absolute_orientation_in
   */
  template <class T>
  sl::rigid_body_map<3,T> absolute_orientation(const std::size_t N,
					       const sl::fixed_size_point<3,T>* a,
					       const sl::fixed_size_point<3,T>* b) {
    return absolute_orientation(N, (const T*)0, a, b);
  }

}

#endif
