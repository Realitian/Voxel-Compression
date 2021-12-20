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
#ifndef SL_FEM_COUPLING_HPP
#define SL_FEM_COUPLING_HPP

#include <sl/fem_basis.hpp>

namespace sl {

  /**
   * Function object representing a constant function returning the same 
   * G_scalar for each parameter of type point
   */
  template <class G_scalar, std::size_t G_dimension>
  class fem_constant {
  public:
    typedef G_scalar                                  value_t;
    typedef sl::fixed_size_point<G_dimension,value_t> arg1_t;
  protected:
    value_t result_;
  public:
    inline fem_constant(const value_t& r): result_(r) {}
    inline const value_t& operator()(const arg1_t&) const { return result_; }
  };

  /**
   *  Set F to the coupling coefficients between two elements src and destination:
   *
   *  F_{dst,\alpha;src,\beta} = \int _{A_{dst}} \phi_{dst,\alpha}(x) \int _{A_{src}} \phi_{src,\beta}(y) K(x,y) dA_{src} dA_{dst}
   *
   */
  template <class G_tensor, class G_scalar, std::size_t G_dimension, class G_discriminant, class G_kernel, class G_src_jacobian, class G_dst_jacobian>
  void fem_coupling_coefficients_in(dense_array<G_tensor, 2, G_discriminant>& F,
				    const G_kernel& K,
				    std::size_t     K_degree,
				    const fem_basis<G_scalar, G_dimension>& src_basis,
				    const G_src_jacobian& src_J,
				    const fem_basis<G_scalar, G_dimension>& dst_basis,
				    const G_dst_jacobian& dst_J) {
    SL_REQUIRE("Good F dimension", F.extent()(0) == dst_basis.count());
    SL_REQUIRE("Good F dimension", F.extent()(1) == src_basis.count());
    
    enum { dimension = G_dimension };
    typedef G_scalar                      value_t;
    typedef G_tensor                      tensor_t;
    typedef G_kernel                      kernel_t;
    typedef fem_basis<value_t, dimension> basis_t;
    typedef typename basis_t::rule_t      cubarule_t;
    typedef typename basis_t::point_t     point_t;

    // x / alpha / k --> dest.  
    // y / beta  / l --> source

    // Choose source and destination cubature rules
    const std::size_t src_degree = K_degree + src_basis.degree();
    const std::size_t dst_degree = dst_basis.degree() + src_degree;

    const cubarule_t& src_cubarule = src_basis.cubature_rule_factory().rule_from_degree(src_degree);
    const cubarule_t& dst_cubarule = dst_basis.cubature_rule_factory().rule_from_degree(dst_degree);
    
    const std::size_t N = dst_cubarule.node_count();
    const std::size_t M = src_cubarule.node_count();
    const std::size_t A = dst_basis.count();
    const std::size_t B = src_basis.count();
    
    // Compute K(xk,yl) for each pair of cubature points
    dense_array<tensor_t, 2, void> Kval(N,M);
    for (std::size_t k=0; k<N; ++k) {
      for (std::size_t l=0; l<M; ++l) {
	const point_t& x_k = dst_cubarule.position(k);
	const point_t& y_l = src_cubarule.position(l);
	Kval(k,l) = K(x_k,y_l);
      }
    }

    // Compute inner integral for each pair of dest. cubature point / src. basis
    dense_array<tensor_t, 2, void> Kinner(N,B);
    for (std::size_t k=0; k<N; ++k) {
      for (std::size_t beta = 0; beta<B; ++beta) {
	std::size_t l=0;
	const point_t& y_l = src_cubarule.position(l);
	tensor_t K_k_beta = 
	  src_cubarule.weight(l) *
	  src_basis(beta)(y_l) *
	  src_J(y_l) *
	  Kval(k,l);
	for (std::size_t l=1; l<M; ++l) {
	  const point_t& y_l = src_cubarule.position(l);
	  K_k_beta += 
	    src_cubarule.weight(l) *
	    src_basis(beta)(y_l) *
	    src_J(y_l) *
	    Kval(k,l);
	}
	Kinner(k,beta) = K_k_beta;
      }
    }

    // Compute outer integral, store in F
    for (std::size_t alpha=0; alpha<A; ++alpha) {
      for (std::size_t beta = 0; beta<B; ++beta) {
	std::size_t k=0;
	const point_t& x_k = dst_cubarule.position(k);
	tensor_t F_alpha_beta = 
	  dst_cubarule.weight(k) *
	  dst_basis(alpha)(x_k) *
	  dst_J(x_k) *
	  Kinner(k,beta);
	for (std::size_t k=1; k<N; ++k) {
	  const point_t& x_k = dst_cubarule.position(k);
	  F_alpha_beta += 
	    dst_cubarule.weight(k) *
	    dst_basis(alpha)(x_k) *
	    dst_J(x_k) *
	    Kinner(k,beta);
	}
	F(alpha,beta) = F_alpha_beta;
      }
    }
  }

}

#endif
