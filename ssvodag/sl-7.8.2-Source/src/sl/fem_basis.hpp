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
#ifndef SL_FEM_BASIS_HPP
#define SL_FEM_BASIS_HPP

#include <sl/math.hpp>
#include <sl/cubature_rule.hpp>
#include <sl/dense_array.hpp> 
#include <sl/smart_pointer.hpp>
#include <sl/affine_map.hpp>
#include <sl/utility.hpp>

namespace sl { 

  /**
   * Function object representing G_scalar multivariate polynomials in
   * G_dimension variables. Quick'n'dirty implementation - should rewrite!
   */
  template <class G_scalar, std::size_t G_dimension>
  class fem_polynomial {
  public:
    enum { dimension = G_dimension };

    typedef G_scalar                                 value_t;
    typedef sl::fixed_size_point<dimension,value_t>  point_t;
  public:

    inline fem_polynomial() {
    }

    inline virtual ~fem_polynomial() {
    }

    /**
     *  The number d such as for each term the sum of the exponents sum
     *  to <= d
     */
    virtual std::size_t maximum_degree() const = 0;

    /// The value of the polynomial at the given parameter value
    virtual value_t operator()(const point_t&) const = 0;

  };
  
#define SL_DECLARE_FEM_POLYNOMIAL_2D(CLASS_NAME_, MAXDEGREE_, EXPR_) \
  template <class G_scalar> \
  class CLASS_NAME_ : public fem_polynomial<G_scalar, 2> { \
  public: \
    enum { dimension = 2 }; \
    typedef G_scalar                                 value_t; \
    typedef sl::fixed_size_point<dimension,value_t>  point_t; \
  public: \
    inline CLASS_NAME_() {} \
    virtual inline std::size_t maximum_degree() const { return MAXDEGREE_; } \
    virtual inline value_t operator()(const point_t& uv) const { { if (0 && uv[0]) {/*use variable*/} } return EXPR_ ; } \
  }

  // Constant to cubic orthonormal basis for the unit square [0,1]^2
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f0, 0, ((G_scalar)( 1.000000000000000)));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f1, 1, ((G_scalar)(-1.732050807568877) + (G_scalar)(  3.464101615137753) * uv[0]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f2, 1, ((G_scalar)(-1.732050807568877) + (G_scalar)(  3.464101615137753) * uv[1]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f3, 2, ((G_scalar)( 3.000000000000003) + (G_scalar)( -6.000000000000006) * uv[0] + (G_scalar)(-6.000000000000009) * uv[1] + (G_scalar)(12.000000000000021) * uv[0]*uv[1]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f4, 2, ((G_scalar)( 2.236067977499749) + (G_scalar)(-13.416407864998552) * uv[0] + (G_scalar)( 13.416407864998591) * uv[0]*uv[0]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f5, 2, ((G_scalar)( 2.236067977499781) + (G_scalar)(-13.416407864998723) * uv[1] + (G_scalar)( 13.416407864998760) * uv[1]*uv[1]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f6, 3, ((G_scalar)(-2.645751311064023) + (G_scalar)( 31.749015732770424) * uv[0] + (G_scalar)(-79.372539331927356) * uv[0]*uv[0] + (G_scalar)(52.915026221285316) * uv[0]*uv[0]*uv[0]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f7, 3, ((G_scalar)(-3.872983346207165) + (G_scalar)( 23.237900077242056) * uv[0] + (G_scalar)(  7.745966692414697) * uv[1] + (G_scalar)(-46.475800154488844) * uv[0]*uv[1] + (G_scalar)(-23.237900077239200) * uv[0]*uv[0] + (G_scalar)(46.475800154488617) * uv[0]*uv[0]*uv[1]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f8, 3, ((G_scalar)(-3.872983346207866) + (G_scalar)(  7.745966692416303) * uv[0] + (G_scalar)( 23.237900077246348) * uv[1] + (G_scalar)(-46.475800154495623) * uv[0]*uv[1] + (G_scalar)(-23.237900077245619) * uv[1]*uv[1] + (G_scalar)(46.475800154491409) * uv[0]*uv[1]*uv[1]));
  SL_DECLARE_FEM_POLYNOMIAL_2D(fem_polynomial_quad01_f9, 3, ((G_scalar)(-2.645751311064409) + (G_scalar)( 31.749015732781054) * uv[1] + (G_scalar)(-79.372539331951486) * uv[1]*uv[1] + (G_scalar)(52.915026221299712) * uv[1]*uv[1]*uv[1]));

  /**
   *  Basis for higher order approximations 
   *  of functions of G_dimension parameters.
   */
  template <class G_scalar, std::size_t G_dimension>
  class fem_basis {
  public:
    enum { dimension = G_dimension };

    typedef G_scalar                                  value_t;
    typedef sl::fixed_size_point<dimension,value_t>   point_t;

    typedef fem_polynomial<value_t, dimension>        polynomial_t;
    typedef cubature_rule_factory<value_t, dimension> cubarule_factory_t;
    typedef typename cubarule_factory_t::rule_t       rule_t;

  protected:
    std::vector<polynomial_t*> basis_array_;
    cubarule_factory_t*        cubarule_factory_;
    std::size_t                degree_;

    void update_degree() {
      degree_ = 0;
      for (std::size_t i = 0; i<basis_array_.size(); ++i) {
	degree_ = max(degree_, basis_array_[i]->maximum_degree());
      }
    }

  public:

    /**
     *  Construct a basis with the specified set of functions.
     *  The cubature rule factory crf is used by default to compute
     *  integrals over the standard domain.
     */
    inline fem_basis(const std::vector<polynomial_t*>& bf,
		     cubarule_factory_t*                 crf) 
      : basis_array_(bf), cubarule_factory_(crf)
    {
      SL_REQUIRE("Good function count", bf.size() >= 1);
      SL_REQUIRE("Cubature rule factory exists", crf);
      update_degree();
    }
      
    /**
     *  The number of functions in the basis
     */
    inline std::size_t count() const {
      return basis_array_.size();
    }

    /**
     *  The degree of the basis (i.e. the maximum degree of number of the basis functions)
     */
    inline std::size_t degree() const {
      return degree_;
    }

    /**
     *  The i-th basis function
     */
    inline const polynomial_t& operator()(std::size_t i) const {
      SL_REQUIRE("Good index", i<count());
      return *basis_array_[i];
    }

    /**
     *  The factory of cubature rules to compute integrals over the 
     *  standard domain.
     */
    inline const cubarule_factory_t& get_cubature_rule_factory() const {
      return *cubarule_factory_;
    }

    /**
     *  The approximated value of a function at position x given the
     *  set of coefficients c. 
     *  \latexonly
     *  The function is approximated by $\sum_\beta c_\beta \phi_\beta(x)$
     *  \endlatexonly
     */
    template <class G_tensor> 
    inline G_tensor value_at(const point_t& x,
			     const sized_raw_array_pointer<const G_tensor> c) const {
      SL_REQUIRE("Good coefficient count", c.count() == count());
      const std::size_t N = count();
      G_tensor result = this->operator()(0)(x) * c[0];
      for (std::size_t i=1; i<N; ++i) {
	result +=  this->operator()(i)(x) * c[i];
      }
      return result;
    }

  };
 
  /**
   *  Set h to the matrix of filter coefficients for push-pull operations 
   *  between a parent and a child in an element hierarchy.
   *  
   *  The filter coefficients are:
   *
   *  H_{\alpha,\beta}^{child} = {1\over A^{child}} \times
   *  \int _{A^{child}} \phi_\beta^{child}(x) \phi_\alpha(x) dA
   *  = \int _S phi_\alpha(u',v') phi_\beta(u,v) du dv
   *
   *  where S is the domain over which the child basis functions are
   *  defined and (u',v') is the point on the parent corresponding to
   *  the point (u,v) on the child.
   * 
   */
  template <class G_scalar, std::size_t G_dimension, class G_discriminant>
  void fem_push_pull_matrix_in(dense_array<G_scalar, 2, G_discriminant>& h,
			       const fem_basis<G_scalar, G_dimension>& parent_basis,
			       const fem_basis<G_scalar, G_dimension>& child_basis,
			       const affine_map<G_dimension, G_scalar>& child_to_parent) {
    SL_REQUIRE("Good matrix size", h.extent()(0) == parent_basis.count());
    SL_REQUIRE("Good matrix size", h.extent()(1) == child_basis.count());
    
    typedef G_scalar value_t;
    typedef sl::fixed_size_point<G_dimension,value_t>  point_t;

    for (std::size_t alpha=0; alpha<parent_basis.count(); ++alpha) {
      for (std::size_t beta=0; beta<child_basis.count(); ++beta) {
	// Choose a cubature rule for integrating parent_basis(alpha)*child_basis(beta)
	std::size_t product_degree = parent_basis(alpha).maximum_degree() + child_basis(beta).maximum_degree();
	const cubature_rule<G_scalar, G_dimension>& cr = parent_basis.get_cubature_rule_factory().rule_from_degree(product_degree);

#if 0
	std::cerr << "h(" << alpha << ", " << beta << ") = ";
#endif
	// Set h_a_b to \int _S phi_\alpha(u',v') phi_\beta(u,v) du dv
	value_t h_a_b = scalar_math<value_t>::zero();
	for (std::size_t k=0; k<cr.node_count(); ++k) {
	  const point_t& x_child  = cr.position(k);
	  const point_t  x_parent = transformation(child_to_parent,x_child);
	  h_a_b += 
	    cr.weight(k) * 
	    parent_basis(alpha)(x_parent) * 
	    child_basis(beta)(x_child);
#if 0
	  std::cerr << " + (" << 
	    cr.weight(k) << " * " <<
	    parent_basis(alpha)(x_parent) << " * " << 
	    child_basis(beta)(x_child) <<
	    ") ";
#endif
	}
	h(alpha,beta) = h_a_b;
#if 0
	std::cerr << " = " << h(alpha,beta) << std::endl;
#endif
      }
    }
  }
  
  /**
   *  Set h to the matrix of filter coefficients for pull operations 
   *  between a parent and a child in an element hierarchy in the case
   *  of potentially overlapping children and arbitrary shapes enclosed in
   *  regular domains (see Holzschuch (2000), Wavelet Radiosity on
   *  Arbitrary Planar Surfaces).
   *  The basic idea is to sample the parent and to use a characteristic
   *  function to restrict sample points to the actual geometry of the child.
   *
   *  Since the relative geometry of the parent and child is already
   *  taken into account by the characteristic function, the area 
   *  ratio should be set to one in the pull operations.
   *
   *  TEMPORARY IMPLEMENTATION, USES NON ADAPTIVE INTEGRATION RULE...
   */
  template <class G_scalar, std::size_t G_dimension, class G_discriminant, class G_characteristic>
  void fem_generalized_pull_matrix_in(dense_array<G_scalar, 2, G_discriminant>& h,
				      const fem_basis<G_scalar, G_dimension>& parent_basis,
				      const fem_basis<G_scalar, G_dimension>& child_basis,
				      const affine_map<G_dimension, G_scalar>& child_to_parent,
				      const G_characteristic& cf) {
    SL_REQUIRE("Good matrix size", h.extent()(0) == parent_basis.count());
    SL_REQUIRE("Good matrix size", h.extent()(1) == child_basis.count());

    // Choose a high precision cubature rule, since the characteristic function has
    // an unknown degree...
    const cubature_rule<G_scalar, G_dimension>& cr = 
      parent_basis.cubature_rule_factory().rule_from_degree(parent_basis.cubature_rule_factory().maximum_degree());

    typedef G_scalar value_t;
    typedef sl::fixed_size_point<G_dimension,value_t>  point_t;

    for (std::size_t alpha=0; alpha<parent_basis.count(); ++alpha) {
      for (std::size_t beta=0; beta<child_basis.count(); ++beta) {
	value_t h_a_b = scalar_math<value_t>::zero();
	for (std::size_t k=0; k<cr.node_count(); ++k) {
	  const point_t& x_parent = cr.position(k);
	  const point_t  x_child  = inverse_transformation(child_to_parent,x_parent);
	  if (cf(x_child)) {
	    h_a_b += 
	      cr.weight(k) * 
	      parent_basis(alpha)(x_parent) * 
	      child_basis(beta)(x_child);
	  }
	}
	h(alpha,beta) = h_a_b;
      }
    }
  }

  /**
   *  Accumulate push coefficients down in an element hierarchy from the
   *  parent element to the child element using the given
   *  push-pull matrix h.
   */
  template <class G_scalar, class G_tensor, class G_discriminant> 
  void fem_accumulate_push_in(sized_raw_array_pointer<G_tensor> child_coefficients,
			      sized_raw_array_pointer<const G_tensor> parent_coefficients,
			      const dense_array<G_scalar, 2, G_discriminant>& h) {
    SL_REQUIRE("Consistent parent", parent_coefficients.count() == h.extent()(0));
    SL_REQUIRE("Consistent child",  child_coefficients.count() == h.extent()(1));

    const std::size_t parent_basis_size = h.extent()(0);
    const std::size_t child_basis_size  = h.extent()(1);
   
    for (std::size_t beta=0; beta<child_basis_size; ++beta) {
      std::size_t alpha=0;
      child_coefficients[beta] += h(alpha, beta) * parent_coefficients[alpha];
      for (alpha=1; alpha<parent_basis_size; ++alpha) {
	child_coefficients[beta] += h(alpha, beta) * parent_coefficients[alpha];
      }
    }
  }

  /**
   *  Push coefficients down in an element hierarchy from the
   *  parent element to the child element using the given
   *  push-pull matrix h. This is equivalent to zeroing out
   *  the child coefficients before calling fem_accumulate_push_in
   */
  template <class G_scalar, class G_tensor, class G_discriminant> 
  void fem_push_in(sized_raw_array_pointer<G_tensor> child_coefficients,
		   sized_raw_array_pointer<const G_tensor> parent_coefficients,
		   const dense_array<G_scalar, 2, G_discriminant>& h) {
    SL_REQUIRE("Consistent parent", parent_coefficients.count() == h.extent()(0));
    SL_REQUIRE("Consistent child",  child_coefficients.count() == h.extent()(1));

    const std::size_t parent_basis_size = parent_coefficients.count();
    const std::size_t child_basis_size  = child_coefficients.count();
   
    for (std::size_t beta=0; beta<child_basis_size; ++beta) {
      std::size_t alpha=0;
      child_coefficients[beta] = h(alpha, beta) * parent_coefficients[alpha];
      for (alpha=1; alpha<parent_basis_size; ++alpha) {
	child_coefficients[beta] += h(alpha, beta) * parent_coefficients[alpha];
      }
    }
  }

  /**
   *  Accumulate pull coefficients up in an element hierarchy from the
   *  child element to the parent element using the given
   *  push-pull matrix h and area ratio A^child / A^parent.
   *
   *  The contribution of the child to the parent is:
   *
   *  pull^{child} = {A^{child} \over A^{parent}} \times
   *			\sum_\beta c_{\beta}^{child} H_{\alpha,\beta}^{child}
   */
  template <class G_scalar, class G_tensor, class G_discriminant> 
  void fem_accumulate_pull_in(sized_raw_array_pointer<G_tensor> parent_coefficients,
			      const dense_array<G_scalar, 2, G_discriminant>& h,
			      sized_raw_array_pointer<const G_tensor> child_coefficients,
			      G_scalar area_ratio = scalar_math<G_scalar>::one()) {
    SL_REQUIRE("Consistent parent", parent_coefficients.count() == h.extent()(0));
    SL_REQUIRE("Consistent child",  child_coefficients.count() == h.extent()(1));
    SL_REQUIRE("Non-negative area ratio", is_non_negative(area_ratio));

    const std::size_t parent_basis_size = h.extent()(0);
    const std::size_t child_basis_size  = h.extent()(1);
   
    for (std::size_t alpha=0; alpha<parent_basis_size; ++alpha) {
      std::size_t beta=0;
      parent_coefficients[alpha] += area_ratio * h(alpha, beta) * child_coefficients[beta];
      for (beta=1; beta<child_basis_size; ++beta) {
	parent_coefficients[alpha] += area_ratio * h(alpha, beta) * child_coefficients[beta];
      }
    }
  }

  /**
   *  Pull coefficients up in an element hierarchy from the
   *  child element to the parent element using the given
   *  push-pull matrix h and area ratio A^child / A^parent.
   *
   *  This is equivalent to zeroint out the parent coefficients
   *  and then call fem_accumulate_pull_in
   */
  template <class G_scalar, class G_tensor, class G_discriminant> 
  void fem_pull_in(sized_raw_array_pointer<G_tensor> parent_coefficients,
		   const dense_array<G_scalar, 2, G_discriminant>& h,
		   sized_raw_array_pointer<const G_tensor> child_coefficients,
		   G_scalar area_ratio = scalar_math<G_scalar>::one()) {
    SL_REQUIRE("Consistent parent", parent_coefficients.count() == h.extent()(0));
    SL_REQUIRE("Consistent child",  child_coefficients.count() == h.extent()(1));
    SL_REQUIRE("Non-negative area ratio", is_non_negative(area_ratio));

    const std::size_t parent_basis_size = h.extent()(0);
    const std::size_t child_basis_size  = h.extent()(1);
   
    for (std::size_t alpha=0; alpha<parent_basis_size; ++alpha) {
      std::size_t beta=0;
      parent_coefficients[alpha] = h(alpha, beta) * child_coefficients[beta];
      for (beta=1; beta<child_basis_size; ++beta) {
	parent_coefficients[alpha] += h(alpha, beta) * child_coefficients[beta];
      }
      parent_coefficients[alpha] *= area_ratio;
    }
  }
  
} // namespace sl

namespace sl {

  /// The type of supported basis functions
  typedef enum { 
    FEM_BASIS_CONSTANT,
    // Non-product basis functions
    FEM_BASIS_LINEAR,  
    FEM_BASIS_QUADRATIC,
    FEM_BASIS_CUBIC,
    // Product-Legendre basis functions
    FEM_BASIS_PRODUCT_LINEAR,
    FEM_BASIS_PRODUCT_QUADRATIC,
    FEM_BASIS_PRODUCT_CUBIC,
    // Size
    FEM_BASIS_FIRST_TYPE = FEM_BASIS_CONSTANT,
    FEM_BASIS_LAST_TYPE = FEM_BASIS_PRODUCT_CUBIC
  } fem_basis_type;
  
  /**
   *  Objects that create bases for higher order 
   *  approximations of functions of G_dimension parameters.
   */
  template <class G_scalar, std::size_t G_dimension>
  class fem_basis_factory {
  public:
    enum { dimension = G_dimension };
    
    typedef G_scalar                                value_t;
    typedef sl::fixed_size_point<dimension,value_t> point_t;
    
    typedef fem_basis<value_t, dimension>           basis_t;
    typedef fem_polynomial<value_t, dimension>      polynomial_t;

  protected:
   
    basis_t* basis_array_[std::size_t(FEM_BASIS_LAST_TYPE)+1];

  protected:

    /**
     *  Construct the object
     */
    inline fem_basis_factory() {
    }
 
  public:

    /**
     *  Destruct the object
     */
    virtual inline ~fem_basis_factory() {
      for (std::size_t i=0; i< std::size_t(FEM_BASIS_LAST_TYPE)+1; ++i) {
	delete basis_array_[i];
      }
    }

    /**
     *  Is tp a valid basis type?
     */
    inline bool good_basis_type(fem_basis_type tp) const {
      return 
	(tp >= FEM_BASIS_FIRST_TYPE) && 
	(tp <  FEM_BASIS_LAST_TYPE);
    }

    /**
     *  The number of basis functions (and thus of coefficients) for
     *  a basis of type tp
     */
    inline std::size_t basis_size(fem_basis_type tp) const {
      SL_REQUIRE("Good type", good_basis_type(tp));
      
      return basis_array_[tp] ? basis_array_[tp]->count() : 0;
    }

    /**
     *  A basis of type tp for the domain handled by this factory
     */
    inline const basis_t& basis(fem_basis_type tp) const {
      SL_REQUIRE("Good type", good_basis_type(tp));
      SL_REQUIRE("Basis exist", basis_size(tp) > 0);

      return *basis_array_[tp];
    }

    /// A new constant basis
    virtual basis_t* new_constant_basis() const = 0;
    /// A new linear basis
    virtual basis_t* new_linear_basis() const = 0;
    /// A new quadratic basis
    virtual basis_t* new_quadratic_basis() const = 0;
    /// A new cubic basis
    virtual basis_t* new_cubic_basis() const = 0;
    /// A new product linear basis
    virtual basis_t* new_product_linear_basis() const = 0;
    /// A new product quadratic basis
    virtual basis_t* new_product_quadratic_basis() const = 0;
    /// A new product cubic basis
    virtual basis_t* new_product_cubic_basis() const = 0;

  }; // fem_basis_factory

  template <fem_basis_type tp>
  struct quad01_fem_basis_coefficient_count {
    enum { value = 0 }; /// 0 = Unsupported basis
  };

  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_CONSTANT>  { enum { value = 1 }; };

  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_LINEAR>    { enum { value = 3 }; };
  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_QUADRATIC> { enum { value = 6 }; };
  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_CUBIC>     { enum { value = 10 }; };

  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_PRODUCT_LINEAR>    { enum { value = 4 }; };
  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_PRODUCT_QUADRATIC> { enum { value = 6 }; };
  template <> struct quad01_fem_basis_coefficient_count<FEM_BASIS_PRODUCT_CUBIC>     { enum { value = 10 }; };


  /**
   *  Objects that create bases for higher order 
   *  approximations of functions  defined on the unit
   *  square [0,1]^2
   */
  template <class G_scalar>
  class quad01_fem_basis_factory: public fem_basis_factory<G_scalar, 2> {
  public:
    enum { dimension = 2 };
    
    typedef G_scalar                                value_t;
    typedef sl::fixed_size_point<dimension,value_t> point_t;
    
    typedef fem_basis<value_t,dimension>            basis_t;
    typedef fem_polynomial<value_t, dimension>      polynomial_t;
    
  public:

    quad01_fem_basis_factory() {
      this->basis_array_[FEM_BASIS_CONSTANT] = new_constant_basis();
      // Non-product basis functions
      this->basis_array_[FEM_BASIS_LINEAR] = new_linear_basis();
      this->basis_array_[FEM_BASIS_QUADRATIC] = new_quadratic_basis();
      this->basis_array_[FEM_BASIS_CUBIC] = new_cubic_basis();
      // Product-Legendre basis functions
      this->basis_array_[FEM_BASIS_PRODUCT_LINEAR] = new_product_linear_basis();
      this->basis_array_[FEM_BASIS_PRODUCT_QUADRATIC] = new_product_quadratic_basis();
      this->basis_array_[FEM_BASIS_PRODUCT_CUBIC] = new_product_cubic_basis();
    }

    /// A new constant basis
    virtual basis_t* new_constant_basis() const {
      std::vector<polynomial_t*> f(1);
      f[0] = new fem_polynomial_quad01_f0<value_t>();
      return new basis_t(f, new quad01_cubature_rule_factory<value_t>());
    }

    /// A new linear basis
    virtual basis_t* new_linear_basis() const {
      std::vector<polynomial_t*> f(3);
      f[0] = new fem_polynomial_quad01_f0<value_t>();
      f[1] = new fem_polynomial_quad01_f1<value_t>();
      f[2] = new fem_polynomial_quad01_f2<value_t>();
      return new basis_t(f, new quad01_cubature_rule_factory<value_t>());
    }
      
    /// A new quadratic basis
    virtual basis_t* new_quadratic_basis() const {
      std::vector<polynomial_t*> f(6);
      f[0] = new fem_polynomial_quad01_f0<value_t>();
      f[1] = new fem_polynomial_quad01_f1<value_t>();
      f[2] = new fem_polynomial_quad01_f2<value_t>();
      f[3] = new fem_polynomial_quad01_f3<value_t>();
      f[4] = new fem_polynomial_quad01_f4<value_t>();
      f[5] = new fem_polynomial_quad01_f5<value_t>();
      return new basis_t(f, new quad01_cubature_rule_factory<value_t>());
    }
      
    /// A new cubic basis
    virtual basis_t* new_cubic_basis() const {
      std::vector<polynomial_t*> f(10);
      f[0] = new fem_polynomial_quad01_f0<value_t>();
      f[1] = new fem_polynomial_quad01_f1<value_t>();
      f[2] = new fem_polynomial_quad01_f2<value_t>();
      f[3] = new fem_polynomial_quad01_f3<value_t>();
      f[4] = new fem_polynomial_quad01_f4<value_t>();
      f[5] = new fem_polynomial_quad01_f5<value_t>();
      f[6] = new fem_polynomial_quad01_f6<value_t>();
      f[7] = new fem_polynomial_quad01_f7<value_t>();
      f[8] = new fem_polynomial_quad01_f8<value_t>();
      f[9] = new fem_polynomial_quad01_f9<value_t>();
      return new basis_t(f, new quad01_cubature_rule_factory<value_t>());
    }

    /// A new product linear basis
    virtual basis_t* new_product_linear_basis() const {
      std::vector<polynomial_t*> f(4);
      f[0] = new fem_polynomial_quad01_f0<value_t>();
      f[1] = new fem_polynomial_quad01_f1<value_t>();
      f[2] = new fem_polynomial_quad01_f2<value_t>();
      f[3] = new fem_polynomial_quad01_f3<value_t>();
      return new basis_t(f, new quad01_cubature_rule_factory<value_t>());
    }
      
    /// A new product quadratic basis
    virtual basis_t* new_product_quadratic_basis() const {
      SL_TRACE_OUT(1) << "Product-Legendre Quadratic Basis not implemented - using non-product one" << std::endl;
      return new_quadratic_basis();
    }
      
    /// A new product cubic basis
    virtual basis_t* new_product_cubic_basis() const {
      SL_TRACE_OUT(1) << "Product-Legendre Cubic Basis not implemented - using non-product one" << std::endl;
      return new_cubic_basis();
    }

  }; // quad01_fem_basis_factory
  

  // TODO: Basis functions for triangular elements

} // namespace sl

#endif
