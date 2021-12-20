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
#ifndef SL_CUBATURE_RULE_HPP
#define SL_CUBATURE_RULE_HPP

#include <sl/math.hpp>
#include <sl/fixed_size_point.hpp>
#include <sl/std_serializer.hpp>
#include <vector>

namespace sl {

  /**
   *  Cubature rule for integrating scalar functions of
   *  G_dimension variables.
   */
  template <class G_scalar, std::size_t G_dimension>
  class cubature_rule {
  public:
    enum { dimension = G_dimension };
    
    typedef G_scalar                       value_t;
    typedef sl::fixed_size_point<dimension,value_t>  point_t;
    typedef std::pair<value_t, point_t>    node_t;

  protected:
    std::size_t         degree_;
    std::vector<node_t> nodes_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << degree_ << nodes_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> degree_ >> nodes_;
    }

  public:
    
    /**
     *  Construct a cubature rule with the given degree and nodes. 
     *  The nodes in a_nodes should be given so that the rule is 
     *  exact for polynomials in x1, x2, .. xn whose 
     *  exponents sum to d <= a_degree.
     */
    inline cubature_rule(const std::size_t a_degree,
			 const std::vector<node_t>& a_nodes):
      degree_(a_degree), nodes_(a_nodes) {
      SL_REQUIRE("Good degree", a_degree > 0);
      SL_REQUIRE("Good node count", a_nodes.size() > 0);
    }

    /**
     *  The degree of the rule
     */
    inline std::size_t degree() const {
      return degree_;
    }

    /**
     *  The number of nodes of the rule
     */
    inline std::size_t node_count() const { 
      return nodes_.size();
    }

    /**
     *  The i-th node of the rule
     */
    inline const node_t& node(std::size_t i) const {
      SL_REQUIRE("Good index", i<node_count());
      return nodes_[i];
    }

    /**
     *  The weight of the i-th node of the rule
     */
     inline const value_t& weight(std::size_t i) const {
      SL_REQUIRE("Good index", i<node_count());
      return nodes_[i].first;
    }

    /**
     *  The position of the i-th node of the rule
     */
     inline const point_t& position(std::size_t i) const {
      SL_REQUIRE("Good index", i<node_count());
      return nodes_[i].second;
    }

  }; // cubature_rule


  /**
   *  The value of the integral of f. G_tensort is the type of an
   *  object defining a copy constructor, the multiplication
   *  with a scalar and the sum with another tensor
   *  with the usual operator notation. G_tensor_function is the type
   *  of an object defining an operator() with parameter point_t
   *  returning a G_tensor.
   */
  template <class G_tensor, class G_scalar, std::size_t dimension, class G_tensor_function>
  inline G_tensor integral(const cubature_rule<G_scalar, dimension>& r, 
			   const G_tensor_function& f) {
    const std::size_t N = r.node_count();
    
    G_tensor result = r.weight(0)*f(r.position(0));
    for (std::size_t i=1; i<N; ++i) {
      result += r.weight(i)*f(r.position(i));
    }
    return result;
  }

} // namespace sl



namespace sl {

  /**
   *  Objects that create cubature rules for
   *  scalar functions of G_dimension variables.
   */
  template <class G_scalar, std::size_t G_dimension>
  class cubature_rule_factory {
  public:
    enum { dimension = G_dimension };
    
    typedef G_scalar                                 value_t;
    typedef sl::fixed_size_point<dimension, value_t> point_t;
    typedef std::pair<value_t, point_t>              node_t;

    typedef cubature_rule<value_t, dimension>        rule_t;

  protected:

    std::vector<rule_t> rules_;

  protected:

    inline cubature_rule_factory() {}
 
  public:

    virtual inline ~cubature_rule_factory() {
    }

    /**
     *  A cubature rule that integrates exactly polynomials of
     *  degree n.
     */
    const rule_t& rule_from_degree(std::size_t degree) const { 
      SL_REQUIRE("Has rules", rules_.size() > 0);
      std::size_t actual_degree = degree;
      if (actual_degree >= rules_.size()) {
	actual_degree = rules_.size() - 1;
	SL_TRACE_OUT(0) << 
	  "Selecting a rule of degree " << 
	  actual_degree << " to integrate polynomials of degree " <<
	  degree <<
	  std::endl;
      }

      const rule_t& result = rules_[actual_degree];
      SL_ENSURE("Consistent", result.degree() >= actual_degree);
      return result;
    }

    /**
     *  The highest degree of a rule returned by this factory
     */
    inline std::size_t maximum_degree() const {
      SL_REQUIRE("Has rules", rules_.size() > 0);
       return rules_.size()-1;
    }

    /**
     *  The lowest degree of a rule returned by this factory
     */
    inline std::size_t minimum_degree() const {
      SL_REQUIRE("Has rules", rules_.size() > 0);
      return rules_[0].degree();
    }

  };


  /**
   *  Objects that create cubature rules for
   *  scalar functions defined over the unit square [0,1]^2
   */
  template <class G_scalar>
  class quad01_cubature_rule_factory: public cubature_rule_factory<G_scalar, 2> {
  public:
    enum { dimension = 2 };
    
    typedef G_scalar                                value_t;
    typedef sl::fixed_size_point<dimension,value_t> point_t;
    typedef std::pair<value_t, point_t>             node_t;

    typedef cubature_rule<value_t, dimension>       rule_t;

  protected:

    /**
     *  Transform the nodes from the [-1,1]^N domain to the
     *  [0,1]^N domain.
     */
    static void remap_1_1_to_0_1(std::vector<node_t>& nodes) {
      static const value_t one      = scalar_math<value_t>::one();
      static const value_t one_half = reciprocal(scalar_math<value_t>::two());
      static const value_t scale    = one_half * one_half; // 1 / Area of -1..1
      for (std::size_t i=0; i<nodes.size(); ++i) {
	nodes[i].first *= scale;
	for (std::size_t j=0; j<dimension; ++j) {
	  nodes[i].second[j] = one_half * (nodes[i].second[j] + one);
	}
      }
    }

  public:

    /**
     *  Construct a factory of cubature rules for
     *  integrating scalar functions over the unit
     *  square
     */
    quad01_cubature_rule_factory() {
      (this->rules_).clear();
      (this->rules_).push_back(d1n1());  // degree = 1 also used for degree = 0
      (this->rules_).push_back(d1n1());  // degree = 1
      (this->rules_).push_back(d2n3());  // degree = 2
      (this->rules_).push_back(d3n4());  // degree = 3
      (this->rules_).push_back(d4n6());  // degree = 4
      (this->rules_).push_back(d5n7());  // degree = 5
      (this->rules_).push_back(d6n10()); // degree = 6
      (this->rules_).push_back(d7n12()); // degree = 7
      (this->rules_).push_back(d8n16()); // degree = 8
      (this->rules_).push_back(d9n17()); // degree = 9
    }

  public: 

    /**
     *  Cubature rule, Quad [0,1]^2, degree 1, 1 node
     */
    rule_t d1n1() const {
      std::vector<node_t> n(1);
      n[0]= node_t(value_t(4.0),point_t(0.0,0.0));
      remap_1_1_to_0_1(n);
      
      return rule_t(1,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 2, 3 nodes
     *  Stroud, "Approximate Calculation of Multiple Integrals", 1971
     */
    rule_t d2n3() const {
      static const value_t w = value_t(4.0/3.0);
      static const value_t u = value_t(std::sqrt(2./3.));
      static const value_t c = value_t(std::cos(2.0*sl::Pi(u)/3.0));
      static const value_t s = value_t(std::sin(2.0*sl::Pi(u)/3.0));
      static const value_t z = value_t(0.0);
      
      std::vector<node_t> n(3);
      n[0] = node_t(w, point_t(u  ,  z));
      n[1] = node_t(w, point_t(u*c,  u*s));
      n[2] = node_t(w, point_t(u*c, -u*s));
      remap_1_1_to_0_1(n);
		  
      return rule_t(2,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 3, 4 nodes
     *  Davis & Rabinowitz, Methods of Numerical Integration, 2nd edition 1984, p. 367.
     */
    rule_t d3n4() const {
      static const value_t w = value_t(1.0);
      static const value_t z = value_t(0.0);
      static const value_t u = value_t(std::sqrt(2./3.));
      
      std::vector<node_t> n(4);
      n[0] = node_t(w, point_t( u,  z));
      n[1] = node_t(w, point_t( z,  u));
      n[2] = node_t(w, point_t(-u,  z));
      n[3] = node_t(w, point_t( z, -u));
      remap_1_1_to_0_1(n);
		  
      return rule_t(3,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 3, 4 nodes
     *  Product Gauss-Legendre
     */
    rule_t d3n4pg() const {
      static const value_t w = value_t(1.0);
      static const value_t u = value_t(std::sqrt(1./3.));

      std::vector<node_t> n(4);
      n[0] = node_t(w, point_t( u,  u));
      n[1] = node_t(w, point_t( u, -u));
      n[2] = node_t(w, point_t(-u,  u));
      n[3] = node_t(w, point_t(-u, -u));
      remap_1_1_to_0_1(n);
      
      return rule_t(3,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 4, 6 nodes
     *  Schmid, "On Cubature Formulae with a Minimal Number of Knots", Numer. Math. Vol 31 (1978) p281
     */
    rule_t d4n6() const {
      std::vector<node_t> n(6);
      n[0] = node_t(value_t(1.286412084888852), point_t(value_t( 0.0              ), value_t(-0.356822089773090)));
      n[1] = node_t(value_t(0.491365692888926), point_t(value_t( 0.0              ), value_t( 0.934172358962716)));
      n[2] = node_t(value_t(0.761883709085613), point_t(value_t( 0.774596669241483), value_t( 0.390885162530071)));
      n[3] = node_t(value_t(0.761883709085613), point_t(value_t(-0.774596669241483), value_t( 0.390885162530071)));
      n[4] = node_t(value_t(0.349227402025498), point_t(value_t( 0.774596669241483), value_t(-0.852765377881771)));
      n[5] = node_t(value_t(0.349227402025498), point_t(value_t(-0.774596669241483), value_t(-0.852765377881771)));
      remap_1_1_to_0_1(n);
      
      return rule_t(4,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 5, 7 nodes
     *  Stroud, "Approximate Calculation of Multiple Integrals", 1971
     */
    rule_t d5n7() const {
      static const value_t w1 = value_t(8./7.);
      static const value_t w2 = value_t(5./9.);
      static const value_t w3 = value_t(20./63.);
      static const value_t r  = value_t(std::sqrt(14./15.));
      static const value_t s  = value_t(std::sqrt(1./3.));
      static const value_t t  = value_t(std::sqrt(3./5.));

      std::vector<node_t> n(7);
      n[0] = node_t(value_t(w1), point_t(value_t(0.), value_t(0.)));
      n[1] = node_t(value_t(w2), point_t(value_t( s), value_t( t)));
      n[2] = node_t(value_t(w2), point_t(value_t( s), value_t(-t)));
      n[3] = node_t(value_t(w2), point_t(value_t(-s), value_t( t)));
      n[4] = node_t(value_t(w2), point_t(value_t(-s), value_t(-t)));
      n[5] = node_t(value_t(w3), point_t(value_t( r), value_t(0.)));
      n[6] = node_t(value_t(w3), point_t(value_t(-r), value_t(0.)));
      remap_1_1_to_0_1(n);
      
      return rule_t(5,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 6, 10 nodes
     *  Wissman, Becker, "Partially Symmetric Cubature Formulas for Even 
     *  Degrees of Exactness", SIAM. J. Numer. Anal., Vol 23 nr 3 (1986), p 676
     */
    rule_t d6n10() const {
      std::vector<node_t> n(10);
      n[0] = node_t(value_t(0.392750590964348), point_t(value_t( 0.0              ), value_t( 0.86983337525005)));
      n[1] = node_t(value_t(0.754762881242610), point_t(value_t( 0.0              ), value_t(-0.479406351612111)));
      n[2] = node_t(value_t(0.206166050588279), point_t(value_t( 0.863742826346154), value_t( 0.802837516207657)));
      n[3] = node_t(value_t(0.206166050588279), point_t(value_t(-0.863742826346154), value_t( 0.802837516207657)));
      n[4] = node_t(value_t(0.689992138489864), point_t(value_t( 0.518690521392592), value_t( 0.262143665508058)));
      n[5] = node_t(value_t(0.689992138489864), point_t(value_t(-0.518690521392592), value_t( 0.262143665508058)));
      n[6] = node_t(value_t(0.260517488732317), point_t(value_t( 0.933972544972849), value_t(-0.363096583148066)));
      n[7] = node_t(value_t(0.260517488732317), point_t(value_t(-0.933972544972849), value_t(-0.363096583148066)));
      n[8] = node_t(value_t(0.269567586086061), point_t(value_t( 0.608977536016356), value_t(-0.896608632762453)));
      n[9] = node_t(value_t(0.269567586086061), point_t(value_t(-0.608977536016356), value_t(-0.896608632762453)));
      remap_1_1_to_0_1(n);
      
      return rule_t(6,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 7, 12 nodes
     *  Stroud, "Approximate Calculation of Multiple Integrals", 1971
     */
    rule_t d7n12() const {
      static const value_t r = value_t(std::sqrt(6./7.));
      static const value_t s = value_t(std::sqrt((114.-3.*std::sqrt(583.))/287.));
      static const value_t t = value_t(std::sqrt((114.+3.*std::sqrt(583.))/287.));
      static const value_t w1 = value_t(49./810.*4.);
      static const value_t w2 = value_t((178981.+2769.*std::sqrt(583.))/1888920.*4.);
      static const value_t w3 = value_t((178981.-2769.*std::sqrt(583.))/1888920.*4.);

      std::vector<node_t> n(12);
      n[0] = node_t(value_t(w1), point_t(value_t( r), value_t(0.)));
      n[1] = node_t(value_t(w1), point_t(value_t(-r), value_t(0.)));
      n[2] = node_t(value_t(w1), point_t(value_t(0.), value_t( r)));
      n[3] = node_t(value_t(w1), point_t(value_t(0.), value_t(-r)));
      n[4] = node_t(value_t(w2), point_t(value_t( s), value_t( s)));
      n[5] = node_t(value_t(w2), point_t(value_t( s), value_t(-s)));
      n[6] = node_t(value_t(w2), point_t(value_t(-s), value_t( s)));
      n[7] = node_t(value_t(w2), point_t(value_t(-s), value_t(-s)));
      n[8] = node_t(value_t(w3), point_t(value_t( t), value_t( t)));
      n[9] = node_t(value_t(w3), point_t(value_t( t), value_t(-t)));
      n[10]= node_t(value_t(w3), point_t(value_t(-t), value_t( t)));
      n[11]= node_t(value_t(w3), point_t(value_t(-t), value_t(-t)));
      remap_1_1_to_0_1(n);
      
      return rule_t(7,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 8, 16 nodes
     *  Wissman, Becker, "Partially Symmetric Cubature Formulas for Even 
     *  Degrees of Exactness", SIAM. J. Numer. Anal., Vol 23 nr 3 (1986), p 676
     */
    rule_t d8n16() const {
      std::vector<node_t> n(16);
      n[0] = node_t(value_t(0.450276776305590),point_t(value_t( 0.0              ), value_t( 0.659560131960342)));
      n[1] = node_t(value_t(0.166570426777813),point_t(value_t( 0.0              ), value_t(-0.949142923043125)));
      n[2] = node_t(value_t(0.098869459933431),point_t(value_t( 0.952509466071562), value_t( 0.765051819557684)));
      n[3] = node_t(value_t(0.098869459933431),point_t(value_t(-0.952509466071562), value_t( 0.765051819557684)));
      n[4] = node_t(value_t(0.153696747140812),point_t(value_t( 0.532327454074206), value_t( 0.936975981088416)));
      n[5] = node_t(value_t(0.153696747140812),point_t(value_t(-0.532327454074206), value_t( 0.936975981088416)));
      n[6] = node_t(value_t(0.396686976072903),point_t(value_t( 0.684736297951735), value_t( 0.333656717735747)));
      n[7] = node_t(value_t(0.396686976072903),point_t(value_t(-0.684736297951735), value_t( 0.333656717735747)));
      n[8] = node_t(value_t(0.352014367945695),point_t(value_t( 0.233143240801405), value_t(-0.079583272377397)));
      n[9] = node_t(value_t(0.352014367945695),point_t(value_t(-0.233143240801405), value_t(-0.079583272377397)));
      n[10]= node_t(value_t(0.189589054577798),point_t(value_t( 0.927683319306117), value_t(-0.272240080612534)));
      n[11]= node_t(value_t(0.189589054577798),point_t(value_t(-0.927683319306117), value_t(-0.272240080612534)));
      n[12]= node_t(value_t(0.375101001147587),point_t(value_t( 0.453120687403749), value_t(-0.613735353398028)));
      n[13]= node_t(value_t(0.375101001147587),point_t(value_t(-0.453120687403749), value_t(-0.613735353398028)));
      n[14]= node_t(value_t(0.125618791640072),point_t(value_t( 0.837503640422812), value_t(-0.888477650535971)));
      n[15]= node_t(value_t(0.125618791640072),point_t(value_t(-0.837503640422812), value_t(-0.888477650535971)));
      remap_1_1_to_0_1(n);
      
      return rule_t(8,n);
    }

    /**
     *  Cubature rule, Quad [0,1]^2, degree 9, 17 nodes
     * Moeller, "Kubaturformeln mit minimaler Knotenzahl, Numer. Math. 25, 185 (1976) 
     */
    rule_t d9n17() const {
      static const value_t  b1 = value_t(0.96884996636197772072);
      static const value_t  b2 = value_t(0.75027709997890053354);
      static const value_t  b3 = value_t(0.52373582021442933604);
      static const value_t  b4 = value_t(0.07620832819261717318);
      static const value_t  c1 = value_t(0.63068011973166885417);
      static const value_t  c2 = value_t(0.92796164595956966740);
      static const value_t  c3 = value_t(0.45333982113564719076);
      static const value_t  c4 = value_t(0.85261572933366230775);
      static const value_t  w0 = value_t(0.52674897119341563786);
      static const value_t  w1 = value_t(0.08887937817019870697);
      static const value_t  w2 = value_t(0.11209960212959648528);
      static const value_t  w3 = value_t(0.39828243926207009528);
      static const value_t  w4 = value_t(0.26905133763978080301);

      std::vector<node_t> n(17);
      n[0] = node_t(w0, point_t( 0., 0.));
      n[1] = node_t(w1, point_t( b1, c1));
      n[2] = node_t(w1, point_t(-b1,-c1));
      n[3] = node_t(w1, point_t(-c1, b1));
      n[4] = node_t(w1, point_t( c1,-b1));
      n[5] = node_t(w2, point_t( b2, c2));
      n[6] = node_t(w2, point_t(-b2,-c2));
      n[7] = node_t(w2, point_t(-c2, b2));
      n[8] = node_t(w2, point_t( c2,-b2));
      n[9] = node_t(w3, point_t( b3, c3));
      n[10]= node_t(w3, point_t(-b3,-c3));
      n[11]= node_t(w3, point_t(-c3, b3));
      n[12]= node_t(w3, point_t( c3,-b3));
      n[13]= node_t(w4, point_t( b4, c4));
      n[14]= node_t(w4, point_t(-b4,-c4));
      n[15]= node_t(w4, point_t(-c4, b4));
      n[16]= node_t(w4, point_t( c4,-b4));
      remap_1_1_to_0_1(n);
      
      return rule_t(9,n);
    }

  }; // quad01_cubature_rule_factory

} // namespace sl

#endif
