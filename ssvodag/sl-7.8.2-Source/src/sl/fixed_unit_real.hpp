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
#ifndef SL_FIXED_UNIT_REAL_HPP
#define SL_FIXED_UNIT_REAL_HPP

# ifdef _MSC_VER
#   pragma warning(disable:4146)
# endif //_MSC_VER

#include <sl/interval.hpp> 
#include <sl/numeric_traits.hpp>
#include <sl/integer.hpp>
#include <sl/math.hpp>
#include <stdio.h>

namespace sl {

  /**
   *  Real numbers represented by (at least) N_Bits quantities linearly 
   *  mapped to the range [0,+1] if B_Signed is false, and [-1,+1) if 
   *  B_Signed is true. Currently, only saturated arithmetic is 
   *  implemented.
   */
  template <std::size_t N_Bits, bool B_Signed>
  class fixed_unit_real {
  public: // Types 

    /// Is the number signed?
    static  const bool     is_signed = B_Signed;

    /// The internal representation of the mantissa
    typedef typename gen_if<B_Signed, typename int_t<N_Bits>::least,
				      typename uint_t<N_Bits>::least>::type mantissa_t; 
    
    /// The unsigned representation of the mantissa
    typedef typename uint_t<N_Bits>::least umantissa_t;

    /// The actual number of bits in the representation (greater or equal than the requested size)
    static  const std::size_t bits = sizeof(mantissa_t)*8;
 
    typedef fixed_unit_real<N_Bits, B_Signed> self_t;

  protected: // Data

    mantissa_t           mantissa_;

  public: // Serialization
    
    inline void store_to(output_serializer& s) const {
      s << mantissa_;
    }
    
    inline void retrieve_from(input_serializer& s) {
      s >> mantissa_;
    }

  protected: // Helpers

    static inline mantissa_t mantissa_msb() {
      return (std::size_t)(1) << (bits-1);
    }

    static inline mantissa_t mantissa_max() {
      return sl::scalar_math<mantissa_t>::finite_upper_bound();
    }

    static inline mantissa_t mantissa_min() {
      return (is_signed) ? mantissa_t(1 << (bits-1)) : mantissa_t(0);
    }
    
    static inline double double_scale_factor() {
      return (is_signed) ? ((double)(umantissa_t)(1 << (bits-1))) : ((double)(mantissa_max()));
    }

    /// The upper bound (double) of this representation
    static inline double double_max() {
      if (is_signed) {
	return 1.0 - 0.5 / double_scale_factor();
      } else {
	return 1.0;
      }
    }

    /// The lower bound (double) of this representation
    static inline double double_min() {
      if (is_signed) {
	return -1.0;
      } else {
	return 0.0;
      }
    }

    static inline mantissa_t saturated(double d) {
      if (d >= double_max()) {
	return mantissa_max();
      } else if (d <= double_min()) {
	return mantissa_min();
      } else {
	double dd = d * double_scale_factor();
	return (dd >= 0)? mantissa_t(dd + 0.5) : mantissa_t(dd - 0.5);
      }
    }
    
  public: // Creation, Copy & Destruction

    /// Explicit init from internal representation v.
    inline explicit fixed_unit_real(const mantissa_t& v = 0) {
      mantissa_ = v;
    }

    /// Explicit init from double precision number v
    inline explicit fixed_unit_real(double v) {
      mantissa_ = saturated(v);
    }

    /// Set this to other
    inline fixed_unit_real(const self_t& other) {
      mantissa_  = other.mantissa();
    }

  public: // Access
    
   inline const mantissa_t& mantissa() const {
      return mantissa_;
    }

    /// The double precision representation 
    inline operator double() const {
      return (1.0 / double_scale_factor()) * mantissa_;
    }

    /// The double precision representation
    inline double value() const {
      return (1.0 / double_scale_factor()) * mantissa_;
    }

  public: // Assignment

    /// Set this to other
    inline self_t& operator= (const self_t& other) {
      mantissa_  = other.mantissa();
      return *this;
    }

    /// Set this to other
    inline self_t& operator= (double v) {
      mantissa_  = saturated(v);
      return *this;
    }

    /// Set this to the upper bound of v
    inline void set_to_sup(double v) {
      mantissa_  = saturated(v);
      if (value() < v) {
	mantissa_t m2 = mantissa_ + mantissa_t(1);
	if (m2 > mantissa_) mantissa_ = m2;
      }
    }
  
    /// Set this to the lower bound of v
    inline void set_to_inf(double v) {
      mantissa_  = saturated(v);
      if (value() > v) {
	mantissa_t m2 = mantissa_ - mantissa_t(1);
	if (m2 < mantissa_) mantissa_ = m2;
      }
    }
    
  public: // Comparison

    /** 
     * is this equal to other?
     */
    inline bool operator== (const self_t& other) const {
      return mantissa_ == other.mantissa_;
    }

    /** 
     * is this less than other?
     */
    inline bool operator< (const self_t& other) const {
      return mantissa_ < other.mantissa_;
    }

    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);
    
  public: // Arithmetic operations

    /// A copy of this
    inline self_t operator+() const {
      return *this;
    }

    /// Zero minus this
    inline self_t operator-() const {
      return (is_signed ? self_t(-mantissa_) : self_t(mantissa_t(0)));
    }

    /// The sum of this and other
    inline self_t operator+(const self_t& other) {
      mantissa_t result = this->mantissa() + other.mantissa();
      if (is_signed) {
	// Signed saturation
	if ((this->mantissa() ^ result) & (other.mantissa() ^ result) & mantissa_msb()) {
	  result = (result>0) ? mantissa_min() : mantissa_max();
	}
      } else {
	// Unsigned saturation
	if ((result <= this->mantissa()) && (other.mantissa() > 0)) {
	  result = mantissa_max();
	}
      }
      return self_t(result);
    }

    /// The subtraction of this and other
    inline self_t operator-(const self_t& other) {
      mantissa_t result = this->mantissa() - other.mantissa();
      if (is_signed) {
	// Signed saturation
	if ((this->mantissa() ^ result) & (-other.mantissa() ^ result) & mantissa_msb()) {
	  result = (result>0) ? mantissa_min() : mantissa_max();
	}
      } else {
	// Unsigned saturation
	if ((this->mantissa() < other.mantissa())) {
	  result = mantissa_min();
	}
      }
      return self_t(result);
    }
    
    /// The multiplication of this and other
    inline self_t operator*(const self_t& other) {
      /// Currently uses floating point, rewrite!
      return self_t(double(*this) * double(other));
    }

    /// The division of this and other
    inline self_t operator/(const self_t& other) {
      /// Currently uses floating point, rewrite!
      return self_t(double(*this) / double(other));
    }
    
    /// Add other to this
    inline self_t& operator+= (const self_t& other) {
      *this = *this + other;
      return *this;
    }

    /// Subtract other from this
    inline self_t& operator-= (const self_t& other) {
      *this = *this - other;
      return *this;
    }

    /// Multiply this by other
    inline self_t& operator*= (const self_t& other) {
      *this = *this * other;
      return *this;
    }
    
    /// Divide this by other
    inline self_t& operator/= (const self_t& other) {
      *this = *this / other;
      return *this;
    }
      
  };

};

// I/O

template <std::size_t N_Bits, bool B_Signed>
std::ostream& operator<< (std::ostream & os, const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
  typedef typename  sl::fixed_unit_real<N_Bits,B_Signed>::umantissa_t umantissa_t;
  os << "( " << intmax_t(x.mantissa()) << " / " << uintmax_t((umantissa_t)(1 << (x.bits-1))) << " )";
  return os;
}

namespace std {

  /// numeric limits for a fixed unit real
  template<std::size_t N_Bits, bool B_Signed> class numeric_limits< sl::fixed_unit_real<N_Bits,B_Signed> > {
    typedef sl::fixed_unit_real<N_Bits,B_Signed> fixreal_t;
    typedef typename fixreal_t::mantissa_t       mantissa_t;
    
  public:
    static const bool is_specialized = true;

    static fixreal_t min() throw() { return fixreal_t(mantissa_t(1)); }
    static fixreal_t max() throw() { return fixreal_t(numeric_limits<mantissa_t>::max()); }

    static const int digits = numeric_limits<mantissa_t>::digits;
    static const int digits10 = numeric_limits<mantissa_t>::digits10;
    static const bool is_signed = numeric_limits<mantissa_t>::is_signed;
    static const bool is_integer = false;
    static const bool is_exact = true;
    static const int radix =  2;
    static fixreal_t epsilon() throw() { return min(); }
    static fixreal_t round_error() throw() { return min(); }

    static const int min_exponent = 0; // ERROR
    static const int min_exponent10 = 0; // ERROR 
    static const int max_exponent = 0; 
    static const int max_exponent10 = 0; 

    static const bool has_infinity = false;
    static const bool has_quiet_NaN =  false;
    static const bool has_signaling_NaN = false;
    static const float_denorm_style has_denorm = (float_denorm_style)numeric_limits<mantissa_t>::has_denorm;
    static const bool has_denorm_loss = false;

    static fixreal_t infinity() throw() { return fixreal_t(numeric_limits<mantissa_t>::infinity()); }
    static fixreal_t quiet_NaN() throw() { return fixreal_t(numeric_limits<mantissa_t>::quiet_NaN()); }
    static fixreal_t signaling_NaN() throw() { return fixreal_t(numeric_limits<mantissa_t>::signaling_NaN()); }
    static fixreal_t denorm_min() throw() { return fixreal_t(numeric_limits<mantissa_t>::denorm_min()); }
 
    static const bool is_iec559 = false;
    static const bool is_bounded = true;
    static const bool is_modulo = false;
    
    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = (float_round_style)numeric_limits<double>::round_style;
  };

};

namespace sl {

  /// numeric traits for a fixed unit real
  template <std::size_t N_Bits, bool B_Signed>
  class numeric_traits< fixed_unit_real<N_Bits,B_Signed> > {
    typedef fixed_unit_real<N_Bits,B_Signed> fixed_unit_real_t;
    typedef typename fixed_unit_real<N_Bits,B_Signed>::mantissa_t mantissa_t;
  public:
    static std::string what() { 
      static char buf[16];
      sprintf(buf,"%d", (int)N_Bits);
      return std::string("fixed_unit_real< ")+
	std::string(buf)+ 
	std::string(",") + 
	std::string(B_Signed ? "true" : "false")+
	std::string(" >"); 
    }
    typedef fixed_unit_real< N_Bits, B_Signed >  T_sumtype; 
    typedef fixed_unit_real< N_Bits, B_Signed >  T_difftype;
    typedef fixed_unit_real< N_Bits, B_Signed >  T_floattype;
    typedef fixed_unit_real< N_Bits, true >      T_signedtype;
    enum { has_trivial_constructor = 0 }; 
    enum { is_specialized = 1 };   
    enum { precision_rank = 10+numeric_traits<mantissa_t>::precision_rank }; 
  };

};
    
namespace sl {

  /// constants and numeric functions of fixed unit real
  template<std::size_t N_Bits, bool B_Signed>
  class scalar_math< fixed_unit_real<N_Bits,B_Signed> >: public scalar_math_base< fixed_unit_real<N_Bits,B_Signed> >  {
    typedef fixed_unit_real<N_Bits,B_Signed> fixed_unit_real_t;
  public:
    static inline fixed_unit_real_t zero             () { return fixed_unit_real_t(); }
    static inline fixed_unit_real_t one              () { return fixed_unit_real_t(scalar_math<double>::one()); }
    static inline fixed_unit_real_t two              () { return fixed_unit_real_t(scalar_math<double>::two()); }
    static inline fixed_unit_real_t E                () { return fixed_unit_real_t(sl::E(double())); }
    static inline fixed_unit_real_t Log2_E           () { return fixed_unit_real_t(sl::Log2_E(double())); }
    static inline fixed_unit_real_t Log10_E          () { return fixed_unit_real_t(sl::Log10_E(double())); }
    static inline fixed_unit_real_t Ln_2             () { return fixed_unit_real_t(sl::Ln_2(double())); }
    static inline fixed_unit_real_t Ln_10            () { return fixed_unit_real_t(sl::Ln_10(double())); }
    static inline fixed_unit_real_t Pi               () { return fixed_unit_real_t(sl::Pi(double())); }
    static inline fixed_unit_real_t Pi_Over_2        () { return fixed_unit_real_t(sl::Pi_Over_2(double())); }
    static inline fixed_unit_real_t Pi_Over_4        () { return fixed_unit_real_t(sl::Pi_Over_4(double())); }
    static inline fixed_unit_real_t One_0ver_Pi      () { return fixed_unit_real_t(sl::One_Over_Pi(double())); }
    static inline fixed_unit_real_t Two_0ver_Pi      () { return fixed_unit_real_t(sl::Two_Over_Pi(double())); }
    static inline fixed_unit_real_t Two_Over_Sqrt_Pi () { return fixed_unit_real_t(sl::Two_Over_Sqrt_Pi(double())); }
    static inline fixed_unit_real_t Sqrt_2           () { return fixed_unit_real_t(sl::Sqrt_2(double())); }
    static inline fixed_unit_real_t One_Over_Sqrt_2  () { return fixed_unit_real_t(sl::One_Over_Sqrt_2(double())); }
  };

} // namespace sl


// basic functions of bounded scalars
namespace sl {

  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> sign(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(sign(x.mantissa()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> non_zero_sign(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(non_zero_sign(x.mantissa()));
  }
    
  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> round_down(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(round_down(x.mantissa()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> round_up(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(round_up(x.mantissa()));
  }
  
  template <std::size_t N_Bits, bool B_Signed> 
  static inline fixed_unit_real<N_Bits,B_Signed> pred(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(pred(x.mantissa()));
  }
  
  template <std::size_t N_Bits, bool B_Signed> 
  static inline fixed_unit_real<N_Bits,B_Signed> succ(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(succ(x.mantissa()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> abs(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(abs(x.mantissa()));
  }
  
  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> integer_round_towards_zero(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(integer_round_towards_zero(x.value()));
  }
  
  template <std::size_t N_Bits, bool B_Signed>
  inline fixed_unit_real<N_Bits,B_Signed> mod(const fixed_unit_real<N_Bits,B_Signed>& x, const fixed_unit_real<N_Bits,B_Signed>& y) {
    return fixed_unit_real<N_Bits,B_Signed>(mod(x.value(),y.value()));
  }
  
  template <std::size_t N_Bits, bool B_Signed> inline fixed_unit_real<N_Bits,B_Signed> sqr(const fixed_unit_real<N_Bits,B_Signed>& x) {
    return fixed_unit_real<N_Bits,B_Signed>(sqr(x.value()));
  }
}

// standard basic functions of bounded scalars
namespace std {

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> exp(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(exp(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> log(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(log(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> log10(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(log10(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> pow(const sl::fixed_unit_real<N_Bits,B_Signed>& x, int n) { 
    return sl::fixed_unit_real<N_Bits,B_Signed>(pow(x.value(),n));
  }    
	
  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> pow(const sl::fixed_unit_real<N_Bits,B_Signed>& x, 
											 const sl::fixed_unit_real<N_Bits,B_Signed>& y) { 
    return sl::fixed_unit_real<N_Bits,B_Signed>(pow(x.value(),y.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> ceil(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(ceil(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> floor(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(floor(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> sqrt(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(sqrt(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> atan(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(atan(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> acos(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(acos(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> asin(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(asin(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed>
  inline sl::fixed_unit_real<N_Bits,B_Signed> atan2(const sl::fixed_unit_real<N_Bits,B_Signed>& y, 
						    const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return 
      sl::fixed_unit_real<N_Bits,B_Signed>(atan2(y.value(),x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> sin(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(sin(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> cos(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(cos(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> tan(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(tan(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> cosh(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(cosh(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> sinh(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(sinh(x.value()));
  }

  template <std::size_t N_Bits, bool B_Signed> inline sl::fixed_unit_real<N_Bits,B_Signed> tanh(const sl::fixed_unit_real<N_Bits,B_Signed>& x) {
    return sl::fixed_unit_real<N_Bits,B_Signed>(tanh(x.value()));
  }
} // namespace std


// Useful typedefs
namespace sl {

  typedef fixed_unit_real< 8, true> fix8_t;
  typedef fixed_unit_real<16, true> fix16_t;
  typedef fixed_unit_real<32, true> fix32_t;

  typedef fixed_unit_real< 8,false> ufix8_t;
  typedef fixed_unit_real<16,false> ufix16_t;
  typedef fixed_unit_real<32,false> ufix32_t;

} // namespace sl

#endif






