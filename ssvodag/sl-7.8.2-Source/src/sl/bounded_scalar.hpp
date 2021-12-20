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
#ifndef SL_BOUNDED_SCALAR_HPP
#define SL_BOUNDED_SCALAR_HPP

#include <sl/interval.hpp> 
#include <sl/numeric_traits.hpp>
#include <sl/utility.hpp>

namespace sl {

  /// Bounded scalars, represented by an expected values with associated scalar bounds
  template <class T>
  class bounded_scalar {

  public: // Types 

    typedef T value_t;
    typedef bounded_scalar<value_t> self_t;

  protected: // Data

    value_t           value_;
    interval<value_t> bounds_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << value_ << bounds_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> value_ << bounds_;
    }

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);

  public: // Invariant

    /// is this in a coherent state?
    bool invariant() {
      bool result = bounds_.contains(value_);
#ifndef NDEBUG
      if (!result) {
        SL_TRACE_OUT(0) << "Invariant failed: " << value_ << " not in " << bounds_ << std::endl;
      }
#endif
      return result;
    }

  public: // Creation, Copy & Destruction
    
    /// Default init (zero)
    inline bounded_scalar()
        :
        value_(scalar_math<value_t>::zero()), bounds_(scalar_math<value_t>::zero(), scalar_math<value_t>::zero()) {
      SL_INVARIANT(invariant());
    }

    /// Explicit init from scalar value v
    inline bounded_scalar(const value_t& v)
        :
        value_(v), bounds_(v,v) {
      SL_INVARIANT(invariant());
    }

    /// Explicit init from scalar value v and bound b. If fix is true, b is "adjusted" to force containment of v.
    inline bounded_scalar(const value_t& v, const interval<T>& b, bool fix = true)
        :
        value_(v), bounds_(b) {
      SL_REQUIRE("Inclusion", fix || b.contains(v));
      if (fix) bounds_.merge(value_);
      SL_INVARIANT(invariant());
    }

#ifdef SL_OBSOLETE
    /// Explicit init from scalar value v and bound b. If fix is true, b is "adjusted" to force containment of v.
    inline bounded_scalar(const value_t& vmin, const value_t& v, const value_t& vmax, bool fix=true)
      :
        value_(v), bounds_(min(vmin,vmax), max(vmin,vmax)) {
      if (fix) bounds_.merge(value_);
      SL_INVARIANT(invariant());
    }
#endif
    
    /**
     * Constructor setting the bounded scalar from a manifest constant.
     * The range of the interval bound is given by the number of digits
     * of precision in the specification. For example, 
     * bounded_scalar("0.314") = [0.314-0.0005, 0.314, 0.314+0.0005]
     */
    bounded_scalar(const std::string& s) {
      SL_REQUIRE("String exists", !s.empty());
      bounds_ = s;
      value_ = bounds_.center();
      SL_INVARIANT(invariant());
    }

    /// Set this to other
    inline bounded_scalar(const self_t& other) {
      value_  = other.value();
      bounds_ = other.bounds();
      SL_INVARIANT(invariant());
    }

    /// Set this to other converting from another numeric representation
    template <class OTHER_T>
    explicit bounded_scalar(const bounded_scalar<OTHER_T>& other) {
      value_  = numeric_cast<T>(other.value());
      bounds_ = interval<T>(other.bounds());
      bounds_.merge(value_);
      SL_INVARIANT(invariant());
    }

  public: // Access

    /// The expected value 
    inline value_t value() const {
      return value_;
    }

    /// The expected value
    inline operator value_t() const {
      return value_;
    }

    /// The bounds
    inline const interval<T>& bounds() const {
      return bounds_;
    }

    /// The lower bound
    inline const value_t& lower_bound() const {
      return bounds_.lower_bound();
    }

    /// The upper bound
    inline const value_t& upper_bound() const {
      return bounds_.upper_bound();
    }

    /// The width of the associated interval
    inline value_t width() const {
      return bounds_.width();
    }

  public: // Comparison

    /** 
     * Is this equal to other?
     */
    inline bool operator== (const self_t& other) const {
      return value_ == other.value_ && bounds_ == other.bounds_;
    }

    /** 
     * Is this less than other?
     */
    inline bool operator< (const self_t& other) const {
      return
        (value_ < other.value_) ||
        ((value_ == other.value()) &&
         ((lower_bound() < other.lower_bound()) ||
          ((lower_bound() == other.lower_bound()) && (upper_bound() < other.upper_bound()))));
    }

    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);
    
  public: // Arithmetic operations

    /// Add other to this
    inline self_t& operator+= (const self_t& other) {
      value_  += other.value_;
      bounds_ += other.bounds_;
      bounds_.merge(value_); // Needed because of optimizer behavior 
      SL_INVARIANT(invariant());
      return *this;
    }

    /// Subtract other from this
    inline self_t& operator-= (const self_t& other) {
      value_  -= other.value_;
      bounds_ -= other.bounds_;
      bounds_.merge(value_); // Needed because of optimizer behavior 
      SL_INVARIANT(invariant());
      return *this;
    }

    /// Multiply this by other
    inline self_t& operator*= (const self_t& other) {
      value_  *= other.value_;
      bounds_ *= other.bounds_;
      bounds_.merge(value_); // Needed because of optimizer behavior 
      SL_INVARIANT(invariant());
      return *this;
    }
    
    /// Divide this by other
    inline self_t& operator/= (const self_t& other) {
      value_  /= other.value_;
      bounds_ /= other.bounds_;
      bounds_.merge(value_); // Needed because of optimizer behavior 
      SL_INVARIANT(invariant());
      return *this;
    }

    SL_OP_FIELD(self_t); // Without type promotion

    /// 0 - this
    inline self_t operator- () const {
      return zero(self_t()) - *this;
    }

    /// the result of the multiplication of this by other
    template <class T2>
    inline bounded_scalar<SL_PROMOTENAME(value_t,T2)> operator*(const bounded_scalar<T2>& y) const { 
      bounded_scalar<SL_PROMOTENAME(value_t,T2)> x(*this); 
      bounded_scalar<SL_PROMOTENAME(value_t,T2)> y2(y); 
      return x *= y2; 
    }

    /// the result of the division of this by other
    template <class T2>
    inline bounded_scalar<SL_PROMOTENAME(SL_FLOATTYPENAME(value_t),SL_FLOATTYPENAME(T2))> operator/(const bounded_scalar<T2>& y) const { 
      bounded_scalar<SL_PROMOTENAME(SL_FLOATTYPENAME(value_t),SL_FLOATTYPENAME(T2))> x(*this); 
      bounded_scalar<SL_PROMOTENAME(SL_FLOATTYPENAME(value_t),SL_FLOATTYPENAME(T2))> y2(y); 
      return x /= y2; 
    }

    /// the result of the sum of this with other
    template <class T2>
    inline bounded_scalar<SL_PROMOTENAME(value_t,T2)> operator+(const bounded_scalar<T2>& y) const { 
      bounded_scalar<SL_PROMOTENAME(value_t,T2)> x(*this); 
      bounded_scalar<SL_PROMOTENAME(value_t,T2)> y2(y); 
      return x += y2; 
    }
    
    /// the result of the subtraction of this with other
    template <class T2>
    inline bounded_scalar<SL_PROMOTENAME(value_t,T2)> operator-(const bounded_scalar<T2>& y) const { 
      bounded_scalar<SL_PROMOTENAME(value_t,T2)> x(*this); 
      bounded_scalar<SL_PROMOTENAME(value_t,T2)> y2(y); 
      return x -= y2; 
    }

    /// linear interpolation
    template <class T_PARAMETER>
    inline self_t lerp(const self_t& other, const T_PARAMETER& t) const {
      return *this + self_t((other - *this) * t);
    }
      
  };

};

// I/O

template <class T>
std::ostream& operator<< (std::ostream & os, const sl::bounded_scalar<T>& x) {
  os << "< " << x.value() << " " << x.bounds() << " >";
  return os;
}

namespace std {

  /// numeric limits for a bounded scalar
  template<class T> class numeric_limits< sl::bounded_scalar<T> > {
  public:
    static const bool is_specialized = true;

    static sl::bounded_scalar<T> min() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::min()); }
    static sl::bounded_scalar<T> max() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::max()); }

    static const int digits = numeric_limits<T>::digits;
    static const int digits10 = numeric_limits<T>::digits10;
    static const bool is_signed = numeric_limits<T>::is_signed;
    static const bool is_integer = numeric_limits<T>::is_integer;
    static const bool is_exact = numeric_limits<T>::is_exact;
    static const int radix =  numeric_limits<T>::radix;
    static sl::bounded_scalar<T> epsilon() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::epsilon()); }
    static sl::bounded_scalar<T> round_error() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::round_error()); }

    static const int min_exponent = numeric_limits<T>::min_exponent;
    static const int min_exponent10 = numeric_limits<T>::min_exponent10;
    static const int max_exponent = numeric_limits<T>::max_exponent;
    static const int max_exponent10 = numeric_limits<T>::max_exponent10;

    static const bool has_infinity = numeric_limits<T>::has_infinity;
    static const bool has_quiet_NaN =  numeric_limits<T>::has_quiet_NaN;
    static const bool has_signaling_NaN = numeric_limits<T>::has_signaling_NaN;
    static const float_denorm_style has_denorm = /*(float_denorm_style)*/numeric_limits<T>::has_denorm;
    static const bool has_denorm_loss = numeric_limits<T>::has_denorm_loss;

    static sl::bounded_scalar<T> infinity() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::infinity()); }
    static sl::bounded_scalar<T> quiet_NaN() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::quiet_NaN()); }
    static sl::bounded_scalar<T> signaling_NaN() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::signaling_NaN()); }
    static sl::bounded_scalar<T> denorm_min() throw() { return sl::bounded_scalar<T>(numeric_limits<T>::denorm_min()); }
 
    static const bool is_iec559 = numeric_limits<T>::is_iec559;
    static const bool is_bounded = numeric_limits<T>::is_bounded;
    static const bool is_modulo = numeric_limits<T>::is_modulo;
    
    static const bool traps = numeric_limits<T>::traps;
    static const bool tinyness_before = numeric_limits<T>::tinyness_before;
    static const float_round_style round_style = /*(float_round_style)*/numeric_limits<T>::round_style;
  };

};

namespace sl {

  /// numeric traits for a bounded scalar
  template <class T>
  class numeric_traits< bounded_scalar<T> > {
  public:
    static std::string what() { return std::string("bounded_scalar< ")+numeric_traits<T>::what()+std::string(" >"); }
    typedef bounded_scalar< typename numeric_traits<T>::T_sumtype >    T_sumtype; 
    typedef bounded_scalar< typename numeric_traits<T>::T_difftype >   T_difftype;
    typedef bounded_scalar< typename numeric_traits<T>::T_floattype >  T_floattype;
    typedef bounded_scalar< typename numeric_traits<T>::T_signedtype > T_signedtype;
    enum { has_trivial_constructor = 0 }; 
    enum { is_specialized = 1 };   
    enum { precision_rank = 10000+numeric_traits<T>::precision_rank }; 
  };

};
    
namespace sl {

  /// constants and numeric functions of bounded scalar type
  template<class T>
  class scalar_math< bounded_scalar<T> >: public scalar_math_base< bounded_scalar<T> >  {
    typedef T value_t;
    typedef bounded_scalar<T> bounded_scalar_t;
    typedef SL_FLOATTYPENAME(T) float_t;
    typedef bounded_scalar< float_t > float_bounded_scalar_t;
  public:
    static inline bounded_scalar_t zero                   () { return bounded_scalar_t(); }
    static inline bounded_scalar_t one                    () { return bounded_scalar_t(scalar_math<value_t>::one()); }
    static inline bounded_scalar_t two                    () { return bounded_scalar_t(scalar_math<value_t>::two()); }
    static inline float_bounded_scalar_t E                () { return float_bounded_scalar_t(sl::E(float_t())); }
    static inline float_bounded_scalar_t Log2_E           () { return float_bounded_scalar_t(sl::Log2_E(float_t())); }
    static inline float_bounded_scalar_t Log10_E          () { return float_bounded_scalar_t(sl::Log10_E(float_t())); }
    static inline float_bounded_scalar_t Ln_2             () { return float_bounded_scalar_t(sl::Ln_2(float_t())); }
    static inline float_bounded_scalar_t Ln_10            () { return float_bounded_scalar_t(sl::Ln_10(float_t())); }
    static inline float_bounded_scalar_t Pi               () { return float_bounded_scalar_t(sl::Pi(float_t())); }
    static inline float_bounded_scalar_t Pi_Over_2        () { return float_bounded_scalar_t(sl::Pi_Over_2(float_t())); }
    static inline float_bounded_scalar_t Pi_Over_4        () { return float_bounded_scalar_t(sl::Pi_Over_4(float_t())); }
    static inline float_bounded_scalar_t One_0ver_Pi      () { return float_bounded_scalar_t(sl::One_Over_Pi(float_t())); }
    static inline float_bounded_scalar_t Two_0ver_Pi      () { return float_bounded_scalar_t(sl::Two_Over_Pi(float_t())); }
    static inline float_bounded_scalar_t Two_Over_Sqrt_Pi () { return float_bounded_scalar_t(sl::Two_Over_Sqrt_Pi(float_t())); }
    static inline float_bounded_scalar_t Sqrt_2           () { return float_bounded_scalar_t(sl::Sqrt_2(float_t())); }
    static inline float_bounded_scalar_t One_Over_Sqrt_2  () { return float_bounded_scalar_t(sl::One_Over_Sqrt_2(float_t())); }
  };

} // namespace sl


// basic functions of bounded scalars
namespace sl {

  template <class T> inline bounded_scalar<T> sign(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(sign(x.value()), sign(x.bounds()));
  }

  template <class T> inline bounded_scalar<T> non_zero_sign(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(non_zero_sign(x.value()), non_zero_sign(x.bounds()));
  }
    
  template <class T> inline bounded_scalar<T> round_down(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(round_down(x.value()), round_down(x.bounds()));
  }

  template <class T> inline bounded_scalar<T> round_up(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(round_up(x.value()), round_up(x.bounds()));
  }
  
  template <class T> 
  static inline bounded_scalar<T> pred(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(pred(x.value()), pred(x.bounds()), true);
  }
  
  template <class T> 
  static inline bounded_scalar<T> succ(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(succ(x.value()), succ(x.bounds()), true);
  }

  template <class T> inline bounded_scalar<T> abs(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(abs(x.value()), abs(x.bounds()));
  }
  
  template <class T> inline bounded_scalar<T> integer_round_towards_zero(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(integer_round_towards_zero(x.value()), 
			     integer_round_towards_zero(x.bounds()));
  }
  
  template <class T1, class T2> 
  inline bounded_scalar<SL_PROMOTENAME(T1,T2)> mod(const bounded_scalar<T1>& x, const bounded_scalar<T2>& y) {
    return bounded_scalar<SL_PROMOTENAME(T1,T2)>(mod(x.value(),y.value()), mod(x.bounds(),y.bounds()), true);
  }
  
  template <class T> inline bounded_scalar<T> sqr(const bounded_scalar<T>& x) {
    return bounded_scalar<T>(sqr(x.value()), sqr(x.bounds()), true);
  }
}

// standard basic functions of bounded scalars
namespace std {

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> exp(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(exp(x.value()), exp(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> log(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(log(x.value()), log(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> log10(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(log10(x.value()), log10(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<T> pow(const sl::bounded_scalar<T>& x, int n) { 
    return sl::bounded_scalar<T>(pow(x.value(),n), pow(x.bounds(),n), true);
  }    
	
  template <class T1, class T2> inline sl::bounded_scalar<SL_PROMOTENAME(SL_FLOATTYPENAME(T1),SL_FLOATTYPENAME(T2))> pow(const sl::bounded_scalar<T1>& x, 
															      const sl::bounded_scalar<T2>& y) { 
    typedef SL_PROMOTENAME(SL_FLOATTYPENAME(T1),SL_FLOATTYPENAME(T2)) T;
    return sl::bounded_scalar<T>(pow(x.value(),y.value()), pow(x.bounds(),y.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<T> ceil(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<T>(ceil(x.value()), ceil(x.bounds()));
  }

  template <class T> inline sl::bounded_scalar<T> floor(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<T>(floor(x.value()), floor(x.bounds()));
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> sqrt(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(sqrt(x.value()), sqrt(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> atan(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(atan(x.value()), atan(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> acos(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(acos(x.value()), acos(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> asin(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(asin(x.value()), asin(x.bounds()), true);
  }

  template <class T1, class T2> 
  inline sl::bounded_scalar<SL_PROMOTENAME(SL_FLOATTYPENAME(T1),
					   SL_FLOATTYPENAME(T2))> atan2(const sl::bounded_scalar<T1>& y, 
									     const sl::bounded_scalar<T2>& x) {
    return 
      sl::bounded_scalar<SL_PROMOTENAME(SL_FLOATTYPENAME(T1),SL_FLOATTYPENAME(T2))>(atan2(y.value(),x.value()), atan2(y.bounds(),x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> sin(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(sin(x.value()), sin(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> cos(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(cos(x.value()), cos(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> tan(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(tan(x.value()), tan(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> cosh(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(cosh(x.value()), cosh(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> sinh(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(sinh(x.value()), sinh(x.bounds()), true);
  }

  template <class T> inline sl::bounded_scalar<SL_FLOATTYPENAME(T)> tanh(const sl::bounded_scalar<T>& x) {
    return sl::bounded_scalar<SL_FLOATTYPENAME(T)>(tanh(x.value()), tanh(x.bounds()), true);
  }
} // namespace std


// Useful typedefs
namespace sl {
  /// A bounded scalar built around a single precision floating point number
  typedef bounded_scalar<float> bounded_scalarf;
  /// A bounded scalar built around a double precision floating point number
  typedef bounded_scalar<double> bounded_scalard;
} // namespace sl

#endif






