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
#ifndef SL_INTERVAL_HPP
#define SL_INTERVAL_HPP

#include <sl/utility.hpp> // For min/max
#include <sl/operators.hpp>
#include <sl/assert.hpp>
#include <sl/cast.hpp>
#include <sl/math.hpp>
#include <sl/numeric_traits.hpp>
#include <sl/serializer.hpp>
#include <iostream>

namespace sl {
 
  namespace tags {
    /**
     * A struct for selecting rounded-out operations
     */
    struct rounded_out { };
    /**
     * A struct for selecting operations without rounding
     */
    struct not_rounded_out { };
  };

  /**
   * Intervals of scalar values.
   */ 
  template <class T>
  class interval {

  public: // Constants and types

    typedef interval<T>               self_t;
    typedef T                         value_t;
    typedef SL_FLOATTYPENAME(value_t) float_t;

  protected: // Constraints

    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<value_t>::is_specialized);

  protected: // Data

    value_t limits_[2];

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << limits_[0] << limits_[1];
    }
    
    void retrieve_from(input_serializer& s) {
      s >> limits_[0] >> limits_[1];
    }

  public: // Creation, Copy & Destruction

    /// Default constructor, initializing to [-0,+0]
    inline interval() {
      limits_[0] = -scalar_math<value_t>::zero();
      limits_[1] =  scalar_math<value_t>::zero();
    }

    /// Constructor setting limits to the parameters, without rounding
    inline interval(const value_t& infimum, 
		    const value_t& supremum) {
      limits_[0] = infimum;
      limits_[1] = supremum;
    }

    /// Constructor setting limits to the parameters, without rounding
    inline interval(const value_t& infimum, 
		    const value_t& supremum,
		    const tags::not_rounded_out) {
      limits_[0] = infimum;
      limits_[1] = supremum;
    }

    /// Constructor setting limits to the parameters, rounded out if needed.
    inline interval(const value_t& infimum, 
		    const value_t& supremum,
		    const tags::rounded_out) {
      limits_[0] = round_down(infimum);
      limits_[1] = round_up(supremum);
    }

    /// Constructor setting both limits to value, no rounding
    inline interval(const value_t& value) {
      limits_[0] = value;
      limits_[1] = value;
    }

    /// Constructor setting both limits to value, no rounding
    inline interval(const value_t& value,
		    const tags::not_rounded_out) {
      limits_[0] = value;
      limits_[1] = value;
    }

    /// Constructor setting both limits to value, rounded out if needed.
    inline interval(const value_t& value,
		    const tags::rounded_out) {
      limits_[0] = round_down(value);
      limits_[1] = round_up(value);
    }

    /**
     * Constructor setting the limits from a manifest constant.
     * The range of the interval is given by the number of digits
     * of precision in the specification. For example, 
     * interval("0.314") = [0.314-0.0005, 0.314+0.0005]
     */
    interval(const char* s) {
      SL_REQUIRE("String exists", s);

      value_t result = zero(value_t());
      int sign, ex = 0;
      // Eat whitespace
      int i = 0;
      while (s[i] == ' ' || s[i] == '\t' || s[i] == '\n') i++;
      // Get sign of mantissa
      switch (s[i]) { 
      case '-': { sign = -1;  i++; break; }
      case '+': { sign =  1; i++; break; }
      default:  { sign = 1; break; }
      }
      // get digits before decimal point
      int n;
      while (n=(s[i])-'0', i++, n>=0 && n<10) {
	result = value_t(10.0f)*result+(value_t)((float)n);
      }
      i--;
      // get digits after decimal point 
      if (s[i] == '.') {  
	i++;
	while (n=(s[i])-'0', i++, n>=0 && n<10) { 
	  result = value_t(10.0f)*result+(value_t)((float)n); --ex; 
	}
	i--;
      }
      // Get range
      value_t result_inf = value_t(10.0f) * result - value_t(5.0f);
      value_t result_sup = value_t(10.0f) * result + value_t(5.0f);
      --ex;

      // Get exponent
      if (s[i] == 'e' || s[i] == 'E') { 
	i++;

	int ex2 = 0, exsign = 1;
	// Get sign of exponent
	switch (s[i]) { 
	case '-': { exsign = -1;  i++; break; }
	case '+': { exsign =  1;  i++; break; }
	default:  { exsign = 1; break; }
	}
	// get digits before decimal point
	while (n=(s[i])-'0', i++, n>=0 && n<10) {
	  ex2 = 10*ex2+n;
	}
	i--; 
	if (exsign < 0) {
	  ex2 = -ex2;
	}
	ex += ex2;
      }
      // Eat whitespace
      while (s[i] == ' ' || s[i] == '\t' || s[i] == '\n') i++;
      SL_CHECK("No garbage in string", s[i] == '\0');

      // exponent adjustment
      while (ex-- > 0) {
	result_inf *= value_t(10.0f);
	result_sup *= value_t(10.0f);
      }

      while (++ex < 0) {
	result_inf /= value_t(10.0f);
	result_sup /= value_t(10.0f);
      }

      if (sign<0) {
	result_inf = -result_sup;
	result_sup = -result_inf;
      }

      limits_[0] = result_inf;
      limits_[1] = result_sup;
      SL_ENSURE("Not empty", !is_empty());
    }                                                                                                                         

  public: // Copy and Conversion

    /// Set this to other.
    inline interval(const self_t& other) {
      limits_[0] = other[0];
      limits_[1] = other[1];
    }

    /// Set this to other, converting limits from another numeric type.
    template <class OTHER_T>
    explicit interval(const interval<OTHER_T>& other) {
      limits_[0] = numeric_cast<T>(other[0]);
      limits_[1] = numeric_cast<T>(other[1]);
      // Check: round-out when reducing precision??
    }

  public: // Access

    /// Access to end-points (0 = inf, 1 = sup)
    inline const value_t& operator[](int i) const {
      SL_REQUIRE("Good index", 0 <= i && i<= 1);
      return limits_[i];
    }

    /// The lower bound
    inline const value_t& lower_bound() const {
      return limits_[0];
    }

    /// The upper bound
    inline const value_t& upper_bound() const {
      return limits_[1];
    }
    
  public: // Comparison

    /// Is this exactly equal to other? 
    inline bool operator== (const self_t& other) const {
      return limits_[0] == other[0] && limits_[1] == other[1];
    }

    SL_OP_EQUALITY_COMPARABLE1(self_t);

    /// Is this certainly less than other?
    inline bool certainly_less (const self_t& other) const {
      return limits_[1] < other[0];
    }

    /// Is this certainly greater than other?
     inline bool certainly_greater (const self_t& other) const {
      return other[1] < limits_[0];
    }

    /// Is this certainly less than or equal to other?
     inline bool certainly_less_or_equal (const self_t& other) const {
      return limits_[1] <= other[0];
    }

    /// Is this certainly greater than or equal to other?
    inline bool certainly_greater_or_equal (const self_t& other) const {
      return other[1] <= limits_[0];
    }

    /// Is this possibly equal to other?
    inline bool possibly_equal (const self_t& other) const {
      return !(certainly_less(other) || certainly_greater(other));
    }

    /// Is this possibly less than other?
    inline bool possibly_less (const self_t& other) const {
      return limits_[0] < other[1];
    }

    /// Is this possibly greater than other?
    inline bool possibly_greater (const self_t& other) const {
      return other[0] < limits_[1];
    }

    /// Is this possibly less than or equal to other?
    inline bool possibly_less_or_equal (const self_t& other) const {
      return limits_[0] <= other[1];
    }

    /// Is this possibly greater than or equal to other?
    inline bool possibly_greater_or_equal (const self_t& other) const {
      return other[0] <= limits_[1];
    }

  public: // Set operations

    /// Merge this with value x
    inline void merge(const value_t& x) {
      limits_[0]=min(limits_[0],x);
      limits_[1]=max(limits_[1],x);
    }
  
    /// Merge this with interval x
    inline void merge(const self_t& other) {
      limits_[0]=min(limits_[0],other[0]);
      limits_[1]=max(limits_[1],other[1]);
    }

    /// Intersect this with interval x
    inline void intersect(const self_t& other) {
      limits_[0]=max(limits_[0],other[0]);
      limits_[1]=min(limits_[1],other[1]);
    }                                                                                               

  public: // Set queries

    /// Is this empty?
    inline bool is_empty() const {
      return limits_[1] < limits_[0];
    }

    /// Is this overlapping other?
    inline bool is_overlapping(const self_t& other) const {
      return !((limits_[1]<=other[0]) || (limits_[0]>=other[1]));
    }

    /// Is this disjoint from other?
    inline bool is_disjoint(const self_t& other) const {
      return !is_overlapping(other);
    }

    /// Does this contain value?
    inline bool contains(const value_t& value) const {
      return (limits_[0]<=value) && (value <= limits_[1]);
    }

    /// Does this fully contain other?
    inline bool contains(const self_t& other) const {
      return contains(other[0]) && contains(other[1]);
    }

  public: // Queries

    /// The center of the interval
    inline value_t center() const {
      SL_REQUIRE("Not empty", !is_empty());
      float_t result = median((limits_[0] + limits_[1])/two(value_t()),
			      limits_[0],
			      limits_[1]);
      SL_ENSURE("Containment", contains(result));
      return result;
    }

    /// The center of the interval, in floating point
    inline float_t floatcenter() const {
      SL_REQUIRE("Not empty", !is_empty());
      float_t result = median(float_t(limits_[0] + limits_[1])/two(float_t()),
			      float_t(limits_[0]),
			      float_t(limits_[1]));
      SL_ENSURE("Containment", contains(result));
      return result;
    }

    /// The width of the interval
    inline value_t width() const {
      value_t result;
      if (is_empty()) {
	result = zero(value_t());
      } else {
	result = round_up(limits_[1] - limits_[0]);
      }
      SL_ENSURE("Non negative", result >= 0);
      return result;
    }
    
  public: // Field Operations

    /// Increment this by other
    inline self_t& operator += (const interval<T>& other) {
      SL_REQUIRE("This not empty", !is_empty());
      SL_REQUIRE("Other not empty", !other.is_empty());

      if (is_zero(other[0])) {
	// No op
      } else if (is_zero(limits_[0])) {
	limits_[0] = other[0];
      } else {
	limits_[0] = round_down(limits_[0]+other[0]);
      }

      if (is_zero(other[1])) {
	// No op
      } else if (is_zero(limits_[1])) {
	limits_[1] = other[1];
      } else {
	limits_[1] = round_up(limits_[1]+other[1]);
      }

      return *this;
    }

    /// The negation of this
    inline self_t operator - () const {
      SL_REQUIRE("This not empty", !is_empty());
      
      self_t result(-limits_[1], -limits_[0], tags::not_rounded_out());
      return result;
    }

    /// Decrement this by other
    inline self_t& operator -= (const self_t& other) {
      SL_REQUIRE("This not empty", !is_empty());
      SL_REQUIRE("Other not empty", !other.is_empty());

      if (is_zero(other[1])) {
	// No op
      } else if (is_zero(limits_[0])) {
	limits_[0] = -other[1];
      } else {
	limits_[0] = round_down(limits_[0]-other[1]);
      }

      if (is_zero(other[0])) {
	// No op
      } else if (is_zero(limits_[1])) {
	limits_[1] = -other[0];
      } else {
	limits_[1] = round_up(limits_[1]-other[0]);
      }

      return *this;
    }
    
    /// Post-multiply this with other
    /*inline*/ self_t& operator *= (const self_t& other) {
      SL_REQUIRE("This not empty", !is_empty());
      SL_REQUIRE("Other not empty", !other.is_empty());

      value_t l,h;      

      if (is_non_negative(limits_[0])) {
	if (is_non_negative(other[0])) {
	  // + +
	  l = limits_[0] * other[0];
	  h = limits_[1] * other[1];
	} else if (is_non_positive(other[1])) {
	  // + -
	  l = limits_[1] * other[0];
	  h = limits_[0] * other[1];
	} else{
	  // + ?
	  l = limits_[1] * other[0];
	  h = limits_[1] * other[1];
	}
      } else if (is_non_negative(limits_[1])) {
	if (is_non_negative(other[0])) {
	  // - +
	  l = limits_[0] * other[1];
	  h = limits_[1] * other[0];
	} else if (is_non_positive(other[1])) {
	  // - -
	  l = limits_[1] * other[1];
	  h = limits_[0] * other[0];
	} else{
	  // - ?
	  l = limits_[0] * other[1];
	  h = limits_[0] * other[0];
	}
      } else {
	// ? ?
	l = min(limits_[0]*other[1], limits_[1]*other[0]);
	h = max(limits_[0]*other[0], limits_[1]*other[1]);
      }
      
      limits_[0] = round_down(l);
      limits_[1] = round_up(h);
		
      return *this;
    }

    /// Divide this by other
    self_t& operator /= (const self_t& other) {
      SL_REQUIRE("This not empty", !is_empty());
      SL_REQUIRE("Other not empty", !other.is_empty());
      SL_REQUIRE("Other does not contain 0", !other.contains(zero(value_t()))); 

      if (other[0] > basic_op_zero_test(value_t())) {
	*this *= self_t(value_t(1) / other[1],
			value_t(1) / other[0],
			tags::rounded_out());
      } else if (other[1] < -basic_op_zero_test(value_t())) {
	*this *= self_t(value_t(1) / other[1],
			value_t(1) / other[0],
			tags::rounded_out());
      } else {
	SL_CHECK_WARNING("Division by zero", true);
	if ((limits_[0] >= 0 && other[0] >= 0) ||
	    (limits_[1] < 0 && other[1] < 0)) {
	  // +/+ -/- => +
	  *this = self_t(zero(value_t()),
			 sl::upper_bound(value_t()),
			 sl::tags::not_rounded_out());
	} else if ((limits_[0] >= 0 && other[1] < 0) ||
		   (limits_[1] < 0 && other[0] >= 0)) {
	  // +/- -/+ -> +
	  *this = self_t(sl::lower_bound(value_t()),
			 zero(value_t()),
			 sl::tags::not_rounded_out());
	} else {
	  // ? -> ?
	  *this = self_t(sl::lower_bound(value_t()), 
			 sl::upper_bound(value_t()),
			 sl::tags::not_rounded_out());
	}
      }
      return *this;
    }

    SL_OP_FIELD(self_t); // Without type promotion

    /// The multiplication of this by other
    template <class T2>
    inline interval<SL_PROMOTENAME(value_t,T2)> operator*(const interval<T2>& y) const { 
      interval<SL_PROMOTENAME(value_t,T2)> x(*this); 
      interval<SL_PROMOTENAME(value_t,T2)> y2(y); 
      return x *= y2; 
    }

    /// The division of this by other
    template <class T2>
    inline interval<SL_PROMOTENAME(SL_FLOATTYPENAME(value_t),SL_FLOATTYPENAME(T2))> operator/(const interval<T2>& y) const { 
      interval<SL_PROMOTENAME(SL_FLOATTYPENAME(value_t),SL_FLOATTYPENAME(T2))> x(*this); 
      interval<SL_PROMOTENAME(SL_FLOATTYPENAME(value_t),SL_FLOATTYPENAME(T2))> y2(y); 
      return x /= y2; 
    }

    /// The sum of this and other
    template <class T2>
    inline interval<SL_PROMOTENAME(value_t,T2)> operator+(const interval<T2>& y) const { 
      interval<SL_PROMOTENAME(value_t,T2)> x(*this); 
      interval<SL_PROMOTENAME(value_t,T2)> y2(y); 
      return x += y2; 
    }
    
    /// The subtraction of this and other
    template <class T2>
    inline interval<SL_PROMOTENAME(value_t,T2)> operator-(const interval<T2>& y) const { 
      interval<SL_PROMOTENAME(value_t,T2)> x(*this); 
      interval<SL_PROMOTENAME(value_t,T2)> y2(y); 
      return x -= y2; 
    }
    
  };
} // namespace sl

// I/O

template <class T>
std::ostream& operator << (std::ostream & os, const sl::interval<T>& x) {
  os << "[ " << x[0] << " .. " << x[1] << " ]";
  return os;
}

namespace std {

  /// numeric limits for a bounded scalar
  template<class T> class numeric_limits< sl::interval<T> > {
  public:
    static const bool is_specialized = true;

    static sl::interval<T> min() throw() { return sl::interval<T>(numeric_limits<T>::min()); }
    static sl::interval<T> max() throw() { return sl::interval<T>(numeric_limits<T>::max()); }

    static const int digits = numeric_limits<T>::digits;
    static const int digits10 = numeric_limits<T>::digits10;
    static const bool is_signed = numeric_limits<T>::is_signed;
    static const bool is_integer = numeric_limits<T>::is_integer;
    static const bool is_exact = numeric_limits<T>::is_exact;
    static const int radix =  numeric_limits<T>::radix;
    static sl::interval<T> epsilon() throw() { return sl::interval<T>(numeric_limits<T>::epsilon()); }
    static sl::interval<T> round_error() throw() { return sl::interval<T>(numeric_limits<T>::round_error()); }

    static const int min_exponent = numeric_limits<T>::min_exponent;
    static const int min_exponent10 = numeric_limits<T>::min_exponent10;
    static const int max_exponent = numeric_limits<T>::max_exponent;
    static const int max_exponent10 = numeric_limits<T>::max_exponent10;

    static const bool has_infinity = numeric_limits<T>::has_infinity;
    static const bool has_quiet_NaN =  numeric_limits<T>::has_quiet_NaN;
    static const bool has_signaling_NaN = numeric_limits<T>::has_signaling_NaN;
    static const float_denorm_style has_denorm = /*(float_denorm_style)*/numeric_limits<T>::has_denorm;
    static const bool has_denorm_loss = numeric_limits<T>::has_denorm_loss;

    static sl::interval<T> infinity() throw() { return sl::interval<T>(numeric_limits<T>::infinity()); }
    static sl::interval<T> quiet_NaN() throw() { return sl::interval<T>(numeric_limits<T>::quiet_NaN()); }
    static sl::interval<T> signaling_NaN() throw() { return sl::interval<T>(numeric_limits<T>::signaling_NaN()); }
    static sl::interval<T> denorm_min() throw() { return sl::interval<T>(numeric_limits<T>::denorm_min()); }
 
    static const bool is_iec559 = numeric_limits<T>::is_iec559;
    static const bool is_bounded = numeric_limits<T>::is_bounded;
    static const bool is_modulo = numeric_limits<T>::is_modulo;
    
    static const bool traps = numeric_limits<T>::traps;
    static const bool tinyness_before = numeric_limits<T>::tinyness_before;
    static const float_round_style round_style = /*(float_round_style)*/numeric_limits<T>::round_style;
  };

};

namespace sl {

  /// Constants and generic functions of interval type
  template<class T>
  class scalar_math< interval<T> >: public scalar_math_base< interval<T> > {
    typedef T value_t;
    typedef interval<T> interval_t;
    typedef SL_FLOATTYPENAME(T) float_t;
    typedef interval< float_t > float_interval_t;
  public:
    static inline interval_t zero                  () { return interval_t(scalar_math<value_t>::zero(),tags::not_rounded_out()); }
    static inline interval_t one                   () { return interval_t(scalar_math<value_t>::one()); }
    static inline interval_t two                   () { return interval_t(scalar_math<value_t>::two()); }
    static inline float_interval_t E                () { return float_interval_t(sl::E(float_t())); }
    static inline float_interval_t Log2_E           () { return float_interval_t(sl::Log2_E(float_t())); }
    static inline float_interval_t Log10_E          () { return float_interval_t(sl::Log10_E(float_t())); }
    static inline float_interval_t Ln_2             () { return float_interval_t(sl::Ln_2(float_t())); }
    static inline float_interval_t Ln_10            () { return float_interval_t(sl::Ln_10(float_t())); }
    static inline float_interval_t Pi               () { return float_interval_t(sl::Pi(float_t())); }
    static inline float_interval_t Pi_Over_2        () { return float_interval_t(sl::Pi_Over_2(float_t())); }
    static inline float_interval_t Pi_Over_4        () { return float_interval_t(sl::Pi_Over_4(float_t())); }
    static inline float_interval_t One_0ver_Pi      () { return float_interval_t(sl::One_Over_Pi(float_t())); }
    static inline float_interval_t Two_0ver_Pi      () { return float_interval_t(sl::Two_Over_Pi(float_t())); }
    static inline float_interval_t Two_Over_Sqrt_Pi () { return float_interval_t(sl::Two_Over_Sqrt_Pi(float_t())); }
    static inline float_interval_t Sqrt_2           () { return float_interval_t(sl::Sqrt_2(float_t())); }
    static inline float_interval_t One_Over_Sqrt_2  () { return float_interval_t(sl::One_Over_Sqrt_2(float_t())); }
  };

} // namespace sl

// Helper for defining postconditions
#define SL_ENSURE_INTERVAL_FN1(_fn_, _x_) \
      SL_ENSURE("Containment", \
                result.contains(_fn_(_x_[0])) && \
                result.contains(_fn_(_x_[1])) && \
		result.contains(_fn_(_x_.center())))

//---- basic functions of intervals
namespace sl {

  /// 1 if x > 0, -1 if x < 0, 0 otherwise 
  template<class T> static inline interval<T> sign(const interval<T>& x) {
    interval<T> result; 
    if (is_positive(x[0])) {
      result = interval<T>(scalar_math<T>::one(), scalar_math<T>::one());
    } else if (is_negative(x[1])) {
      result = interval<T>(-scalar_math<T>::one(), -scalar_math<T>::one());
    } else if (is_zero(x)) {
      result = interval<T>(-scalar_math<T>::zero(), scalar_math<T>::zero());
    } else {
      result = interval<T>(-scalar_math<T>::one(), scalar_math<T>::one());
    }
    return result;
  }

  /// 1 if x >= 0, -1 if otherwise 
  template<class T> static inline interval<T> non_zero_sign(const interval<T>& x) {
    interval<T> result; 
    if (is_positive(x[0])) {
      result = interval<T>(scalar_math<T>::one(), scalar_math<T>::one());
    } else if (is_negative(x[1])) {
      result = interval<T>(-scalar_math<T>::one(), -scalar_math<T>::one());
    } else if (is_zero(x)) {
      result = interval<T>(scalar_math<T>::one(), scalar_math<T>::one());
    } else {
      result = interval<T>(-scalar_math<T>::one(), scalar_math<T>::one());
    } 
    return result;
  }

  /// The previous representable number before x
  template <class T> 
  inline interval<T> pred(const interval<T>& x) {
    return interval<T>(pred(x[0]), pred(x[1]));
  }
  
  /// The next representable number after x
  template <class T> 
  inline interval<T> succ(const interval<T>& x) {
    return interval<T>(succ(x[0]),succ(x[1]));
  }

  /// Directed rounding towards minus infinity of infimum
  template <class T> inline interval<T> round_down(const interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    interval<T> result = interval<T>(round_down(x[0]), x[1]);
    return result;
  }

  /// Directed rounding towards infinity of supremum
  template <class T> inline interval<T> round_up(const interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    interval<T> result = interval<T>(x[0], round_up(x[1]));
    return result;
  }

  /// The remainder of x divided by two Pi
  template <class T> interval<SL_FLOATTYPENAME(T)> mod_two_pi(const interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    
    const float_t TwoPi = two(float_t()) * Pi(float_t());
    float_t q;
    float_t x_inf, x_sup;
    
    x_inf = x[0];
    x_sup = x[1];
  
    if (0 <= x_inf && x_sup < TwoPi) {
      // Nothing to do
    } else {
      q = x[0] / TwoPi;
      q = std::floor(pred(q));
      q *= TwoPi;
      q = succ(q);
      if (q > 0) {
	q = round_up(q); // error estimation: 1 unit in LSB
	q = succ(q);
      }
      x_inf = x[0] - q;
      x_inf = pred(x_inf);
      if (x_inf < 0) {
	x_inf += TwoPi;
      x_inf = max(zero(T()),pred(x_inf));
      }
      q = x[1] / TwoPi;
      q = std::floor(pred(q));
      q *= TwoPi;
      q = pred(q);
      if (q < 0) {
	q = round_down(q); // error estimation: 1 unit in LSB
	q = pred(q);
      }
      x_sup = x[1] - q;
      x_sup = succ(x_sup);
      if (x_sup > TwoPi) {
	x_sup -= TwoPi;
	x_sup = min(zero(float_t()),succ(x_sup));
      }
    }
    
    interval<float_t> result = interval<float_t>(min(x_inf,x_sup), 
							 max(x_inf,x_sup));
    
    SL_ENSURE("Not empty", !result.is_empty());
    SL_ENSURE("Inside [0..2Pi]", zero(float_t()) <= result[0] && result[1] <= TwoPi); 
    
    return result;
  }

  /// Absolute value of x
  template <class T> inline interval<T> abs(const interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    interval<T> result;
    
    if (x[0] > 0) {
      result = x;
    } else if (x[1] < 0) {
      result = interval<T>(-x[1], -x[0]);
    } else {
      result = interval<T>(zero(T()), max(abs(x[0]), abs(x[1])));
    }
    
    SL_ENSURE("Non negative", result[0] >= 0);
    SL_ENSURE_INTERVAL_FN1(abs,x);
    return result;
  }


  /// the integral value between zero and x which is closest to x
  template <class T> inline interval<T> integer_round_towards_zero(const interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    interval<T> result;
    if (x[0]>=0) {
      result = interval<T>(std::floor(x[0]), std::floor(x[1]));
    } else if (x[1]<0) {
      result = interval<T>(std::ceil(x[0]), std::ceil(x[1]));
    } else {
      result = interval<T>(std::ceil(x[0]), std::floor(x[1]));
    }
    SL_ENSURE_INTERVAL_FN1(integer_round_towards_zero,x);
    return result;
  }

  /// the remainder of x/y
  template <class T1, class T2>
  inline interval<SL_PROMOTENAME(T1,T2)> mod(const interval<T1>& x, 
					     const interval<T2>& y) {
    SL_REQUIRE("Not empty", !x.is_empty());
    interval<SL_PROMOTENAME(T1,T2)> result = x - integer_round_towards_zero(x / y) * y;
    SL_ENSURE("Containment", result.contains(mod(x[0],y[0])) && result.contains(mod(x[1],y[0])) && result.contains(mod(x.center(),y[0])));
    return result;
  }

  /// x^2
  template <class T> inline interval<T> sqr(const interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    interval<T> abs_x = abs(x);
    interval<T> result = abs_x * abs_x;
    SL_ENSURE("Non negative", result[0] >= 0);
    SL_ENSURE_INTERVAL_FN1(sqr,x);
    return result;
  }
  
} // namespace sl

namespace std {

  /// e^x
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> exp(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    sl::interval<SL_FLOATTYPENAME(T)> result = sl::interval<SL_FLOATTYPENAME(T)>(sl::round_down(std::exp(x[0])),
										 sl::round_up(std::exp(x[1])));

    SL_ENSURE("Non negative", result[0] >= 0);
    SL_ENSURE_INTERVAL_FN1(exp,x);
    return result;
  }

  /// log_e(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> log(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    SL_ENSURE("Positive argument", x[0] > 0);
    sl::interval<SL_FLOATTYPENAME(T)> result = sl::interval<SL_FLOATTYPENAME(T)>(std::log(x[0]), std::log(x[1]), sl::tags::rounded_out()); 
    SL_ENSURE_INTERVAL_FN1(log,x);
    return result;
  }

  /// log_10(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> log10(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    SL_ENSURE("Positive argument", x[0] > 0);
    sl::interval<SL_FLOATTYPENAME(T)> result = sl::interval<SL_FLOATTYPENAME(T)>(std::log10(x[0]), std::log10(x[1])); 
    SL_ENSURE_INTERVAL_FN1(log10,x);
    return result;
  }

  /// x^n
  template <class T> sl::interval<T> pow(const sl::interval<T>& x, int n) { 
    SL_REQUIRE("Not empty", !x.is_empty());

    sl::interval<T> result;

    int absn = (n < 0) ? (-n) : n;

    if (absn == 0) {
      result = sl::one(sl::interval<T>());
    } else {
      sl::interval<T> z = x;
      result = x;
      absn--;
      do {
	int halfn = absn / 2;
	if (halfn + halfn != absn) {
	  result *= z;
	}
	absn = halfn;
	if (absn > 0) {
	  z = sl::sqr(z);
	}
      } while (absn > 0);

      if (n < 0) {
	result = sl::reciprocal(result);
      }
    }
    
    SL_ENSURE("Containment", result.contains(std::pow(x[0],n)) && result.contains(std::pow(x[1],n)) && result.contains(std::pow(x.center(),n)));
    return result;
  }    
	
  /// x^y
  template <class T, class U> sl::interval<SL_PROMOTENAME(SL_FLOATTYPENAME(T),SL_FLOATTYPENAME(U))> pow(const sl::interval<T>& x, const sl::interval<U>& y) { 
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_PROMOTENAME(SL_FLOATTYPENAME(T),SL_FLOATTYPENAME(U)) float_t;
    sl::interval<float_t> result = std::exp(y*std::log(x));
    //--- Should check numeric accuracy SL_ENSURE("Containment", result.contains(std::pow(x[0],y[0])) && result.contains(std::pow(x[1],y[0])) && result.contains(std::pow(x.center(),y[0])));
    return result;
  }

  /// smallest integral value not less than x 
  template <class T> inline sl::interval<T> ceil(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    sl::interval<T> result = sl::interval<T>(std::ceil(x[0]), std::ceil(x[1]));
    SL_ENSURE_INTERVAL_FN1(ceil,x);
    return result;
  }
  /// largest integral value not greater than x      
  template <class T> inline sl::interval<T> floor(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    sl::interval<T> result = sl::interval<T>(std::floor(x[0]), std::floor(x[1]));
    SL_ENSURE_INTERVAL_FN1(floor,x);
    return result;
  }

  /// sqrt(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> sqrt(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    SL_REQUIRE("Non negative x", x[0] >= sl::zero(T())); 
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = sl::interval<float_t>(sl::round_down(std::sqrt(static_cast<float_t>(x[0]))), sl::round_up(std::sqrt(static_cast<float_t>(x[1])))); 
    SL_ENSURE_INTERVAL_FN1(sqrt,x);
    return result;
  }

  /// atan(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> atan(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = sl::interval<float_t>(sl::max(sl::round_down(std::atan(static_cast<float_t>(x[0]))),-sl::Pi_Over_2(float_t())), 
							 sl::min(sl::round_up(std::atan(static_cast<float_t>(x[1]))), sl::Pi_Over_2(float_t())));
    SL_ENSURE_INTERVAL_FN1(atan,x);
    return result;
  }

  /// atan2(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> acos(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    SL_REQUIRE("Good range", x[0] >= -sl::one(T()) && x[1] <= sl::one(T()));
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = sl::interval<float_t>(sl::max(sl::round_down(std::acos(static_cast<float_t>(x[1]))), sl::zero(float_t())), 
							 sl::min(sl::round_up(std::acos(static_cast<float_t>(x[0]))), sl::Pi(float_t())));
    SL_ENSURE_INTERVAL_FN1(acos,x);
    return result;
  }

  /// asin(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> asin(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    SL_REQUIRE("Good range", x[0] >= -sl::one(T()) && x[1] <= sl::one(T()));
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = sl::interval<float_t>(sl::max(sl::round_down(std::asin(static_cast<float_t>(x[0]))),-sl::Pi_Over_2(float_t())),
							 sl::min(sl::round_up(std::asin(static_cast<float_t>(x[1]))),   sl::Pi_Over_2(float_t())));
    SL_ENSURE_INTERVAL_FN1(asin,x);
    return result;
  }

  /// arc tangent of y/x, with the sign of both values used to determine the result
  template <class T1, class T2> sl::interval<SL_PROMOTENAME(SL_FLOATTYPENAME(T1),SL_FLOATTYPENAME(T2))> atan2(const sl::interval<T1>& y, 
														   const sl::interval<T2>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());

    typedef SL_FLOATTYPENAME(T1) float_t1;
    typedef SL_FLOATTYPENAME(T2) float_t2;
    typedef SL_PROMOTENAME(float_t1, float_t2) float_t;
    sl::interval<float_t> result;

    if (x.certainly_greater(sl::zero(T1()))) {
      result = sl::two(float_t())*std::atan(y/(x + std::sqrt(sl::sqr(x) + sl::sqr(y))));
    } else if (!y.contains(sl::zero(T2()))) {
      result = sl::two(float_t())*std::atan((std::sqrt(sl::sqr(x) + sl::sqr(y)) - x)/ y);
    } else {
      SL_FAIL("Not implemented");
    }
    ("Containment", result.contains(std::atan2(y[0],x[0])) && result.contains(std::atan2(y[1],x[0])) && result.contains(std::atan2(y.center(),x[0])));

    return result;
  }

  /// sin(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> sin(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    // Most of this code has been freely adapted from the BIAS package

    typedef SL_FLOATTYPENAME(T) float_t;
 
    sl::interval<float_t> result = sl::interval<float_t>(-sl::one(float_t()), sl::one(float_t()));

    const float_t Pi = sl::Pi(float_t());
    const float_t Two_Pi = sl::two(float_t()) * Pi;

    if (x.width() >= Two_Pi) {
      // -1 .. 1
    } else {
      // Normalize the interval to 0 -- 2 Pi

      sl::interval<float_t> xx = sl::mod_two_pi(x);
      float_t x_inf = xx[0];
      float_t x_sup = xx[1];

      int q_inf = static_cast<int>(x_inf / sl::Pi_Over_2(float_t()));
      int q_sup = static_cast<int>(x_sup / sl::Pi_Over_2(float_t()));

      if ((q_inf == q_sup) && (x_sup > x_inf + Pi)) { 
	// -1 ... 1
      } else {

	float_t y_inf = sl::zero(float_t());
	float_t y_sup = sl::zero(float_t());

	switch((q_sup << 2) + q_inf)
	  {
	  case 0:
	  case 3:
	  case 15:
	    y_inf = sl::round_down(std::sin(x_inf));
	    y_sup = sl::round_up(std::sin(x_sup));
	    break;
	  case 1:
	  case 14:
	    y_inf = -sl::one(float_t());
	    x_inf = sl::round_up(std::sin(x_inf));
	    x_sup = sl::round_up(std::sin(x_sup));
	    y_sup = sl::max(x_inf, x_sup);
	    break;
	  case 2:
	    y_inf = -sl::one(float_t());
	    y_sup = sl::round_up(std::sin(x_sup));
	    break;
	  case 4:
	  case 11:
	    y_sup = sl::one(float_t());
	    x_inf = sl::round_down(std::sin(x_inf));
	    x_sup = sl::round_down(std::sin(x_sup));
	    y_inf = sl::min(x_inf, x_sup);
	    break;
	  case 5:
	  case 9:
	  case 10:
	    y_inf = sl::round_down(std::sin(x_sup));
	    y_sup = sl::round_up(std::sin(x_inf));
	    break;
	  case 6:
	  case 12:
	    // -1 .. 1
	    break;
	  case 7:
	    y_sup = sl::one(float_t());
	    y_inf = sl::round_down(std::sin(x_inf));
	    break;
	  case 8:
	    y_sup = sl::one(float_t());
	    y_inf = sl::round_down(std::sin(x_sup));
	    break;
	  case 13:
	    y_inf = -sl::one(float_t());
	    y_sup = sl::round_up(std::sin(x_inf));
	    break;
	  }
	if (y_inf < -sl::one(float_t())) y_inf = -sl::one(float_t());	
	if (y_sup > sl::one(float_t())) y_sup = sl::one(float_t());	
	result = sl::interval<T>(y_inf, y_sup);
      }
    }

    SL_ENSURE_INTERVAL_FN1(sin,x);
    SL_ENSURE("Good range", 
	      result[0] >= -sl::one(float_t()) && 
	      result[1] <= sl::one(float_t()));

    return result;
  }

  /// cos(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> cos(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = std::sin(x + sl::interval<float_t>(sl::Pi_Over_2(float_t())));

    // -- Should check numeric accuracy SL_ENSURE_INTERVAL_FN1(cos,x);
    return result;
  }
 
  /// tan(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> tan(const sl::interval<T>& x) {   
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result(sl::lower_bound(float_t()),
				 sl::upper_bound(float_t()));

    const float_t Pi    = sl::Pi(float_t());
    const float_t Two_Pi = sl::two(float_t()) * Pi;
 
    if (x.width() >= Two_Pi) {
      // Infinity 
    } else {
      sl::interval<float_t> xx = sl::mod_two_pi(x);
      float_t x_inf = xx[0];
      float_t x_sup = xx[1];
      if (x_sup >= x_inf + Two_Pi) { // security
	// Infinity
      } else {
	int q_inf = static_cast<int>(x_inf / sl::Pi_Over_2(float_t()));
	int q_sup = static_cast<int>(x_sup / sl::Pi_Over_2(float_t()));

	if ((q_inf == q_sup) && (x_sup > x_inf + Pi)) {
	  // Infinity 
	} else {
	  switch ((q_sup << 2) + q_inf) {
	  case 0:
	  case 3:
	  case 5:
	  case 9:
	  case 10:
	  case 15:
	    result = sl::interval<float_t>(sl::round_down(std::tan(x_inf)), sl::round_up(std::tan(x_sup)));
	    break;
	  default:
	    // infinity
	    break;
	  }
	}
      }
    }

    SL_ENSURE_INTERVAL_FN1(tan,x);
    return result;
  }

  /// cosh(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> cosh(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = (std::exp(x) + std::exp(-x)) / sl::two(sl::interval<float_t>());
    // -- Check accuracy! SL_ENSURE_INTERVAL_FN1(cosh,x);
    return result;
  }

  /// sinh(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> sinh(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = (std::exp(x) - std::exp(-x)) / sl::two(sl::interval<float_t>());
    // -- Check accuracy! SL_ENSURE_INTERVAL_FN1(sinh,x);
    return result;
  }

  /// std::tanh(x)
  template <class T> sl::interval<SL_FLOATTYPENAME(T)> tanh(const sl::interval<T>& x) {
    SL_REQUIRE("Not empty", !x.is_empty());
    typedef SL_FLOATTYPENAME(T) float_t;
    sl::interval<float_t> result = std::sinh(x)/std::cosh(x);
    // -- Check accuracy! SL_ENSURE_INTERVAL_FN1(tanh,x);
    return result;
  }

} // namespace std


namespace sl {
  /// An interval of double precision floating points
  typedef interval<double> intervald;
  /// An interval of single precision floating points
  typedef interval<float> intervalf;
} // namespace sl

#endif




