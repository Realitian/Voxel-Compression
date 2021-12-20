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
#ifndef SL_MATH_HPP
#define SL_MATH_HPP

#include <sl/assert.hpp>

#include <cmath>

#include <limits>

#include <sl/numeric_traits.hpp> // for SL_FLOATTYPE
#include <sl/type_traits.hpp> // for is_same ...

namespace sl {

  /// Base class for scalar math generic features
  template <class T>
  class scalar_math_base {

    SL_COMPILE_TIME_CHECK("Numeric value", std::numeric_limits<T>::is_specialized);

  public:

    typedef SL_FLOATTYPENAME(T) float_t;

  public: // Useful constants

    /// The constant 0
    static inline T zero() { return static_cast<T>(0); }
    /// The constant 1
    static inline T one()  { return static_cast<T>(1); }
    /// The constant 2
    static inline T two()  { return one()+one(); }

    /// The constant e
    static inline float_t E() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant log_2(e)
    static inline float_t Log2_E() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant log_10(e)
    static inline float_t Log10_E() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant ln(2)
    static inline float_t Ln_2() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant ln(10)
    static inline float_t Ln_10() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant Pi = 3.14...
    static inline float_t Pi() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant Pi/2
    static inline float_t Pi_Over_2() { return Pi() / two(); }
    /// The constant Pi/4
    static inline float_t Pi_Over_4() { return Pi() / (two()+two()); }
    /// The constant 1/Pi
    static inline float_t One_Over_Pi() { return one() / Pi(); }
    /// The constant 2/Pi
    static inline float_t Two_Over_Pi() { return two() / Pi(); }
    /// The constant 2/sqrt(Pi)
    static inline float_t Two_Over_Sqrt_Pi() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant sqrt(2)
    static inline float_t Sqrt_2() { SL_FAIL("Not implemented"); return float_t(zero()); }
    /// The constant 1/sqrt(2)
    static inline float_t One_Over_Sqrt_2() { return one() / Sqrt_2(); }

  public: // 

    /// max { x in T, x != oo }
    static inline T finite_upper_bound() { 
      return std::numeric_limits<T>::max();
    }

    /// min { x in T, x != oo }
    static inline T finite_lower_bound() { 
      if (std::numeric_limits<T>::is_signed) {
	return - finite_upper_bound();
      } else {
	return zero();
      }
    }

    /// max { x in T }
    static inline T upper_bound() { 
      if (std::numeric_limits<T>::has_infinity) {
	return std::numeric_limits<T>::infinity();
      } else {
	return finite_upper_bound();
      }
    }

    /// min { x in T }
    static inline T lower_bound() { 
      if (std::numeric_limits<T>::is_signed) {
	return - upper_bound();
      } else {
	return zero();
      }
    }

    /// min { x >= 0 : 1 + x > 1 }  
    static inline T epsilon() {   
      return std::numeric_limits<T>::epsilon(); 
    } 

    /// min { x in T : x > 0 }
    static inline T eta() {
      return std::numeric_limits<T>::min();
    }

    static inline T digits() {        
      return std::numeric_limits<T>::digits();
    }
    static inline int digits10() {        
      return std::numeric_limits<T>::digits10();
    }
    static inline int max_exponent() {        
      return std::numeric_limits<T>::max_exponent();
    }
    static inline int min_exponent() {        
      return std::numeric_limits<T>::min_exponent();
    }
    static inline int max_exponent10() {        
      return std::numeric_limits<T>::max_exponent10();
    }
    static inline int min_exponent10() {        
      return std::numeric_limits<T>::min_exponent10();
    }
    static inline int precision() {        
      return std::numeric_limits<T>::digits10();
    }
    static inline int radix() {        
      return std::numeric_limits<T>::radix;
    }

    static inline bool is_signed(T) {
      return std::numeric_limits<T>::is_signed;
    }

    static inline bool is_integer() {
      return std::numeric_limits<T>::is_integer;
    }

    static inline bool is_exact() {
      return std::numeric_limits<T>::is_integer;
    }

    static inline T round_error() {
      return std::numeric_limits<T>::round_error();
    }    

    static inline bool has_infinity() {
      return std::numeric_limits<T>::has_infinity;
    }
    
    static inline bool has_quiet_NaN() {
      return std::numeric_limits<T>::has_quiet_NaN;
    }

    static inline bool has_signaling_NaN() {
      return std::numeric_limits<T>::has_signaling_NaN;
    }

    static inline bool has_denorm() {
      return std::numeric_limits<T>::has_denorm;
    }

    static inline bool has_denorm_loss() {
      return std::numeric_limits<T>::has_denorm_loss;
    }

    static inline T infinity()  {
      return std::numeric_limits<T>::infinity();
    }
    
    static inline T quiet_NaN() {
      return std::numeric_limits<T>::quiet_NaN();
    }
    
    static inline T signaling_NaN()  {
      return std::numeric_limits<T>::signaling_NaN();
    }

    static inline T denorm_min()  {
      return std::numeric_limits<T>::denorm_min();
    }

    static inline bool is_iec559() {
      return std::numeric_limits<T>::is_iec559;
    }

    static inline bool is_bounded() {
      return std::numeric_limits<T>::is_bounded;
    }

    static inline bool is_modulo() {
      return std::numeric_limits<T>::is_modulo;
    }

    static inline bool traps() {
      return std::numeric_limits<T>::traps;
    }

    static inline bool tinyness_before() {
      return std::numeric_limits<T>::tinyness_before;
    }

    static inline std::float_round_style round_style() {
      return std::numeric_limits<T>::round_style;
    }

  public: // Basic operation precision

    static inline int basic_op_max_error_ulp() {
      return (std::numeric_limits<T>::is_integer || std::numeric_limits<T>::is_exact)? 0 : 1;
    }

    static inline T   basic_op_max_error() {
      return T(basic_op_max_error_ulp()) * epsilon();
    }

    static inline T   basic_op_zero_test() {
      return eta() / (one() - T(3) * basic_op_max_error());
    }

    static inline T   one_minus_basic_op_max_error() {
      return one() - basic_op_max_error();
    }

    static inline T   one_plus_basic_op_max_error() {
      return one() + basic_op_max_error();
    }

  };

}; // namespace sl

// --------------------------------------------------------------------------
// scalar_math<T>
// --------------------------------------------------------------------------

namespace sl {

  template <class T>
  class scalar_math: public scalar_math_base<T> {
  public:

    typedef SL_FLOATTYPENAME(T) float_t;

    static const bool t_is_float_t = is_same<T,float_t>::value;

  public:

    static inline float_t E                () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::E(); 
    }

    static inline float_t Log2_E           () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Log2_E(); 
    }

    static inline float_t Log10_E          () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Log10_E(); 
    }

    static inline float_t Ln_2             () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Ln_2(); 
    }

    static inline float_t Ln_10            () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Ln_10(); 
    }

    static inline float_t Pi               () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Pi(); 
    }

    static inline float_t Pi_Over_2        () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Pi_Over_2(); 
    }

    static inline float_t Pi_Over_4        () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Pi_Over_4(); 
    }

    static inline float_t One_Over_Pi      () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::One_Over_Pi(); 
    }

    static inline float_t Two_Over_Pi      () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Two_Over_Pi(); 
    }

    static inline float_t Two_Over_Sqrt_Pi () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::Two_Over_Sqrt_Pi(); 
    }

    static inline float_t Sqrt_2           () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t);
      return scalar_math<float_t>::Sqrt_2(); 
    }
    
    static inline float_t One_Over_Sqrt_2  () { 
      SL_COMPILE_TIME_CHECK("Default for non float types", !t_is_float_t); 
      return scalar_math<float_t>::One_Over_Sqrt_2(); 
    }
  };

}; // namespace sl

// --------------------------------------------------------------------------
// -- scalar_math<float>
// --------------------------------------------------------------------------

namespace sl {

  template<> class scalar_math<float>: public scalar_math_base<float> {
  public:
    static inline float zero             () { return  0.0f; }
    static inline float one              () { return  1.0f; }
    static inline float two              () { return  2.0f; }

    static inline float E                () { return  2.7182818284590452354f; } 
    static inline float Log2_E           () { return  1.4426950408889634074f; } 
    static inline float Log10_E          () { return  0.43429448190325182765f; }
    static inline float Ln_2             () { return  0.69314718055994530942f; }
    static inline float Ln_10            () { return  2.30258509299404568402f; }
    static inline float Pi               () { return  3.14159265358979323846f; }
    static inline float Pi_Over_2        () { return  1.57079632679489661923f; }
    static inline float Pi_Over_4        () { return  0.78539816339744830962f; }
    static inline float One_Over_Pi      () { return  0.31830988618379067154f; }
    static inline float Two_Over_Pi      () { return  0.63661977236758134308f; }
    static inline float Two_Over_Sqrt_Pi () { return  1.12837916709551257390f; }
    static inline float Sqrt_2           () { return  1.41421356237309504880f; }
    static inline float One_Over_Sqrt_2  () { return  0.70710678118654752440f; }
  };

}; // namespace sl

// --------------------------------------------------------------------------
// -- scalar_math<double>
// --------------------------------------------------------------------------

namespace sl {

  template<> class scalar_math<double>: public scalar_math_base<double> {
  public:
    static inline double zero             () { return  0.0; }
    static inline double one              () { return  1.0; }
    static inline double two              () { return  2.0; }

    static inline double E                () { return  2.7182818284590452354; } 
    static inline double Log2_E           () { return  1.4426950408889634074; } 
    static inline double Log10_E          () { return  0.43429448190325182765; }
    static inline double Ln_2             () { return  0.69314718055994530942; }
    static inline double Ln_10            () { return  2.30258509299404568402; }
    static inline double Pi               () { return  3.14159265358979323846; }
    static inline double Pi_Over_2        () { return  1.57079632679489661923; }
    static inline double Pi_Over_4        () { return  0.78539816339744830962; }
    static inline double One_Over_Pi      () { return  0.31830988618379067154; }
    static inline double Two_Over_Pi      () { return  0.63661977236758134308; }
    static inline double Two_Over_Sqrt_Pi () { return  1.12837916709551257390; }
    static inline double Sqrt_2           () { return  1.41421356237309504880; }
    static inline double One_Over_Sqrt_2  () { return  0.70710678118654752440; }
  };

}; // namespace sl

// --------------------------------------------------------------------------
// Directed Rounding 
// --------------------------------------------------------------------------
    
// Directed rounding is implemented in a reasonably
// transportable way by adjusting an interval's
// endpoints after an operation's result is computed
// using standard arithmetic. We do not use mode 
// switching in the FPU here because it introduces
// a shared resource management problem in parallel
// programs

namespace sl {

  /// Directed rounding towards minus infinity
  template <class T>
  inline T round_down(const T& x) {
    T result = - scalar_math<T>::basic_op_zero_test();

    if (x <= result) {  
      result = x * scalar_math<T>::one_plus_basic_op_max_error();
    } else if (x >= -result) {
      result = x * scalar_math<T>::one_minus_basic_op_max_error();
    } else if (x < scalar_math<T>::zero()) {
      // nothing to do
    } else {
      result = scalar_math<T>::zero();
    }
    
    SL_ENSURE("Rounded down: ", result <= x);
    SL_ENSURE("Consistent sign: ", (x <= scalar_math<T>::zero()) == (result <= scalar_math<T>::zero()) );
    
    return result;
  }
 
  /// Directed rounding towards plus infinity
  template <class T>
  inline T round_up(const T& x) {
    T result = scalar_math<T>::basic_op_zero_test();

    if (x >= result) {
      result = x * scalar_math<T>::one_plus_basic_op_max_error();
    } else if (x <= -result) {  
      result = x * scalar_math<T>::one_minus_basic_op_max_error();
    } else if (x > scalar_math<T>::zero()) {
      // Nothing to do
    } else {
      result = scalar_math<T>::zero();
    }
    
    SL_ENSURE("Rounded up: ", result >= x);
    SL_ENSURE("Consistent sign: ", (x >= scalar_math<T>::zero()) == (result >= scalar_math<T>::zero()) );

    return x; 
  }

  /// The previous representable number before x
  template <class T> 
  inline T pred(const T& x) {
    return (std::numeric_limits<T>::is_integer) ? (x-scalar_math<T>::one()) : round_down(x - std::numeric_limits<T>::min());
  }
  
  /// The next representable number after x
  template <class T> 
  inline T succ(const T& x) {
    return (std::numeric_limits<T>::is_integer) ? (x+scalar_math<T>::one()) : round_up(x + std::numeric_limits<T>::min());
  }

}; // namespace sl


// Manual instantiation for common basic
// types. This helps compilers
// such as gcc

namespace sl {

  template<>
  inline float round_down(const float& x) {
    const float zero_neg = - (std::numeric_limits<float>::min() / (1.0f - 3.0f * std::numeric_limits<float>::epsilon()));
    float result;
    if (x <= zero_neg) {  
      result = x * (1.0f + std::numeric_limits<float>::epsilon());
    } else if (x >= -zero_neg) {
      result = x * (1.0f - std::numeric_limits<float>::epsilon());
    } else if (x < 0.0f) {
      result = zero_neg;
    } else {
      result = 0.0f;
    }
    
    SL_ENSURE("Rounded down: ", result <= x);
    SL_ENSURE("Consistent sign: ", (x <= 0.0f) == (result <= 0.0f) );
    
    return result;
  }
 
  template <>
  inline float round_up(const float& x) {
    const float zero_pos = (std::numeric_limits<float>::min() / (1.0f - 3.0f * std::numeric_limits<float>::epsilon()));
    float result;
    if (x >= zero_pos) {
      result = x * (1.0f + std::numeric_limits<float>::epsilon());
    } else if (x <= -zero_pos) {  
      result = x * (1.0f - std::numeric_limits<float>::epsilon());
    } else if (x > 0.0f) {
      result = zero_pos;
    } else {
      result = 0.0f;
    }
    
    SL_ENSURE("Rounded up: ", result >= x);
    SL_ENSURE("Consistent sign: ", (x >= 0.0f) == (result >= 0.0f) );

    return result; 
  }

  template<>
  inline double round_down(const double& x) {
    const double zero_neg = - (std::numeric_limits<double>::min() / (1.0 - 3.0 * std::numeric_limits<double>::epsilon()));
    double result;
    if (x <= zero_neg) {  
      result = x * (1.0 + std::numeric_limits<double>::epsilon());
    } else if (x >= -zero_neg) {
      result = x * (1.0 - std::numeric_limits<double>::epsilon());
    } else if (x < 0.0) {
      result = zero_neg;
    } else {
      result = 0.0;
    }
    
    SL_ENSURE("Rounded down: ", result <= x);
    SL_ENSURE("Consistent sign: ", (x <= 0.0) == (result <= 0.0) );
    
    return result;
  }
 
  template <>
  inline double round_up(const double& x) {
    const double zero_pos = (std::numeric_limits<double>::min() / (1.0 - 3.0 * std::numeric_limits<double>::epsilon()));
    double result;

    if (x >= zero_pos) {
      result = x * (1.0 + std::numeric_limits<double>::epsilon());
    } else if (x <= -zero_pos) {  
      result = x * (1.0 - std::numeric_limits<double>::epsilon());
    } else if (x > 0.0) {
      result = zero_pos;
    } else {
      result = 0.0;
    }
    
    SL_ENSURE("Rounded up: ", result >= x);
    SL_ENSURE("Consistent sign: ", (x >= 0.0) == (result >= 0.0) );

    return result; 
  }


}; // namespace sl

// --------------------------------------------------------------------------
// -- Generic functions for inquiring numeric properties
// --------------------------------------------------------------------------

namespace sl {

  /// is x > 0?
  template<class _T_> static inline bool is_positive(const _T_& x) {
    return scalar_math<_T_>::zero() < x;
  }

  /// is x < 0?
  template<class _T_> static inline bool is_negative(const _T_& x) {
    return x < scalar_math<_T_>::zero();
  }

  /// is x == 0?
  template<class _T_> static inline bool is_zero(const _T_& x) {
    return x == scalar_math<_T_>::zero();
  }

  /// is x == 1?
  template<class _T_> static inline bool is_one(const _T_& x) {
    return x == scalar_math<_T_>::one();
  }

  /// is x >= 0?
  template<class _T_> static inline bool is_non_negative(const _T_& x) {
    return is_zero(x) || is_positive(x);
  }

  /// is x <= 0?
  template<class _T_> static inline bool is_non_positive(const _T_& x) {
    return is_zero(x) || is_negative(x);
  }

  template <class T, bool is_signed>
  struct absolute_value_helper {
  };

  template <class T>
  struct absolute_value_helper<T, true> {
     static inline T value(const T& x) { return is_positive(x) ? x : -x; }
  };

  template <class T>
  struct absolute_value_helper<T, false> {
     static inline T value(const T& x) { return x; }
  };

  /// x if x >= 0, -x otherwise 
  template<class _T_> static inline _T_ absolute_value(const _T_& x) {
    return absolute_value_helper<_T_, std::numeric_limits<_T_>::is_signed>().value(x);
  }

  /// 1 if x > 0, -1 if x < 0, 0 otherwise 
  template<class _T_> static inline _T_ sign(const _T_& x) {
    return 
      is_positive(x) ? 
        (scalar_math<_T_>::one()) : 
        (is_negative(x) ? -(scalar_math<_T_>::one()) : (scalar_math<_T_>::zero()));
  }

  /// 1 if x >= 0, -1 if otherwise 
  template<class _T_> static inline _T_ non_zero_sign(const _T_& x) {
    return 
      is_negative(x) ? 
        -(scalar_math<_T_>::one()) : 
         (scalar_math<_T_>::one());
  }

  /// The reciprocal of x, i.e. 1/x
  template<class _T_> static inline _T_ reciprocal(const _T_& x) {
    SL_REQUIRE("Nonzero", !is_zero(x));
    return scalar_math<_T_>::one() / x;
  }

  /// The reciprocal of sqrt(x), i.e. 1/sqrt(x)
  template<class _T_> static inline _T_ reciprocal_sqrt(const _T_& x) {
    SL_REQUIRE("Positive", is_positive(x));
    return scalar_math<_T_>::one() / std::sqrt(x);
  }
  

#define SL_DECLARE_NUMERIC_PROPERTY_GF(_fn_) \
  template<class _T_> static inline _T_ _fn_ (const _T_ &) { return scalar_math<_T_> :: _fn_ (); }
#define SL_DECLARE_NUMERIC_PROPERTY_GF_2(_fn_,_TRES_) \
  template<class _T_> static inline _TRES_ _fn_ (const _T_ &) { return scalar_math<_T_> :: _fn_ (); }
  
  SL_DECLARE_NUMERIC_PROPERTY_GF(zero)
  SL_DECLARE_NUMERIC_PROPERTY_GF(one)
  SL_DECLARE_NUMERIC_PROPERTY_GF(two)
  SL_DECLARE_NUMERIC_PROPERTY_GF(E)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Log2_E)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Log10_E)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Ln_2)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Ln_10)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Pi)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Pi_Over_2)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Pi_Over_4)
  SL_DECLARE_NUMERIC_PROPERTY_GF(One_Over_Pi)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Two_Over_Pi)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Two_Over_Sqrt_Pi)
  SL_DECLARE_NUMERIC_PROPERTY_GF(Sqrt_2)
  SL_DECLARE_NUMERIC_PROPERTY_GF(One_Over_Sqrt_2)
  SL_DECLARE_NUMERIC_PROPERTY_GF(finite_upper_bound)
  SL_DECLARE_NUMERIC_PROPERTY_GF(finite_lower_bound)
  SL_DECLARE_NUMERIC_PROPERTY_GF(upper_bound)
  SL_DECLARE_NUMERIC_PROPERTY_GF(lower_bound)
  SL_DECLARE_NUMERIC_PROPERTY_GF(epsilon)
  SL_DECLARE_NUMERIC_PROPERTY_GF(eta)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(basic_op_max_error_ulp,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF(basic_op_max_error)
  SL_DECLARE_NUMERIC_PROPERTY_GF(basic_op_zero_test)
  SL_DECLARE_NUMERIC_PROPERTY_GF(one_minus_basic_op_max_error)
  SL_DECLARE_NUMERIC_PROPERTY_GF(one_plus_basic_op_max_error)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(digits,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(digits10,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(max_exponent,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(min_exponent,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(max_exponent10,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(min_exponent10,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(precision,int)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(is_signed,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(is_integer,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(is_exact,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF(round_error)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(has_infinity,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(has_quiet_NaN,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(has_signaling_NaN,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(has_denorm,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(has_denorm_loss,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF(infinity)
  SL_DECLARE_NUMERIC_PROPERTY_GF(quiet_NaN)
  SL_DECLARE_NUMERIC_PROPERTY_GF(signaling_NaN)
  SL_DECLARE_NUMERIC_PROPERTY_GF(denorm_min)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(is_bounded,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(is_modulo,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(traps,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(tinyness_before,bool)
  SL_DECLARE_NUMERIC_PROPERTY_GF_2(round_style,std::float_round_style)

} // namespace sl

// --------------------------------------------------------------------------
// -- Extra Math operations for general types
// --------------------------------------------------------------------------

namespace sl {

  /// x raised to the n-th integer power
  template <class T>
  inline T ipow(const T& x, int n) {
    int  e    = (n >= 0) ? (n) : (-n);
    T result = scalar_math<T>::one();
    T factor = x;
    while (e > 0) {
      if ((e % 2) == 1) {
	result *= factor;
      }
      e /= 2;
      factor *= factor;
    }
    return (n >= 0) ? result : reciprocal(result);
  }
  
} // namespace sl

// --------------------------------------------------------------------------
// -- Extra Math operations for basic types
// --------------------------------------------------------------------------

namespace sl {

  /// Absolute value
  inline double abs(double x) {
    return x >= 0.0 ? x : -x;
  }

  /// Remainder of x/y
  inline double mod(double x, double y) {
    return std::fmod(x,y);
  }

  /// the integral value between zero and x which is closest to x
  inline double integer_round_towards_zero(double x) {
    return (x>0.0) ? std::floor(x) : std::ceil(x);
  }

  /// x^2
  inline double sqr(double x) {
    return x*x;
  }

}; // namespace sl

namespace sl {

  /// Absolute value
  inline float abs(float x) {
    return x >= 0.0f ? x : -x;
  }

  /// Remainder of x/y
  inline float mod(float x, float y) {
    return std::fmod(x,y);
  }

  /// the integral value between zero and x which is closest to x
  inline float integer_round_towards_zero(float x) {
    return (x>0.0f) ? std::floor(x) : std::ceil(x);
  }

  // x^2
  inline float sqr(float x) {
    return x*x;
  }

}; // namespace sl

#endif
