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
#ifndef SL_UTILITY_HPP
#define SL_UTILITY_HPP

#include <sl/assert.hpp>
#include <sl/time_duration.hpp>
#include <algorithm>
#include <string>    // std::string
#include <sstream>   // std::to_string
#include <cstdio>    // sprintf
#include <sl/cstdint.hpp>   

namespace sl {

  namespace tags {

    /**
     *  An empty struct representing the value initialized, for use
     *  in template specializations.
     */
    struct initialized { };

    /**
     *  An empty struct representing the value not initialized, for use
     *  in template specializations.
     */
    struct not_initialized { };

    /**
     *  An empty struct representing the value shared, for use
     *  in template specializations.
     */
    struct shared_state { };
    
    /**
     *  An empty struct representing the value not shared, for use
     *  in template specializations.
     */
    struct not_shared_state { };

  };

} // namespace sl

namespace sl {

  /// The characters used for directory separation in path names
  extern std::string pathname_directory_separators();

  /// The characters used for extension separation in path names
  extern std::string pathname_extension_separators();

  /// Is c a directory separator?
  extern bool pathname_is_directory_separator(const char& c);

  /// Is s an absolute path?
  extern bool pathname_is_absolute(const std::string& s);

  /// The path name before last directory separator
  extern std::string pathname_directory(const std::string& s);

  /// The path name after last directory separator
  extern std::string pathname_base(const std::string& s);

  /// The base name after the last extension separator
  extern std::string pathname_extension(const std::string& s);

  /// The path name before the last '.'
  extern std::string pathname_without_extension(const std::string& s);
  
  /**
   * Matches a string with a pattern. Special characters:
   *	    '*'  Matches any string, including the null string.
   *        '?'  Matches any single character.
   *        '\'  Escapes next special character
   */
  extern bool matches(const char* str, const char* regexp);
  
  /**
   * Matches a string with a pattern. Special characters:
   *	    '*'  Matches any string, including the null string.
   *        '?'  Matches any single character.
   *        '\'  Escapes next special character
   */
  extern bool matches(const std::string& str, const std::string& regexp);

  /** 
   *  The object of type T represented by string s. Type
   *  T must have an operator >> on streams defined.
   *
   *  Example: float x = from_string<float>("1.23");
   */
  template<class T>
  extern T from_string(const std::string& s) {
    std::istringstream is(s);
    T t;
    is >> t;
    return t;
  }

  /**
   *  The string representation of t. Type T
   *  must have an operator << on streams defined.
   *
   *  Example std::string s = to_string(1.25);
   */
  template<class T>
  extern std::string to_string(const T& t) {
    std::ostringstream s;
    s << t;
    return s.str();
  }

  /**
   * Split string into components separated by sep
   * Example:
   *  std::vector<std::string> words;
   *  sl::split(words, " \t.,;:!?", std::back_inserter(sentence));
   */
  template<class OutIt>
  void split(const std::string& s, const std::string& sep, OutIt dest) {
    std::string::size_type left = s.find_first_not_of( sep );
    std::string::size_type right = s.find_first_of( sep, left );
    while( left < right ) {
      *dest = s.substr( left, right-left );
      ++dest;
      left = s.find_first_not_of( sep, right );
      right = s.find_first_of( sep, left );
    }
  }

} // namespace sl


namespace sl {
  
  extern std::string human_readable_duration(const time_duration& t);

  extern std::string human_readable_quantity(const uint64_t q,
					     const char*    unit="",
					     const uint64_t K=1000); 
  
  extern std::string human_readable_size(const uint64_t sz);

  extern std::string human_readable_percent(const double percent);

}

namespace sl {

#ifdef _WIN32
#undef min
#undef max
#endif

#if HAVE_CXX_STD_MINMAX

  using std::min;
  using std::max;

#else

  /**
   *  The minimum of two values.
   */
  template <class T>
  inline const T& min(const T& a, const T& b) {
    if (a < b) {
      return a;
    } else {
      return b;
    }
  }

  /**
   *  The maximum of two values.
   */
  template <class T>
  inline const T& max(const T& a, const T& b) {
    if (a < b) {
      return b;
    } else {
      return a;
    }
  }

#endif

  /**
   *  The median of three values.
   */
  template <class T>
  inline const T& median(const T& a, const T& b, const T& c) {
    if (a < b)
      if (b < c)
	return b;
      else if (a < c)
	return c;
      else
	return a;
    else if (a < c)
      return a;
    else if (b < c)
      return c;
    else
      return b;
  }                                                                                                                         
#if 1
  using std::swap;
#else
  /**
   *  Swap the values referenced by the parameters.
   */
  template <class T> 
  inline void swap(T& a, T& b) {
    T tmp = a;
    a = b;
    b = tmp;
  }
#endif

  /**
   *  Greatest common divisor of a and b
   */
  template <class T>
  inline T gcd(T a, T b) {
    SL_REQUIRE("Positive a", a>0);
    SL_REQUIRE("Positive b", b>0);
    do {
      const T tmp(b);
      b = a % b;
      a = tmp;
    } while (b != T(0));
    
    return a;
  }

  /**
   *  Greatest common multiple of a and b
   */
  template <class T>
  inline T lcm(const T& a, const T& b) {
    SL_REQUIRE("Positive a", a>0);
    SL_REQUIRE("Positive b", b>0);
    
    return (a / gcd(a,b)) * b;
  }
  
  /*
   * Utility class for binding two individual
   * variables in a single std::pair<>. Inspired
   * by boost.
   */
  template <class A, class B>
  class tied {
  public:
    inline tied(A& a, B& b) : a_(a), b_(b) { 
    }

    template <class U, class V>
    inline tied& operator=(const std::pair<U,V>& p) {
      a_ = p.first;
      b_ = p.second;
      return *this;
    }
  protected:
    A& a_;
    B& b_;
  };


  /**
   *  Emulation of C++11 addressof operator to obtains the actual address of the
   *  object or function arg, even in presence of overloaded operator&
   */
  template <class T>
  inline T* addressof(T& arg) {
    return reinterpret_cast<T*>(&const_cast<char&>(reinterpret_cast<const volatile char&>(arg)));
  }

  /**
   *  Tie two values in a pair (inspired by boost).
   *  This is a utility function that makes it more convenient to 
   *  work with a function which returns a std::pair<>.
   */
  template <class A, class B>
  inline tied<A,B> tie(A& a, B& b) { 
    return tied<A,B>(a, b); 
  } 

  /**
   *  The cyclic successor of x in the interval [x_first, x_last]
   */
  template <class T>
  inline T cyclic_succ(T x, const T& x_first, const T& x_last) {
    if (x == x_last) {
      x = x_first;
    } else {
      ++x;
    }
    return x;
  }

  /**
   *  The cyclic predecessor of x in the interval [x_first, x_last]
   */
  template <class T>
  inline T cyclic_pred(T x, const T& x_first, const T& x_last) {
    if (x == x_first) {
      x = x_last;
    } else {
      --x;
    }
    return x;
  }

} // namespace sl

namespace sl {

  /// The i-th prime number (0<=i<=255)
  extern std::size_t small_prime(const unsigned char i);

}

namespace sl {
 
#define SL_DECLARE_GENERIC_SUPERCLASS_DATA(G_USERDEFINED, NAME_USERDEFINED) \
  G_USERDEFINED NAME_USERDEFINED

#define SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(G_DERIVED) \
  inline const G_DERIVED* derived() const {  \
    return static_cast<const G_DERIVED*>(this);  \
  }  \
  inline G_DERIVED* derived() {  \
    return static_cast<G_DERIVED*>(this);  \
  }  \
  inline const G_DERIVED& derived_ref() const {  \
    return *(static_cast<const G_DERIVED*>(this));  \
  }  \
  inline G_DERIVED& derived_ref() {  \
    return *(static_cast<G_DERIVED*>(this)); \
  }
  
  /** 
   *  Generic superclass interface
   */
  template <class I, class W>
  class generic_superclass {
  public:
    /// The type of user defined data
    typedef W                           user_defined_data_t;
    /// The type of the derived class
    typedef I                           derived_t;

  protected:

    /// User defined data. Declared here to remove the space overhead associated to empty base classes.
    SL_DECLARE_GENERIC_SUPERCLASS_DATA(user_defined_data_t,data_);

  public:

    /// Default init - no operation by default
    inline generic_superclass() {}

    /// User-defined data init
    inline generic_superclass(const user_defined_data_t& data): data_(data) {}

  public:
    
    /// Features to access this as a derived_t pointer
    SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);

  };

} // namespace sl

namespace sl {

  /** 
   * A simple visitor counting the number of times
   * operator() is called
   */
  class counting_visitor {
  protected:
    std::size_t count_;
  public:
    
    inline counting_visitor(): count_(0) {
    }
    
    template <typename T1>
    inline void operator()(const T1&) {
      ++count_;
    }

    template <typename T1, typename T2>
    inline void operator()(const T1&, const T2&) {
      ++count_;
    }

    template <typename T1, typename T2, typename T3>
    inline void operator()(const T1&, const T2&, const T3&) {
      ++count_;
    }

    template <typename T1, typename T2, typename T3, typename T4>
    inline void operator()(const T1&, const T2&, const T3&, const T4&) {
      ++count_;
    }

    inline std::size_t count() const {
      return count_;
    }
  }; // class counting_visitor

}

#endif



