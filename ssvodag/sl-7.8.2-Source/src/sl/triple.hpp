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
#ifndef SL_TRIPLE_HPP
#define SL_TRIPLE_HPP

#include <sl/config.hpp>
#include <sl/serializer.hpp>
#include <sl/hash.hpp>

namespace sl {

  /** 
   * Triple of values, base on a straightforward 
   * modification of std::pair
   */
  template <class T1, class T2, class T3>
  struct triple {
    typedef T1 first_type;    
    typedef T2 second_type;   
    typedef T3 third_type;    
 
    T1 first;                 
    T2 second;                
    T3 third;                

    triple() : first(T1()), second(T2()), third (T3()) {}
 
    triple(const T1& a, const T2& b, const T3& c) : first(a), second(b), third(c) {}
    
    template <class _U1, class _U2, class _U3>
    triple(const triple<_U1, _U2, _U3>& p) : first(p.first), second(p.second), third(p.third) {}
  };
 
  template <class T1, class T2, class T3>
  inline bool operator==(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) { 
    return x.first == y.first && x.second == y.second && x.third == y.third;
  }
 
  template <class T1, class T2, class T3>
  inline bool operator<(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) { 
    return x.first < y.first || 
      (!(y.first < x.first) && x.second < y.second) ||
      ((!(y.first < x.first) && !(y.second < x.second) && x.third < y.third)); 
  }
 
  template <class T1, class T2, class T3>
  inline bool operator!=(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) {
    return !(x == y);
  }
 
  template <class T1, class T2, class T3>
  inline bool operator>(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) {
    return y < x;
  }
 
  template <class T1, class T2, class T3>
  inline bool operator<=(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) {
    return !(y < x);
  }
 
  template <class T1, class T2, class T3>
  inline bool operator>=(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) {
    return !(x < y);
  }
 
  template <class T1, class T2, class T3>
  inline triple<T1, T2, T3> make_triple(const T1& x, const T2& y, const T3& z) {
    return triple<T1, T2, T3>(x, y, z);
  }

  template <class A, class B, class C>
  class tied_triple {
  protected:
    A& a_;
    B& b_;
    C& c_;
  public:
    inline tied_triple(A& a, B& b, C& c) : a_(a), b_(b), c_(c) { 
    }

    template <class U, class V, class W>
    inline tied_triple& operator=(const sl::triple<U,V,W>& p) {
      a_ = p.first;
      b_ = p.second;
      c_ = p.third;
      return *this;
    }
  };
 
  /**
   *  Tie three values in a triple (inspired by boost).
   *  This is a utility function that makes it more convenient to 
   *  work with a function which returns a sl::triple<>.
   */
  template <class A, class B, class C>
  inline tied_triple<A,B,C> tie(A& a, B& b, C& c) { 
    return tied_triple<A,B,C>(a, b,c); 
  } 

  /// Serialization traits for sl::triple container
  template <class T1, class T2, class T3>
  struct serialization_traits< sl::triple<T1, T2,T3> > {
    typedef sl::triple<T1, T2, T3> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return s << x.first << x.second << x.third;
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return s >> x.first >> x.second >> x.third;
    }
  };

  /// Hashing 
  template <class T1, class T2, class T3>
  inline std::size_t hash_value(sl::triple<T1, T2, T3> const& v) {
    std::size_t seed = 0;
    hash_combine(seed, v.first);
    hash_combine(seed, v.second);
    hash_combine(seed, v.third);
    return seed;
  }

} // namespace sl

#endif
