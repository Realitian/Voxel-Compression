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
#ifndef SL_OPERATORS_HPP
#define SL_OPERATORS_HPP

#include <sl/config.hpp>

// -- Macros to define global operators in terms of a set 
// -- of fundamental ones. These are *not* defined in classes
// -- to be inherited because current C++ compilers would
// -- cause the resulting classes to be much larger than they 
// -- should be (inheriting from empty classes imposes
// -- has an overhead!)

// -- Because of a bug in overloading resolution (in at least
// -- g++), we do not use friend declarations

#include <sl/cast.hpp>

#define SL_OP_NO_TEMPLATE 
#define SL_OP_TEMPLATE1(A) A
#define SL_OP_TEMPLATE2(A,B) A , B
#define SL_OP_TEMPLATE3(A,B,C) A , B , C
#define SL_OP_TEMPLATE4(A,B,C,D) A , B , C , D
#define SL_OP_TEMPLATE5(A,B,C,D,E) A , B ,  C , D ,  E


#define SL_OP_COMPARABLE1(T) \
     inline bool operator>(const T& y)  const { return y < *(static_cast<const T*>(this)); } \
     inline bool operator<=(const T& y) const { return !(y < *(static_cast<const T*>(this))); } \
     inline bool operator>=(const T& y) const { return !(*(static_cast<const T*>(this)) < y); } 

#define SL_OP_COMPARABLE2(T, U) \
     inline bool operator<=(const U& y) const { return !(*(static_cast<const T*>(this)) > y); } \
     inline bool operator>=(const U& y) const { return !(*(static_cast<const T*>(this)) < y); }

#define SL_OP_COMPARABL2_OVERLOADS(TEMPLATE_DECL,T,U) \
     TEMPLATE_DECL inline bool operator>(const U& x, const T& y)  { return y < x; } \
     TEMPLATE_DECL inline bool operator<(const U& x, const T& y)  { return y > x; } \
     TEMPLATE_DECL inline bool operator<=(const U& x, const T& y) { return !(y < x); } \
     TEMPLATE_DECL inline bool operator>=(const U& x, const T& y) { return !(y > x); }


#define SL_OP_EQUALITY_COMPARABLE1(T) \
     inline bool operator!=(const T& y) const { return !(*(static_cast<const T*>(this)) == y); }

#define SL_OP_EQUALITY_COMPARABLE2(T, U) \
     inline bool operator!=(const T& y, const U& x) const { return !(y == x); }

#define SL_OP_EQUALITY_COMPARABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \
     TEMPLATE_DECL inline bool operator==(const U& y, const T& x) { return x == y; } \
     TEMPLATE_DECL inline bool operator!=(const U& y, const T& x) { return !(x == y); } 


#define SL_OP_MULTIPLIABLE1(T) \
     inline T operator*(const T& y) const { T x(*(static_cast<const T*>(this))); return x *= y; }

#define SL_OP_MULTIPLIABLE2(T, U) \
     inline T operator*(const U& y) const { T x(*(static_cast<const T*>(this))); return x *= y; } 

// Check: commutative?
#define SL_OP_MULTIPLIABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \
     TEMPLATE_DECL inline T operator*(const U& y, T x) { return x *= y; }


#define SL_OP_ADDABLE1(T) \
     inline T operator+(const T& y) const { T x(*(static_cast<const T*>(this))); return x += y; }

#define SL_OP_ADDABLE2(T,U) \
     inline T operator+(const U& y) const { T x(*(static_cast<const T*>(this))); return x += y; }

#define SL_OP_ADDABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \
     TEMPLATE_DECL inline T operator+(const U& y, T x) { return x += y; }


#define SL_OP_SUBTRACTABLE1(T) \
     inline T operator-(const T& y) const { T x(*(static_cast<const T*>(this))); return x -= y; }

#define SL_OP_SUBTRACTABLE2(T, U) \
     inline T operator-(const U& y) const { T x(*(static_cast<const T*>(this))); return x -= y; }

#define SL_OP_SUBTRACTABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \


#define SL_OP_DIVIDABLE1(T) \
     inline T operator/(const T& y) const { T x(*(static_cast<const T*>(this))); return x /= y; }

#define SL_OP_DIVIDABLE2(T,U) \
     inline T operator/(const U& y) const { T x(*(static_cast<const T*>(this))); return x /= y; }

#define SL_OP_DIVIDABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \


#define SL_OP_MODABLE1(T) \
     inline T operator%(const T& y) const { T x(*(static_cast<const T*>(this))); return x %= y; }

#define SL_OP_MODABLE2(T, U) \
     inline T operator%(const U& y) const { T x(*(static_cast<const T*>(this))); return x %= y; }

#define SL_OP_MODABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \


#define SL_OP_XORABLE1(T) \
     inline T operator^(const T& y) const { T x(*(static_cast<const T*>(this))); return x ^= y; }

#define SL_OP_XORABLE2(T, U) \
     inline T operator^(const U& y) const { T x(*(static_cast<const T*>(this))); return x ^= y; } 

#define SL_OP_XORABLE2_OVERLOADS(TEMPLATE_DECL,T, U) \
     TEMPLATE_DECL inline T operator^(const U& y, T x) { return x ^= y; }

#define SL_OP_ANDABLE1(T) \
     inline T operator&(const T& y) const { T x(*(static_cast<const T*>(this))); return x &= y; }

#define SL_OP_ANDABLE2(T,U) \
     inline T operator&(const U& y) const { T x(*(static_cast<const T*>(this))); return x &= y; } 

#define SL_OP_ANDABLE2_OVERLOADS(TEMPLATE_DECL,T,U) \
     TEMPLATE_DECL inline T operator&(const U& y, T x) { return x &= y; }


#define SL_OP_ORABLE1(T) \
     inline T operator|(const T& y) const { T x(*(static_cast<const T*>(this))); return x |= y; }

#define SL_OP_ORABLE2(T,U) \
     inline T operator|(const U& y) const { T x(*(static_cast<const T*>(this))); return x |= y; }

#define SL_OP_ORABLE2_OVERLOADS(TEMPLATE_DECL,T,U) \
     TEMPLATE_DECL inline T operator|(const U& y, T x) { return x |= y; }


#define SL_OP_INCREMENTABLE(T) \
     inline T operator++(int) { T tmp(*(static_cast<const T*>(this))); ++(*this); return tmp; }

#define SL_OP_DECREMENTABLE(T) \
     inline T operator--(int) { T tmp(*(static_cast<const T*>(this))); --(*this); return tmp; }

//
// -- Combinations of operators defining common algebraic
// -- structures
//

// Algebraic structures over a single set
//
// MULTIPLICATIVE_SEMIGROUP: Associative multiplicative law (operator *)
// MULTIPLICATIVE_MONOID   : SEMIGROUP + Identity (one)
// MULTIPLICATIVE_GROUP    : MONOID + Inverse (operator /)
//
// ADDITIVE_SEMIGROUP      : Associative + Commutative law (operator +)
// ADDITIVE_MONOID         : ADDITIVE_SEMIGROUP + Identity (zero)
// ADDITIVE_GROUP          : ADDITIVE_MONOID + Inverse (operator -)
//
// RING                    : MULTIPLICATIVE_SEMIGROUP + ADDITIVE_GROUP
// RING_WITH_UNIT          : MULTIPLICATIVE_MONOID    + ADDITIVE_GROUP
// FIELD                   : MULTIPLICATIVE_GROUP     + ADDITIVE_GROUP

#define SL_OP_MULTIPLICATIVE_SEMIGROUP(T) \
  SL_OP_MULTIPLIABLE1(T)

#define SL_OP_MULTIPLICATIVE_MONOID(T) \
  SL_OP_MULTIPLIABLE1(T)

#define SL_OP_MULTIPLICATIVE_GROUP(T) \
  SL_OP_MULTIPLIABLE1(T) \
  SL_OP_DIVIDABLE1(T)

#define SL_OP_ADDITIVE_SEMIGROUP(T) \
  SL_OP_ADDABLE1(T)

#define SL_OP_ADDITIVE_MONOID(T) \
  SL_OP_ADDABLE1(T)

#define SL_OP_ADDITIVE_GROUP(T) \
  SL_OP_ADDABLE1(T) \
  SL_OP_SUBTRACTABLE1(T)

#define SL_OP_RING(T) \
  SL_OP_MULTIPLICATIVE_SEMIGROUP(T) \
  SL_OP_ADDITIVE_GROUP(T)

#define SL_OP_RING_WITH_UNIT(T) \
  SL_OP_MULTIPLICATIVE_MONOID(T) \
  SL_OP_ADDITIVE_GROUP(T)

#define SL_OP_FIELD(T) \
  SL_OP_MULTIPLICATIVE_GROUP(T) \
  SL_OP_ADDITIVE_GROUP(T)

// Algebraic structures over two sets
//
// LINEAR_SPACE            : Vectors: ADDITIVE_GROUP Scalars: FIELD
// LINEAR_ALGEBRA          : Vectors: RING           Scalars: FIELD
// LINEAR_ALGEBRA_WITH_UNIT: Vectors: RING_WITH_UNIT Scalars: FIELD
// DIVISION_ALGEBRA        : Vectors: FIELD          Scalars: FIELD

#define SL_OP_LINEAR_SPACE(V,S) \
  SL_OP_ADDITIVE_GROUP(V) \
  SL_OP_MULTIPLIABLE2(V,S) \
  SL_OP_DIVIDABLE2(V,S)

#define SL_OP_LINEAR_SPACE_OVERLOADS(TEMPLATE_DECL,V,S)  \
     TEMPLATE_DECL inline V operator*(const S& y, V x) { return x *= y; }


#define SL_OP_LINEAR_ALGEBRA(V,S) \
  SL_OP_RING(V) \
  SL_OP_MULTIPLIABLE2(V,S) \
  SL_OP_DIVIDABLE2(V,S)

#define SL_OP_LINEAR_ALGEBRA_OVERLOADS(TEMPLATE_DECL,V,S) \
     TEMPLATE_DECL inline V operator*(const S& y, V x) { return x *= y; }

#define SL_OP_LINEAR_ALGEBRA_WITH_UNIT(V,S) \
  SL_OP_RING_WITH_UNIT(V) \
  SL_OP_MULTIPLIABLE2(V,S) \
  SL_OP_DIVIDABLE2(V,S)

#define SL_OP_LINEAR_ALGEBRA_WITH_UNIT_OVERLOADS(TEMPLATE_DECL,V,S) \
     TEMPLATE_DECL inline V operator*(const S& y, V x) { return x *= y; }


#define SL_OP_DIVISION_ALGEBRA(V,S) \
  SL_OP_FIELD(V) \
  SL_OP_MULTIPLIABLE2(V,S) \
  SL_OP_DIVIDABLE2(V,S)

#define SL_OP_DIVISION_ALGEBRA_OVERLOADS(TEMPLATE_DECL,V,S) \
     TEMPLATE_DECL inline V operator*(const S& y, V x) { return x *= y; }


#endif


