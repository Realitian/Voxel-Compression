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
#ifndef SL_FASTEST_HPP
#define SL_FASTEST_HPP

#include <sl/utility.hpp>
#include <sl/operators.hpp>

namespace sl {

  /**
   * Fixed Algorithm Size Template Extended Software Toolkit 
   *  
   * This library provides a number of linear algorithms 
   * that are customized for sizes known at compile time.
   * The design is inspired by the FAST library by 
   * Lumsdaine et al., but the implementation is considerably
   * different. In particular:
   * template classes rather than template functions,
   * to avoid the need for creating zero-sized objects
   * at run-time;
   * template expansion based on binary subdivision
   * (rather than linear). This makes it possible to
   * supports sizes up to 2^N (rather than N, where N is
   * the compiler limit (16 for ISO C++);
   * customizable loop unrolling behavior, to reduce
   * code bloat and maximize cache coherence.
   *
   * If SL_NO_TEMPLATE_METAPROGRAMS is nonzero, the
   * implementation forces the unroll level to 1 (i.e. 
   * reverts to simple loops).
   */
  namespace fastest {

#   define SL_FASTEST_HEAD(N_)          ((((N_)/2)>=1)?((N_)/2):(N_))
#   define SL_FASTEST_TAIL(N_)          ((((N_)/2)>=1)?((N_)-(N_)/2):(0))
#   define SL_FASTEST_REST(N_,NMAX_)    ((N_)%(NMAX_))
#   define SL_FASTEST_BATCH(N_,NMAX_)   ((N_)>=(NMAX_)?(NMAX_):(0))

#if SL_NO_TEMPLATE_METAPROGRAMS

#   define SL_FASTEST_APPLY_TAIL_RECURSION(CLASS_ID_,N_,ROUTINE_CALL_) \

#   define SL_FASTEST_APPLY_BATCH_RECURSION(CLASS_ID_,N_,ROUTINE_CALL_) \
	{ \
          for (size_t __loop_count__ = (N_); __loop_count__ != 0; --__loop_count__) { \
	    ROUTINE_CALL_ ; \
	  } \
	}

    /// Max unroll level for loops
    const size_t default_max_unroll_level = 1;

#else

#   define SL_FASTEST_APPLY_TAIL_RECURSION(CLASS_ID_,N_,ROUTINE_CALL_) \
	{ \
          const size_t HEAD_ITERS_     = SL_FASTEST_HEAD(N_); \
          const size_t TAIL_ITERS_     = SL_FASTEST_TAIL(N_); \
	  ::sl::fastest:: CLASS_ID_ < HEAD_ITERS_, (HEAD_ITERS_ == 0 ? 0 : HEAD_ITERS_) >:: ROUTINE_CALL_ ; \
          ::sl::fastest:: CLASS_ID_ < TAIL_ITERS_, (TAIL_ITERS_ == 0 ? 0 : TAIL_ITERS_) >:: ROUTINE_CALL_ ; \
        }

#   define SL_FASTEST_APPLY_BATCH_RECURSION(CLASS_ID_,N_,ROUTINE_CALL_) \
	{ \
          const size_t BATCH_ITERS_    = SL_FASTEST_BATCH(N_,N_MAX); \
          const size_t REST_ITERS_     = SL_FASTEST_REST(N_,N_MAX); \
          for (size_t __loop_count__ = (N_)/N_MAX; __loop_count__ != 0; --__loop_count__) { \
	    ::sl::fastest:: CLASS_ID_ < BATCH_ITERS_, (BATCH_ITERS_ == 0 ? 0 : BATCH_ITERS_) >:: ROUTINE_CALL_ ; \
	  } \
	  ::sl::fastest:: CLASS_ID_ < REST_ITERS_, (REST_ITERS_ == 0 ? 0 : REST_ITERS_) >:: ROUTINE_CALL_ ; \
	}

    /// Max unroll level for loops
    const size_t default_max_unroll_level = 16;

#endif

    /// Copy(x,y): y <- x
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class copy {
    public:
      template <class InIter, class OutIter>
      static inline void apply_sub(InIter& first, OutIter& result) {
	*result = *first;
	result++; first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(copy,N-1,apply_sub(first, result));
      }
    public:
      template <class InIter, class OutIter>
      static inline void apply(InIter first, OutIter result) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  result[i] = first[i];
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(copy,N,apply_sub(first, result));
#endif
      }

    }; // copy<N>
    
    template <>
    class copy<0,0> {
    public:
      template <class InIter, class OutIter>
      static inline void apply_sub(InIter&, OutIter&) {
      }
    public: 
      template <class InIter, class OutIter>
      static inline void apply(InIter, OutIter result) {
        SL_USEVAR(result);
      }
    }; // copy<0,0>


    /// Self-transform: x = op(x) or y = op(x,y) or y = op(y,x)  
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class self_transform {
    public:
      template <class InOutIter, class UnaryOp>
      static inline void apply_sub(InOutIter& first_and_result, UnaryOp op) {
	*first_and_result = op(*first_and_result);
	first_and_result++;
	SL_FASTEST_APPLY_TAIL_RECURSION(self_transform,N-1,apply_sub(first_and_result, op));
      }
      template <class InIter1, class InOutIter, class BinaryOp>
      static inline void apply_sub(InIter1& first1, InOutIter& first2_and_result, BinaryOp op) {
	*first2_and_result = op(*first1, *first2_and_result);
	first1++; first2_and_result++;
	SL_FASTEST_APPLY_TAIL_RECURSION(self_transform,N-1,apply_sub(first1, first2_and_result, op));
      }

    public:
      /// Self unary op: x = op(x)
      template <class InOutIter, class UnaryOp>
      static inline void apply(InOutIter first_and_result, UnaryOp op) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  first_and_result[i] = op(first_and_result[i]);
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(self_transform,N,apply_sub(first_and_result, op));
#endif
      }
      /// Self Binary op: y = op(x,y)
      template <class InIter1, class InOutIter, class BinaryOp>
      static inline void apply(InIter1 first1, InOutIter first2_and_result, BinaryOp op) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  first2_and_result[i] = op(first1[i], first2_and_result[i]);
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(self_transform,N,apply_sub(first1, first2_and_result, op));
#endif
      }
    }; // self_transform<N>

    template <>
    class self_transform<0,0> {
    public:
      template <class InOutIter, class UnaryOp>
      static inline void apply_sub(InOutIter&, UnaryOp) {
      }
      template <class InIter1, class InOutIter, class BinaryOp>
      static inline void apply_sub(InIter1&, InOutIter&, BinaryOp) {
      }
    public:
      // Self unary op
      template <class InOutIter, class UnaryOp>
      static inline void apply(InOutIter first_and_result, UnaryOp) {
        SL_USEVAR(first_and_result);
      }
      // Self binary op
      template <class InIter, class InOutIter, class BinaryOp>
      static inline void apply(InIter, InOutIter first2_and_result, BinaryOp) {
        SL_USEVAR(first2_and_result);
      }
    }; // self_transform<0,N_MAX>

    /// Transform: y = op(x) or z= op(x,y)
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class transform {
    public:
      template <class InIter, class OutIter, class UnaryOp>
      static inline void apply_sub(InIter& first, OutIter& result, UnaryOp op) {
	*result = op(*first);
	result++; first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(transform,N-1,apply_sub(first, result, op));
      }

      template <class InIter1, class InIter2, class OutIter, class BinaryOp>
      static inline void apply_sub(InIter1& first1, InIter2& first2, OutIter& result, BinaryOp op) {
	*result = op(*first1, *first2);
	result++; first1++; first2++;
	SL_FASTEST_APPLY_TAIL_RECURSION(transform,N-1,apply_sub(first1, first2, result, op));
      }

    public:
      // Unary op
      template <class InIter, class OutIter, class UnaryOp>
      static inline void apply(InIter first, OutIter result, UnaryOp op) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  result[i] = op(first[i]);
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(transform,N,apply_sub(first, result,op));
#endif
      }

      // Binary op
      template <class InIter1, class InIter2, class OutIter, class BinaryOp>
      static inline void apply(InIter1 first1, InIter2 first2, OutIter result, BinaryOp op) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  result[i] = op(first1[i],first2[i]);
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(transform,N,apply_sub(first1,first2, result,op));
#endif
      }
    }; // transform<N>

    template <>
    class transform<0,0> {
    public:
      template <class InIter, class OutIter, class UnaryOp>
      static inline void apply_sub(InIter&, OutIter&, UnaryOp) {
      }
      template <class InIter1, class InIter2, class OutIter, class BinaryOp>
      static inline void apply_sub(InIter1&, InIter2&, OutIter&, BinaryOp) {
      }
    public:
      // Unary op
      template <class InIter, class OutIter, class UnaryOp>
      static inline void apply(InIter, OutIter result, UnaryOp) {
        SL_USEVAR(result);
      }
      // Binary op
      template <class InIter1, class InIter2, class OutIter, class BinaryOp>
      static inline void apply(InIter1, InIter2, OutIter result, BinaryOp) {
        SL_USEVAR(result);
      }
    }; // transform<0,0>


    /// Fill(x,v): x <- v
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class fill {
    public:
      template <class OutIter, class T>
      static inline void apply_sub(OutIter& first, const T& value) {
	*first = value;
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(fill,N-1,apply_sub(first, value));
      }

    public:
      template <class OutIter, class T>
      static inline void apply(OutIter first, const T& value) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  first[i] = value;
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(fill,N,apply_sub(first, value));
#endif
      }
    }; // fill<N>

    template <>
    class fill<0,0> {
    public:
      template <class OutIter, class T>
      static inline void apply_sub(OutIter&, const T&) {
      }
    public:
      template <class OutIter, class T>
      static inline void apply(OutIter first, const T&) {
        SL_USEVAR(first);
      }
    }; // fill<0,0>

    /// Swap Ranges
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class swap_ranges {
    public:
      template <class Iter1, class Iter2>
      static inline void apply_sub(Iter1& first1, Iter2& first2) {
	sl::swap(*first2,*first1);
	first2++; first1++;
	SL_FASTEST_APPLY_TAIL_RECURSION(swap_ranges,N-1,apply_sub(first1, first2));
      }
  
    public:
      template <class Iter1, class Iter2>
      static inline Iter2 apply(Iter1 first1, Iter2 first2) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  sl::swap(first1[i],first2[i]);
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(swap_ranges,N,apply_sub(first1, first2));
#endif
      }
    }; // swap_ranges<N>

    template <>
    class swap_ranges<0,0> {
    public:
      template <class Iter1, class Iter2>
      static inline void apply_sub(Iter1&, Iter2&) {
      }
    public:
      template <class Iter1, class Iter2>
      static inline Iter2 apply(Iter1, Iter2 first2) {
        SL_USEVAR(first2);
      }
    }; // swap_ranges<0,N_MAX>

    /// Accumulate: res = v0 + sum(x_i) or res = v0 op x0 op x1 ... op xn  
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class accumulate {
    public:
      template <class Iter, class T>
      static inline void apply_sub(Iter& first, T& value) {
	value += *first;
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(accumulate,N-1,apply_sub(first, value));
      }
      template <class Iter, class T, class BinaryOp>
      static inline void apply_sub(Iter& first, T& value, BinaryOp op) {
	value = op(value,*first);
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(accumulate,N-1,apply_sub(first, value, op));
      }
    public:
      template <class Iter, class T>
      static inline T apply(Iter first, T value) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  value += first[i];
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(accumulate,N,apply_sub(first, value));
#endif
	return value;
      }
      // Binary op
      template <class Iter, class T, class BinaryOp>
      static inline T apply(Iter first, T value, BinaryOp op) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  value = op(value,first[i]);
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(accumulate,N,apply_sub(first, value, op));
#endif
	return value;
      }
    }; // accumulate<N>

    template <>
    class accumulate<0,0> {
    public:
      template <class Iter, class T>
      static inline void apply_sub(Iter&, T) {
      }
      template <class Iter, class T, class BinaryOp>
      static inline void apply_sub(Iter&, T, BinaryOp) {
      }
    public:
      // Sum
      template <class Iter, class T>
      static inline T apply(Iter, T value) {
	return value;
      }
      // Binary op
      template <class Iter, class T, class BinaryOp>
      static inline T apply(Iter, T value, BinaryOp) {
	return value;
      }
    }; // accumulate<0,0>

    /// Inner Product = v0 + sum x_i * y_i or v0 op1 (x0 op2 y0) op1 (x1 op2 y1) ... 
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class inner_product {
    public:

      template <class Iter1, class Iter2, class T>
      static inline void apply_sub(Iter1& first1, Iter2& first2, T& value) {
	value += *first1 * *first2;
	first1++; first2++;
	SL_FASTEST_APPLY_TAIL_RECURSION(inner_product,N-1,apply_sub(first1, first2, value));
      }

      template <class Iter1, class Iter2, class T, class BinaryOp1, class BinaryOp2>
      static inline void apply_sub(Iter1& first1, Iter2& first2, T& value, BinaryOp1 op1, BinaryOp2 op2) {
	value = op1 ( value, op2( *first1, *first2) );
	first1++; first2++;
	SL_FASTEST_APPLY_TAIL_RECURSION(inner_product,N-1,apply_sub(first1, first2, value, op1, op2));
      }

    public:
      // Standard ( + / * )
      template <class Iter1, class Iter2, class T>
      static inline T apply(Iter1 first1, Iter2 first2, T value) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  value += first1[i] * first2[i];
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(inner_product,N,apply_sub(first1, first2, value));
#endif
	return value;
      }
      // Binary op1 / op2
      template <class Iter1, class Iter2, class T, class BinaryOp1, class BinaryOp2>
      static inline T apply(Iter1 first1, Iter2 first2, T value, BinaryOp1 op1, BinaryOp2 op2) {
#if SL_NO_TEMPLATE_METAPROGRAMS
	for (size_t i=0; i<N; ++i) {
	  value = op1 ( value, op2( first1[i], first2[i]) );
	}
#else
	SL_FASTEST_APPLY_BATCH_RECURSION(inner_product,N,apply_sub(first1, first2, value, op1, op2));
#endif
	return value;
      }
    }; // inner_product<N>

    template <size_t N_MAX>
    class inner_product<0,N_MAX> {
    public:

      template <class Iter1, class Iter2, class T>
      static inline void apply_sub(Iter1&, Iter2&, T&) {
      }
      template <class Iter1, class Iter2, class T, class BinaryOp1, class BinaryOp2>
      static inline void apply_sub(Iter1&, Iter2&, T&, BinaryOp1, BinaryOp2) {
      }
    
    public:
      // Standard ( + / * )
      template <class Iter1, class Iter2, class T>
      static inline T apply(Iter1, Iter2, T value) {
	return value;
      }
      // Binary op1 / op2
      template <class Iter1, class Iter2, class T, class BinaryOp1, class BinaryOp2>
      static inline T apply(Iter1, Iter2, T value, BinaryOp1, BinaryOp2) {
	return value;
      }
    }; // inner_product<0,N_MAX>


    // Min_Element(x): result = min x_i or result = min_op x_i 
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class min_element {
    public:
      template <class Iter>
      static inline void apply_sub(Iter& first, Iter& value) {
	if (*first < *value) {
	  value = first;
	}
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(min_element,N-1,apply_sub(first, value));
      }
      template <class Iter, class Compare>
      static inline void apply_sub(Iter& first, Iter& value, Compare comp) {
	if (comp(*first , *value)) {
	  value = first;
	}
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(min_element,N-1,apply_sub(first, value, comp));
      }
    public:
      // standard comparison
      template <class Iter>
      static inline Iter apply(Iter first, Iter value) {
	SL_FASTEST_APPLY_BATCH_RECURSION(min_element,N,apply_sub(first, value));
	return value;
      }

      template <class Iter>
      static inline Iter apply(Iter first) {
	Iter value = first;
	first++;
	const size_t go = N==1?0:1;
	return ::sl::fastest::min_element<(N-1),(go*N_MAX)>::apply(first, value);
      }

      // user-defined comparison

      template <class Iter, class Compare>
      static inline Iter apply(Iter first, Iter value, Compare comp) {
	SL_FASTEST_APPLY_BATCH_RECURSION(min_element,N,apply_sub(first, value, comp));
	return value;
      }

      template <class Iter, class Compare>
      static inline Iter apply(Iter first, Compare comp) {
	Iter value = first;
	first++;
	const size_t go = N==1?0:1;
	return ::sl::fastest::min_element<N-1,(go*N_MAX)>::apply(first, value, comp);
      }
    }; // min_element<N>

    template <>
    class min_element<0,0> {
    public:
      template <class Iter>
      static inline void apply_sub(Iter&, Iter&) {
      }
      template <class Iter, class Compare>
      static inline void apply_sub(Iter&, Iter&, Compare) {
      }
    public:
      // standard comparison

      template <class Iter>
      static inline Iter apply(Iter, Iter value) {
	return value;
      }
      template <class Iter>
      static inline Iter apply(Iter first) {
	return first;
      }

      // user-defined comparison

      template <class Iter, class Compare>
      static inline Iter apply(Iter, Iter value, Compare) {
	return value;
      }
      template <class Iter, class Compare>
      static inline Iter apply(Iter first, Compare) {
	return first;
      }
    }; // min_element<0,N_MAX>

    // Max_Element(x): result = max x_i or result = max_op x_i 
    template <size_t N, size_t N_MAX = default_max_unroll_level>
    class max_element {
    public:
      template <class Iter>
      static inline void apply_sub(Iter& first, Iter& value) {
	if (*value < *first) {
	  value = first;
	}
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(max_element,N-1,apply_sub(first, value));
      }
      template <class Iter, class Compare>
      static inline void apply_sub(Iter& first, Iter& value, Compare comp) {
	if (comp(*first , *value)) {
	  value = first;
	}
	first++;
	SL_FASTEST_APPLY_TAIL_RECURSION(max_element,N-1,apply_sub(first, value, comp));
      }
    public:
      // standard comparison
      template <class Iter>
      static inline Iter apply(Iter first, Iter value) {
	SL_FASTEST_APPLY_BATCH_RECURSION(max_element,N,apply_sub(first, value));
	return value;
      }

      template <class Iter>
      static inline Iter apply(Iter first) {
	Iter value = first;
	first++;
	const size_t go = N==1?0:1;
	return ::sl::fastest::max_element<N-1,(go*N_MAX)>::apply(first, value);
      }

      // user-defined comparison

      template <class Iter, class Compare>
      static inline Iter apply(Iter first, Iter value, Compare comp) {
	SL_FASTEST_APPLY_BATCH_RECURSION(max_element,N,apply_sub(first, value, comp));
	return value;
      }

      template <class Iter, class Compare>
      static inline Iter apply(Iter first, Compare comp) {
	Iter value = first;
	first++;
	const size_t go = N==1?0:1;
	return ::sl::fastest::max_element<N-1,(go*N_MAX)>::apply(first, value, comp);
      }
    }; // max_element<N>

    template <>
    class max_element<0,0> {
    public:
      template <class Iter>
      static inline void apply_sub(Iter&, Iter&) {
      }
      template <class Iter, class Compare>
      static inline void apply_sub(Iter&, Iter&, Compare) {
      }
    public:
      // standard comparison

      template <class Iter>
      static inline Iter apply(Iter, Iter value) {
	return value;
      }
      template <class Iter>
      static inline Iter apply(Iter first) {
	return first;
      }

      // user-defined comparison

      template <class Iter, class Compare>
      static inline Iter apply(Iter, Iter value, Compare) {
	return value;
      }
      template <class Iter, class Compare>
      static inline Iter apply(Iter first, Compare) {
	return first;
      }
    }; // max_element<0,N_MAX>

#   undef SL_FASTEST_HEAD
#   undef SL_FASTEST_TAIL
#   undef SL_FASTEST_REST
#   undef SL_FASTEST_BATCH
#   undef SL_FASTEST_APPLY_TAIL_RECURSION
#   undef SL_FASTEST_APPLY_BATCH_RECURSION

  }; // namespace fastest

}; // namespace sl

#endif
