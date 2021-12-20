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
#ifndef SL_ITERATOR_HPP
#define SL_ITERATOR_HPP

#include <sl/assert.hpp>
#include <sl/operators.hpp>
#include <sl/generative_types.hpp>
#include <sl/type_traits.hpp>
#include <iterator>
#include <cstddef>

# if 0
namespace std {
  /**
   *  Iterator traits for restricted (unaliased) pointers.
   */
  template <class _Tp>
  struct iterator_traits<_Tp* restrict> {
    typedef random_access_iterator_tag  iterator_category;
    typedef _Tp                         value_type;
    typedef ptrdiff_t                   difference_type;
    typedef _Tp* restrict               pointer;
    typedef _Tp& restrict               reference;
  };
};
# endif

namespace sl {

  /**
   *  A wrapper class for constructing a const_iterator from an iterator
   */
  template <class Wrapped_Iterator>
  class const_random_access_iterator {
    //    SL_COMPILE_TIME_CHECK("Random access iterator",
    //			  is_same<
    //                        typename std::iterator_traits<Wrapped_Iterator>::iterator_category,
    //			    std::random_access_iterator_tag
    //                      >::value);
  public:
    typedef Wrapped_Iterator iterator_type;
    typedef const_random_access_iterator<Wrapped_Iterator> self_t;

    typedef typename std::iterator_traits<Wrapped_Iterator>::iterator_category iterator_category;
    typedef typename std::iterator_traits<Wrapped_Iterator>::value_type        base_value_type;
    typedef const base_value_type                                       value_type;

    typedef typename std::iterator_traits<Wrapped_Iterator>::difference_type   difference_type;
    typedef const value_type*                                           pointer;
    typedef const value_type&                                           reference;

  protected:
    
    iterator_type   base_iterator_;
    
  public: // Creating

    inline const_random_access_iterator(const iterator_type& i) : base_iterator_(i) {
    }
	
  public: // Moving
    
    inline self_t& operator ++() {
      base_iterator_ ++;
      return (*this);
    }
    
    inline self_t& operator --() {
      base_iterator_ --;
      return (*this);
    }
    
    SL_OP_INCREMENTABLE(self_t);
    SL_OP_DECREMENTABLE(self_t);
    
    inline self_t& operator +=(difference_type n) {
      base_iterator_ += n ;
      return (*this);
    }
      
    inline self_t& operator -=(difference_type n) {
      base_iterator_ -= n ;
      return (*this);
      }
    
    SL_OP_ADDABLE2(self_t, difference_type);
    SL_OP_SUBTRACTABLE2(self_t, difference_type);
    
    inline difference_type operator -(const self_t& other) const {
      return (base_iterator_ - other.base_iterator_);
    }

  public: // deref

    inline reference operator* () const { 
      return * base_iterator_;
    }

    inline reference operator [](difference_type n) const {
      return *(*this + n);
    }
				
  public: // comparison

    inline bool operator<(const self_t& other) const {
      return base_iterator_ < other.base_iterator_;
    }

    inline bool operator==(const self_t& other) const {
      return base_iterator_ == other.base_iterator_;
    }
    
    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

  }; // class const_random_access_iterator

}; // namespace sl

SL_OP_ADDABLE2_OVERLOADS(SL_OP_TEMPLATE1(template <class Iterator_>),
			 SL_OP_TEMPLATE1(sl::const_random_access_iterator<Iterator_>),
			 SL_OP_TEMPLATE1(typename sl::const_random_access_iterator<Iterator_>::difference_type));
SL_OP_SUBTRACTABLE2_OVERLOADS(SL_OP_TEMPLATE1(template <class Iterator_>),
			      SL_OP_TEMPLATE1(sl::const_random_access_iterator<Iterator_>),
			      SL_OP_TEMPLATE1(typename sl::const_random_access_iterator<Iterator_>::difference_type));

namespace sl {

  /** 
   *  A class for selecting a "good" const_iterator from an iterator.
   */
  template <class Wrapped_Iterator>
  class const_random_access_iterator_selector {
    //    SL_COMPILE_TIME_CHECK("Random access iterator",
    //			  is_same<
    //			  typename std::iterator_traits<Wrapped_Iterator>::iterator_category,
    //			  std::random_access_iterator_tag
    //			  >::value);
  public:
    typedef typename std::iterator_traits<Wrapped_Iterator>::value_type iter_value_type;

    typedef typename gen_if< 
                        is_const < iter_value_type >::value,
			Wrapped_Iterator,
			typename gen_if< is_same< Wrapped_Iterator, iter_value_type* >::value,
					 const iter_value_type*,
                                         const_random_access_iterator<Wrapped_Iterator> >::type >::type type;
  };
  
}; // namespace sl

namespace sl {

  /**
   *  One-dimensional strided iterators. They wrap random access iterators,
   *  jumping stride-length positions at each step.
   */
  template <class Wrapped_Iterator>
  class strided_1d_iterator {
  public:
    typedef typename std::iterator_traits<Wrapped_Iterator>::iterator_category iterator_category;
    typedef typename std::iterator_traits<Wrapped_Iterator>::value_type value_type;
    typedef typename std::iterator_traits<Wrapped_Iterator>::pointer pointer;
    typedef typename std::iterator_traits<Wrapped_Iterator>::reference reference;

    typedef Wrapped_Iterator iterator_type;
    typedef strided_1d_iterator<Wrapped_Iterator> self_t;
    typedef ptrdiff_t                             difference_type;

  protected:
    
    iterator_type   base_iterator_;
    difference_type stride_length_;
    
  public: // Creating

    inline strided_1d_iterator(const iterator_type& i, 
			       const difference_type& stride_length): base_iterator_(i), stride_length_(stride_length) {
      SL_REQUIRE("Non-null stride length", stride_length != 0);
    }
	
  public: // Moving
    
    inline self_t& operator ++() {
      base_iterator_ += stride_length_;
      return (*this);
    }
    
    inline self_t& operator --() {
      base_iterator_ -= stride_length_;
      return (*this);
    }
    
    SL_OP_INCREMENTABLE(self_t);
    SL_OP_DECREMENTABLE(self_t);
    
    inline self_t& operator +=(difference_type n) {
      base_iterator_ += n * stride_length_;
      return (*this);
    }
      
    inline self_t& operator -=(difference_type n) {
      base_iterator_ -= n * stride_length_;
      return (*this);
      }
    
    SL_OP_ADDABLE2(self_t, difference_type);
    SL_OP_SUBTRACTABLE2(self_t, difference_type);
    
    inline difference_type operator -(const self_t& other) const {
      return (base_iterator_ - other.base_iterator_)/stride_length_;
    }

  public: // deref

    inline reference operator* () const { 
      return * base_iterator_;
    }

    inline reference operator [](difference_type n) const {
      return *(*this + n);
    }
				
  public: // comparison

    inline bool operator<(const self_t& other) const {
      return base_iterator_ < other.base_iterator_;
    }

    inline bool operator==(const self_t& other) const {
      return base_iterator_ == other.base_iterator_;
    }
    
    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

  }; // strided_1d_iterator

}; // namespace sl

SL_OP_ADDABLE2_OVERLOADS(SL_OP_TEMPLATE1(template <class Iterator_>),
			 SL_OP_TEMPLATE1(sl::strided_1d_iterator<Iterator_>),
			 ptrdiff_t);
SL_OP_SUBTRACTABLE2_OVERLOADS(SL_OP_TEMPLATE1(template <class Iterator_>),
			      SL_OP_TEMPLATE1(sl::strided_1d_iterator<Iterator_>),
			      ptrdiff_t);


// --------------------------------------------------------------------
// -- fixed_strided_1d_iterator<Wrapped_Iterator>
// --------------------------------------------------------------------

namespace sl {

  /**
   *  One-dimensional strided iterators, with stride length selected
   *  at compile-time. They wrap random access iterators,
   *  jumping stride length positions at each step.
   */
  template <class Wrapped_Iterator, ptrdiff_t stride_length>
  class fixed_strided_1d_iterator {
  public:

    SL_COMPILE_TIME_CHECK("Non-null stride length", stride_length != 0);

    enum { stride_length_ = stride_length };
 
    typedef typename std::iterator_traits<Wrapped_Iterator>::iterator_category iterator_category;
    typedef typename std::iterator_traits<Wrapped_Iterator>::value_type value_type;
    typedef typename std::iterator_traits<Wrapped_Iterator>::pointer pointer;
    typedef typename std::iterator_traits<Wrapped_Iterator>::reference reference;

    typedef Wrapped_Iterator iterator_type;
    typedef fixed_strided_1d_iterator<Wrapped_Iterator, stride_length> self_t;
    typedef ptrdiff_t difference_type;

  protected:
    
    iterator_type   base_iterator_;
    
  public: // Creating

    inline fixed_strided_1d_iterator(const iterator_type& i): base_iterator_(i) {
    }
	
  public: // Moving
    
    inline self_t& operator ++() {
      base_iterator_ += stride_length_;
      return (*this);
    }
    
    inline self_t& operator --() {
      base_iterator_ -= stride_length_;
      return (*this);
    }
    
    SL_OP_INCREMENTABLE(self_t);
    SL_OP_DECREMENTABLE(self_t);
    
    inline self_t& operator +=(difference_type n) {
      base_iterator_ += n * stride_length_;
      return (*this);
    }
      
    inline self_t& operator -=(difference_type n) {
      base_iterator_ -= n * stride_length_;
      return (*this);
      }
    
    SL_OP_ADDABLE2(self_t, difference_type);
    SL_OP_SUBTRACTABLE2(self_t, difference_type);
    
    inline difference_type operator -(const self_t& other) const {
      return (base_iterator_ - other.base_iterator_)/stride_length_;
    }

  public: // deref

    inline reference operator* () const { 
      return * base_iterator_;
    }

    inline reference operator [](difference_type n) const {
      return *(*this + n);
    }
				
  public: // comparison

    inline bool operator<(const self_t& other) const {
      return base_iterator_ < other.base_iterator_;
    }

    inline bool operator==(const self_t& other) const {
      return base_iterator_ == other.base_iterator_;
    }
    
    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

  }; // fixed_strided_1d_iterator

}; // namespase sl

SL_OP_ADDABLE2_OVERLOADS(SL_OP_TEMPLATE2(template <class Iterator_, ptrdiff_t stride_length>),
			 SL_OP_TEMPLATE2(sl::fixed_strided_1d_iterator<Iterator_, stride_length>),
			 ptrdiff_t);
SL_OP_SUBTRACTABLE2_OVERLOADS(SL_OP_TEMPLATE2(template <class Iterator_, ptrdiff_t stride_length>),
			      SL_OP_TEMPLATE2(sl::fixed_strided_1d_iterator<Iterator_, stride_length>),
			      ptrdiff_t);


namespace sl {

  /**
   *  Strided one-dimensional accessors, i.e. objects that implement
   *  STL vector access interface using strided iterators.
   */
  template <class Wrapped_Iterator>
  class strided_1d_accessor {
  public: // constants & types

    typedef Wrapped_Iterator                                                      base_iterator;
    typedef strided_1d_iterator<Wrapped_Iterator>                                 iterator;
    typedef typename const_random_access_iterator_selector< iterator >::type const_iterator;

    typedef typename std::iterator_traits<iterator>::value_type            value_type;
    typedef typename std::iterator_traits<iterator>::pointer               pointer;
    typedef typename std::iterator_traits<iterator>::reference             reference;

    typedef typename std::iterator_traits<const_iterator>::value_type      const_value_type;
    typedef typename std::iterator_traits<const_iterator>::pointer         const_pointer;
    typedef typename std::iterator_traits<const_iterator>::reference       const_reference;

    typedef ptrdiff_t difference_type;

  protected: // Implementation

    iterator start_;
    difference_type max_step_count_;
    
  public: // Creation & Destruction

    inline strided_1d_accessor(const base_iterator& start,
			       const difference_type& stride_length,
			       const difference_type& max_step_count): 
      start_(start,stride_length), max_step_count_(max_step_count) {
    }

  public: // Iteration
 
    inline iterator begin() {
      return start_;
    }

    inline const_iterator begin() const {
      return const_iterator(start_);
    }

    inline iterator end() {
      return start_ + max_step_count_;
    }

    inline const_iterator end() const {
      return const_iterator(start_ + max_step_count_);
    }

    inline size_t size() const {
      return size_t(end() - begin());
    }

    inline bool empty() const {
      return size() == 0;
    }

    inline reference       operator[](size_t __n) { return *(begin() + __n); }
    inline const_reference operator[](size_t __n) const { return *(begin() + __n); }

  }; // strided_1d_accessor

}; // namespace sl

namespace sl {

  /**
   *  Fixed strided one-dimensional accessors, with stride-length and size defined
   *  at compile time. 
   */
  template <class Wrapped_Iterator, ptrdiff_t stride_length, size_t max_step_count>
  class fixed_strided_1d_accessor {
  public: // constants & types

    enum { max_step_count_ = max_step_count };

    typedef Wrapped_Iterator                                                      base_iterator;

    typedef typename gen_if<stride_length == 1,
                            Wrapped_Iterator,
                            fixed_strided_1d_iterator<Wrapped_Iterator, stride_length> >::type iterator;

    typedef typename const_random_access_iterator_selector< iterator >::type const_iterator;

    typedef typename std::iterator_traits<iterator>::value_type            value_type;
    typedef typename std::iterator_traits<iterator>::pointer               pointer;
    typedef typename std::iterator_traits<iterator>::reference             reference;

    typedef typename std::iterator_traits<const_iterator>::value_type      const_value_type;
    typedef typename std::iterator_traits<const_iterator>::pointer         const_pointer;
    typedef typename std::iterator_traits<const_iterator>::reference       const_reference;

    typedef ptrdiff_t difference_type;

  protected: // Implementation

    iterator start_;
    
  public: // Creation & Destruction

    inline fixed_strided_1d_accessor(const base_iterator& start):
      start_(start) {
    }

  public: // Iteration
 
    inline iterator begin() {
      return start_;
    }

    inline const_iterator begin() const {
      return const_iterator(start_);
    }

    inline iterator end() {
      return start_ + max_step_count_;
    }

    inline const_iterator end() const {
      return const_iterator(start_ + max_step_count_);
    }

    inline size_t size() const {
      return size_t(end() - begin());
    }

    inline bool empty() const {
      return size() == 0;
    }

    inline reference       operator[](size_t __n) { return *(begin() + __n); }
    inline const_reference operator[](size_t __n) const { return *(begin() + __n); }

  }; // fixed_strided_1d_accessor
    
}; // namespace sl

namespace sl {

  /**
   * Strided 2D iterators, with step lengths defined at compile time. Objects
   * of this class provide 2D access to a 1D random access iterator.
   */
  template 
  < class Base_Iterator,
    ptrdiff_t twod_step_length,
    ptrdiff_t oned_stride_length,
    size_t oned_max_step_count >
  class fixed_strided_2d_iterator {

    SL_COMPILE_TIME_CHECK("Non-null 2D step length", twod_step_length != 0);
    SL_COMPILE_TIME_CHECK("Non-null 1D step length", oned_stride_length != 0);
    SL_COMPILE_TIME_CHECK("Non-null 1D max step count", oned_max_step_count != 0);

  public: // constants & types

    enum { twod_step_length_ = twod_step_length };

    typedef Base_Iterator                                                                        base_iterator;

    typedef fixed_strided_2d_iterator< base_iterator, twod_step_length, oned_stride_length, oned_max_step_count > self_t;

    typedef fixed_strided_1d_accessor< base_iterator, oned_stride_length, oned_max_step_count >  oned_accessor;

    typedef typename std::iterator_traits<base_iterator>::iterator_category  iterator_category;
    typedef oned_accessor                                                    value_type;
    typedef ptrdiff_t                                                        difference_type;
    typedef value_type*                                                      pointer;
    typedef value_type&                                                      reference;

  protected: // Implementation
    
    base_iterator   base_iterator_;
    
  public: // Creating

    inline fixed_strided_2d_iterator(const base_iterator& i): base_iterator_(i) {
    }
	
  public: // Moving
    
    inline self_t& operator ++() {
      base_iterator_ += twod_step_length_;
      return (*this);
    }
    
    inline self_t& operator --() {
      base_iterator_ -= twod_step_length_;
      return (*this);
    }
    
    SL_OP_INCREMENTABLE(self_t);
    SL_OP_DECREMENTABLE(self_t);
    
    inline self_t& operator +=(difference_type n) {
      base_iterator_ += n * twod_step_length_;
      return (*this);
    }
      
    inline self_t& operator -=(difference_type n) {
      base_iterator_ -= n * twod_step_length_;
      return (*this);
      }
    
    SL_OP_ADDABLE2(self_t, difference_type);
    SL_OP_SUBTRACTABLE2(self_t, difference_type);
    
    inline difference_type operator -(const self_t& other) const {
      return (base_iterator_ - other.base_iterator_)/twod_step_length_;
    }

  public: // deref

    inline value_type operator* () const { 
      return value_type(base_iterator_);
    }

    inline value_type operator [](difference_type n) const {
      return *(*this + n);
    }
				
  public: // comparison

    inline bool operator<(const self_t& other) const {
      return base_iterator_ < other.base_iterator_;
    }

    inline bool operator==(const self_t& other) const {
      return base_iterator_ == other.base_iterator_;
    }
    
    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

  }; // fixed_strided_2d_iterator

}; // namespace sl

SL_OP_ADDABLE2_OVERLOADS(SL_OP_TEMPLATE4(template <class Base_Iterator,
					 ptrdiff_t twod_step_length,
					 ptrdiff_t oned_stride_length,
					 size_t oned_max_step_count>),
			 SL_OP_TEMPLATE4(sl::fixed_strided_2d_iterator<Base_Iterator, 
					 twod_step_length, 
					 oned_stride_length, 
					 oned_max_step_count>),
			 ptrdiff_t); 

SL_OP_SUBTRACTABLE2_OVERLOADS(SL_OP_TEMPLATE4(template <class Base_Iterator,
					 ptrdiff_t twod_step_length,
					 ptrdiff_t oned_stride_length,
					 size_t oned_max_step_count>),
			      SL_OP_TEMPLATE4(sl::fixed_strided_2d_iterator<Base_Iterator, 
					      twod_step_length, 
					      oned_stride_length, 
					      oned_max_step_count>),
			      ptrdiff_t);

namespace sl {

  /**
   *  Two-dimensional accessors, i.e. containers accessed by 2D iterators.
   */
  template <class Base_Iterator,
    ptrdiff_t twod_step_length,
    size_t twod_max_step_count,
    ptrdiff_t oned_stride_length,
    size_t oned_max_step_count >
  class fixed_strided_2d_accessor {

    SL_COMPILE_TIME_CHECK("Non-null 2D step length", twod_step_length != 0);
    SL_COMPILE_TIME_CHECK("Non-null 2D max step count", twod_max_step_count != 0);
    SL_COMPILE_TIME_CHECK("Non-null 1D step length", oned_stride_length != 0);
    SL_COMPILE_TIME_CHECK("Non-null 1D max step count", oned_max_step_count != 0);

  public: // constants & types

    enum { twod_max_step_count_ = twod_max_step_count };

    typedef Base_Iterator                                                  base_iterator;
    typedef fixed_strided_2d_iterator<Base_Iterator, 
                                      twod_step_length, 
                                      oned_stride_length, 
                                      oned_max_step_count>                 iterator;
    typedef typename const_random_access_iterator_selector< iterator >::type const_iterator;

    typedef typename std::iterator_traits<iterator>::value_type            value_type;
    typedef typename std::iterator_traits<iterator>::pointer               pointer;
    typedef typename std::iterator_traits<iterator>::reference             reference;

    typedef typename std::iterator_traits<const_iterator>::value_type      const_value_type;
    typedef typename std::iterator_traits<const_iterator>::pointer         const_pointer;
    typedef typename std::iterator_traits<const_iterator>::reference       const_reference;

    typedef ptrdiff_t difference_type;

  protected: // Implementation

    iterator start_;
    
  public: // Creation & Destruction

    inline fixed_strided_2d_accessor(const base_iterator& start):
      start_(start) {
    }

  public: // Iteration
 
    inline iterator begin() {
      return start_;
    }

    inline const_iterator begin() const {
      return const_iterator(start_);
    }

    inline iterator end() {
      return start_ + twod_max_step_count_;
    }

    inline const_iterator end() const {
      return const_iterator(start_ + twod_max_step_count_);
    }

    inline size_t size() const {
      return size_t(end() - begin());
    }

    inline bool empty() const {
      return size() == 0;
    }

    inline value_type       operator[](size_t __n) { return *(begin() + __n); }
    inline const value_type operator[](size_t __n) const { return *(begin() + __n); }

  }; // fixed_strided_2d_accessor
    
}; // namespace sl

#endif 

