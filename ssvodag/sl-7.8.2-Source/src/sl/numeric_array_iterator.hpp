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
#ifndef SL_NUMERIC_ARRAY_ITERATOR_HPP
#define SL_NUMERIC_ARRAY_ITERATOR_HPP

#include <sl/generative_types.hpp>
#include <sl/math.hpp>

namespace sl {

  //----------------------------------------------------------
  /**
   *  Helper for defining 1D forward iterators that keep
   *  track of the current array index.
   *  Subclasses should define all features marked
   *  SUBCLASS RESPONSIBILITY.
   */
  template <class P, class I, class W> 
  class forward_1d_iterator_base {
  public:
    typedef I                            sub_t;
    typedef P                            value_t;

    typedef value_t                      value_type;
    typedef int                          difference_type;
    typedef const value_type*            pointer;    
    typedef const value_type&            reference;

  protected:

    /// User-defined data, declared here to remove empty class overhead 
    W data_;

  protected:

    /// Data init
    inline explicit forward_1d_iterator_base(const W& data): data_(data) {
    }

    /// Default init
    inline explicit forward_1d_iterator_base() {
    }

  public:
      
    /// This converted to the actual subclass type
    inline sub_t* me() {
      return static_cast<sub_t>(this);
    }

    /// This converted to the actual subclass type
    inline const sub_t* me() const {
      return static_cast<sub_t>(this);
    }

  public: // Subclass responsibility features

    /// Element access (SUBCLASS RESPONSIBILITY)
    inline value_t operator() () const {
      return me()->operator();
    }

    /// Current index (SUBCLASS RESPONSIBILITY)
    inline int index() const {
      return me()->index();
    }

    /// Pre-increment (SUBCLASS RESPONSIBILITY)
    inline sub_t& operator++() {
      return me()->operator++();
    }
    
    /// Increment by n steps (SUBCLASS RESPONSIBILITY)
    inline sub_t& operator+=(difference_type n) {
      return me()->operator+= (n);
    }
 
    /// Is other pointing to an element with an index lower than this? (SUBCLASS RESPONSIBILITY)
    inline bool operator<(const sub_t& other) const {
      return me()->operator<(other);
    }

    /// Is other pointing to the same element? (SUBCLASS RESPONSIBILITY)
    inline bool operator==(const sub_t& other) const {
      return me()->operator=(other);
    }

  public: // Features inherited without change

    /// Element access (alias of operator())
    inline value_t value() const {
      return operator();
    }

    /// Post-increment
    inline sub_t operator++(int) { 
      sub_t tmp(this); ++(*me()); return tmp; 
    }

    /// This incremented n times
    inline sub_t& operator+(difference_type n) {
      return sub_t(me()) += n;
    }
    
    /// Is other pointing to an element with an index higher than this?
    inline bool operator> (const sub_t& other) const { 
      return other < *me(); 
    }

    /// Is other pointing to an element with an index less than or equal to this?
    inline bool operator<=(const sub_t& other) const { 
      return !(other < *me()); 
    }

    /// Is other pointing to an element with an index greater than or equal to this?
    inline bool operator>=(const sub_t& other) const { 
      return !(*me() < other); 
    }

    /// Is other pointing to a different element?
    inline bool operator!=(const sub_t& other) const { 
      return !(*me() == other); 
    }
  };

  //----------------------------------------------------------
  /**
   *  Helper for defining 2D forward iterators that keep
   *  track of the current array index.
   *  Subclasses should define all features marked
   *  SUBCLASS RESPONSIBILITY.
   */
  template <class P, class I, class W> 
  class forward_2d_iterator_base {
  public:
    typedef I                            sub_t;
    typedef P                            value_t;

    typedef value_t                      value_type;
    typedef int                          difference_type;
    typedef const value_type*            pointer;    
    typedef const value_type&            reference;

  protected:

    /// User-defined data, declared here to remove empty class overhead 
    W data_;

  protected:

    /// Data init
    inline explicit forward_2d_iterator_base(const W& data): data_(data) {
    }

    /// Default init
    inline explicit forward_2d_iterator_base() {
    }

  public:
      
    /// This converted to the actual subclass type
    inline sub_t* me() {
      return static_cast<sub_t>(this);
    }

    /// This converted to the actual subclass type
    inline const sub_t* me() const {
      return static_cast<sub_t>(this);
    }

  public: // Subclass responsibility features

    /// Element access (SUBCLASS RESPONSIBILITY)
    inline value_t operator() () const {
      return me()->operator();
    }

    /// Current row index (SUBCLASS RESPONSIBILITY)
    inline int row_index() const {
      return me()->row_index();
    }

    /// Current column index (SUBCLASS RESPONSIBILITY)
    inline int column_index() const {
      return me()->column_index();
    }

    /// Pre-increment (SUBCLASS RESPONSIBILITY)
    inline sub_t& operator++() {
      return me()->operator++();
    }
    
    /// Increment by n steps (SUBCLASS RESPONSIBILITY)
    inline sub_t& operator+=(difference_type n) {
      return me()->operator+= (n);
    }
 
    /// Is other pointing to an element with an index lower than this? (SUBCLASS RESPONSIBILITY)
    inline bool operator<(const sub_t& other) const {
      return me()->operator<(other);
    }

    /// Is other pointing to the same element? (SUBCLASS RESPONSIBILITY)
    inline bool operator==(const sub_t& other) const {
      return me()->operator=(other);
    }

  public: // Features inherited without change

    /// Element access (alias of operator())
    inline value_t value() const {
      return operator();
    }

    /// Pre-increment
    inline sub_t operator++(int) { 
      sub_t tmp(this); ++(*me()); return tmp; 
    }

    /// This incremented n times
    inline sub_t& operator+(difference_type n) {
      return sub_t(me()) += n;
    }
    
    /// Is other pointing to an element with an index higher than this?
    inline bool operator> (const sub_t& other) const { 
      return other < *me(); 
    }

    /// Is other pointing to an element with an index less than or equal to this?
    inline bool operator<=(const sub_t& other) const { 
      return !(other < *me()); 
    }

    /// Is other pointing to an element with an index greater than or equal to this?
    inline bool operator>=(const sub_t& other) const { 
      return !(*me() < other); 
    }

    /// Is other pointing to a different element?
    inline bool operator!=(const sub_t& other) const { 
      return !(*me() == other); 
    }
  };



  //----------------------------------------------------------
  /**
   *  Forward const iterators on dense 1D numeric arrays.
   */
  template <class A> 
  class dense_1d_const_iterator: 
    public forward_1d_iterator_base< typename A::value_t, dense_1d_const_iterator<A>, const A * > {
  public:
    typedef A                                                  array_t;
    typedef dense_1d_const_iterator< array_t >                 self_t;
    typedef forward_1d_iterator_base< typename A::value_t, self_t, const array_t *> super_t;

    typedef typename array_t::value_t                          value_t;

    typedef typename super_t::value_type                       value_type;
    typedef typename super_t::difference_type                  difference_type;
    typedef typename super_t::pointer                          pointer;
    typedef typename super_t::reference                        reference;

  protected:

    int            i_;
    
  public:
    
    inline dense_1d_const_iterator(const array_t *ref = NULL, int i=0): super_t(ref), i_(i) {
    }
    
    inline dense_1d_const_iterator(const self_t& other): super_t(other.data_), i_(other.i_) {
    }
    
  public: // Subclass responsibility features

    /// Element access
    inline value_t operator() () const {
      SL_REQUIRE("Vector exists", data_ != NULL);
      SL_REQUIRE("Not off", (i_>=0 && i_<data_->size())); 
      return (*data_)(i_);
    }

    /// Current index
    inline int index() const {
      return i_;
    }

    /// Pre-increment
    inline self_t& operator++() {
      ++i_;
      return (*this);
    }
    
    /// Increment by n steps
    inline self_t& operator+=(difference_type n) {
      i_ += n ;
      return (*this);
    }
 
    /// Is other pointing to an element with an index lower than this?
    inline bool operator<(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ < other.i_;
    }

    /// Is other pointing to an element with an index lower than this?
    inline bool operator==(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ == other.i_;
    }

  };

  //----------------------------------------------------------
  /**
   *  Forward iterators on elements of dense 1D numeric arrays.
   */
  template <class A> 
  class dense_1d_iterator: 
    public forward_1d_iterator_base< typename A::value_t, dense_1d_iterator<A>, A *  > {
  public:
    typedef A                                                     array_t;
    typedef dense_1d_iterator< array_t >                          self_t;
    typedef forward_1d_iterator_base< typename A::value_t, self_t, array_t * > super_t;

    typedef typename array_t::value_t                             value_t;

    typedef typename super_t::value_type                          value_type;
    typedef typename super_t::difference_type                     difference_type;
    typedef typename super_t::pointer                             pointer;
    typedef typename super_t::reference                           reference;

  protected:

    int            i_;
    
  public:
    
    inline dense_1d_iterator(array_t *ref = NULL, int i=0): super_t(ref), i_(i) {
    }
    
    inline dense_1d_iterator(const self_t& other): super_t(other.data_), i_(other.i_) {
    }
    
  public: // Subclass responsibility features

    /// Element access
    inline value_t operator() () const {
      SL_REQUIRE("Vector exists", data_ != NULL);
      SL_REQUIRE("Not off", (i_>=0 && i_<data_->size())); 
      return (*data_)(i_);
    }

    /// Element write access
    inline value_t& operator() () {
      SL_REQUIRE("Vector exists", data_ != NULL);
      SL_REQUIRE("Not off", i_>=0 && i_<data_->size()); 
      return (*data_)(i_);
    }

    /// Current index
    inline int index() const {
      return i_;
    }

    /// Pre-increment
    inline self_t& operator++() {
      ++i_;
      return (*this);
    }
    
    /// Increment by n steps
    inline self_t& operator+=(difference_type n) {
      i_ += n ;
      return (*this);
    }
 
    /// Is other pointing to an element with an index lower than this?
    inline bool operator<(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ < other.i_;
    }

    /// Is other pointing to an element with an index lower than this?
    inline bool operator==(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ == other.i_;
    }

  };


  //----------------------------------------------------------
  /**
   *  Forward const iterators on non zero elements of dense 1D numeric arrays.
   */
  template <class A> 
  class dense_1d_non_zero_const_iterator: 
    public forward_1d_iterator_base< typename A::value_t, dense_1d_non_zero_const_iterator<A>, const A * > {
  public:
    typedef A                                                   array_t;
    typedef dense_1d_non_zero_const_iterator< array_t >         self_t;
    typedef forward_1d_iterator_base< typename A::value_t, self_t, const array_t *> super_t;

    typedef typename array_t::value_t                          value_t;

    typedef typename super_t::value_type                       value_type;
    typedef typename super_t::difference_type                  difference_type;
    typedef typename super_t::pointer                          pointer;
    typedef typename super_t::reference                        reference;

  protected:

    int            i_;
    
  public:
    
    inline dense_1d_non_zero_const_iterator(const array_t *ref = NULL, int i=0): super_t(ref), i_(i) {
      if (data_ && i_>=0 && i_<data_->size() && !is_zero(data_(i))) {
	this->operator++();
      }
      SL_ENSURE("empty or off or zero", ref == NULL || !(i_>=0 && i_<data_->size()) || !is_zero(operator->()));
    }
    
    inline dense_1d_non_zero_const_iterator(const self_t& other): super_t(other.data_), i_(other.i_) {
    }
    
  public: // Subclass responsibility features

    /// Element access
    inline value_t operator() () const {
      SL_REQUIRE("Vector exists", data_ != NULL);
      SL_REQUIRE("Not off", (i_>=0 && i_<data_->size())); 
      return (*data_)(i_);
    }

    /// Current index
    inline int index() const {
      return i_;
    }

    /// Pre-increment
    inline self_t& operator++() {
      SL_REQUIRE("Not off", (i_>=0 && i_<data_->size())); 
      do {
	++i_;
      } while (i_ < data_->size() && !is_zero(data_(i)));
      return (*this);
    }
    
    /// Increment by n steps
    inline self_t& operator+=(difference_type n) {
      SL_REQUIRE("Non-negative delta", i>=0);
      for (int i=0; i<n; i++) {
	++(*this);
      }	
      return (*this);
    }
 
    /// Is other pointing to an element with an index lower than this?
    inline bool operator<(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ < other.i_;
    }

    /// Is other pointing to an element with an index lower than this?
    inline bool operator==(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ == other.i_;
    }

  };

  //----------------------------------------------------------
  /**
   *  Forward iterators on non zero elements of dense 1D numeric arrays.
   */
  template <class A> 
  class dense_1d_non_zero_iterator: 
    public forward_1d_iterator_base< typename A::value_t, dense_1d_non_zero_iterator<A>, A * > {
  public:
    typedef A                                                  array_t;
    typedef dense_1d_non_zero_iterator< array_t >               self_t;
    typedef forward_1d_iterator_base< typename A::value_t, self_t, array_t *>       super_t;

    typedef typename array_t::value_t                          value_t;

    typedef typename super_t::value_type                       value_type;
    typedef typename super_t::difference_type                  difference_type;
    typedef typename super_t::pointer                          pointer;
    typedef typename super_t::reference                        reference;

  protected:

    int            i_;
    
  public:
    
    inline dense_1d_non_zero_iterator(const array_t *ref = NULL, int i=0): super_t(ref), i_(i) {
      if (data_ && i_>=0 && i_<data_->size() && !is_zero(data_(i))) {
	this->operator++();
      }
      SL_ENSURE("empty or off or zero", ref == NULL || !(i_>=0 && i_<data_->size()) || !is_zero(operator->()));
    }
    
    inline dense_1d_non_zero_iterator(const self_t& other): super_t(other.data_), i_(other.i_) {
    }
    
  public: // Subclass responsibility features

    /// Element access
    inline value_t operator() () const {
      SL_REQUIRE("Vector exists", data_ != NULL);
      SL_REQUIRE("Not off", (i_>=0 && i_<data_->size())); 
      return (*data_)(i_);
    }

    /// Element access
    inline value_t& operator() () {
      SL_REQUIRE("Vector exists", data_ != NULL);
      SL_REQUIRE("Not off", i_>=0 && i_<data_->size()); 
      return (*data_)(i_);
    }

    /// Current index
    inline int index() const {
      return i_;
    }

    /// Pre-increment
    inline self_t& operator++() {
      SL_REQUIRE("Not off", (i_>=0 && i_<data_->size())); 
      do {
	++i_;
      } while (i_ < data_->size() && !is_zero(data_(i)));
      return (*this);
    }
    
    /// Increment by n steps
    inline self_t& operator+=(difference_type n) {
      SL_REQUIRE("Non-negative delta", i>=0);
      for (int i=0; i<n; i++) {
	++(*this);
      }	
      return (*this);
    }
 
    /// Is other pointing to an element with an index lower than this?
    inline bool operator<(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ < other.i_;
    }

    /// Is other pointing to an element with the same index as this?
    inline bool operator==(const self_t& other) const {
      SL_REQUIRE("Same base", data_ == other.data_);
      return i_ == other.i_;
    }

  };

  //----------------------------------------------------------
  /**
   *  Forward const iterators that map another iterator's 
   *  values using an unary function
   */
  template <class P, class T_iter, bool T_iter_is_dense, class T_op> 
  class unary_mapped_1d_iterator: 
    public forward_1d_iterator_base<P, unary_mapped_1d_iterator<P,T_iter,T_iter_is_dense,T_op>, T_iter > {
  public:
    typedef T_iter                                                  iter_t;
    typedef T_op                                                    op_t;
    typedef unary_mapped_1d_iterator<P,T_iter,T_iter_is_dense,T_op> self_t;
    typedef forward_1d_iterator_base<P,  self_t, iter_t>               super_t;

    typedef P                                                       value_t;

    typedef typename super_t::value_type                            value_type;
    typedef typename super_t::difference_type                       difference_type;
    typedef typename super_t::pointer                               pointer;
    typedef typename super_t::reference                             reference;

    enum { is_zero_preserving = op_t::is_zero_preserving };
    enum { have_dense_subiter = T_iter_is_dense };

  protected:

    int            i_;
    
  public:
    
    inline unary_mapped_1d_iterator(const iter_t &ref, int i=0): super_t(ref), i_(i) {
    }
    
    inline unary_mapped_1d_iterator(const self_t& other): super_t(other.data_), i_(other.i_) {
    }
    
  public: // Subclass responsibility features

    /// Element access
    inline value_t operator() () const {
      typedef typename op_t::arg_t arg_t;
      return 
	(have_dense_subiter || i_ == data_.index()) ? 
	static_cast<value_t>(op_t::apply(static_cast<arg_t>(data_()))) : 
	static_cast<value_t>(op_t::apply(scalar_math<arg_t>::zero()));
    }

    /// Current index
    inline int index() const {
      return i_;
    }

    /// Pre-increment
    inline self_t& operator++() {
      if (is_zero_preserving || have_dense_subiter) {
	++data_; i_ = data_.index();
      } else if (i_ == data_.index()) {
	++i_; ++data_;
      } else {
	++i_;
      }
      return (*this);
    }
    
    /// Increment by n steps
    inline self_t& operator+=(difference_type n) {
      SL_REQUIRE("Non-negative delta", i>=0);
      if (is_zero_preserving || have_dense_subiter) {
	data_ += n; i_ = data_.index();
      } else {
	for (int i=0; i<n; i++) {
	  ++(*this);
	}	
      }
      return (*this);
    }
 
    /// Is other pointing to an element with an index lower than this?
    inline bool operator<(const self_t& other) const {
      return i_ < other.i_;
    }

    /// Is other pointing to an element with the same index as this?
    inline bool operator==(const self_t& other) const {
      return i_ == other.i_;
    }

  };

  //----------------------------------------------------------
  /**
   *  Forward const iterators that map other two iterator's 
   *  values using a binary function
   */
  template <class P, class T1_iter, class T2_iter, bool T_iters_are_dense, class T_op> 
  class binary_mapped_1d_iterator: 
    public forward_1d_iterator_base<P, binary_mapped_1d_iterator<P,T1_iter,T2_iter,T_iters_are_dense,T_op>, T1_iter > {
  public:
    typedef T1_iter                                            iter1_t;
    typedef T2_iter                                            iter2_t;
    typedef T_op                                               op_t;
    typedef binary_mapped_1d_iterator<P,T1_iter,T2_iter,T_iters_are_dense,T_op>  self_t;
    typedef forward_1d_iterator_base<P, self_t, iter1_t>       super_t;

    typedef P                                                  value_t;

    typedef typename super_t::value_type                       value_type;
    typedef typename super_t::difference_type                  difference_type;
    typedef typename super_t::pointer                          pointer;
    typedef typename super_t::reference                        reference;

    enum { is_zero_preserving = op_t::is_zero_preserving };
    enum { have_dense_subiters = T_iters_are_dense };

  protected:

    T2_iter        data2_;

    int            i_;
    
  public:
    
    inline binary_mapped_1d_iterator(const iter1_t &ref1,
				     const iter2_t &ref2,
				     int i=0): super_t(ref1), data2_(ref2), i_(i) {
    }
    
    inline binary_mapped_1d_iterator(const self_t& other): super_t(other.data_), data2_(other.data2_), i_(other.i_) {
    }
    
  public: // Subclass responsibility features

    /// Element access
    inline value_t operator() () const {
      typedef typename op_t::arg1_t arg1_t;
      typedef typename op_t::arg2_t arg2_t;
      return 
	(have_dense_subiters) ?
	(static_cast<value_t>(op_t::apply(static_cast<arg1_t>(data_()),static_cast<arg2_t>(data2_())))) : 
	((i_ == data_.index()) ? 
	 ((i_ == data2_.index()) ?
	  static_cast<value_t>(op_t::apply(static_cast<arg1_t>(data_()),static_cast<arg2_t>(data2_()))) : 
	  static_cast<value_t>(op_t::apply(static_cast<arg1_t>(data_()),scalar_math<arg2_t>::zero()))) :
	 ((i_ == data2_.index()) ?
	  static_cast<value_t>(scalar_math<arg1_t>::zero(),static_cast<arg2_t>(data2_())) : 
	  static_cast<value_t>(scalar_math<arg1_t>::zero(),scalar_math<arg2_t>::zero())));
    }

    /// Current index
    inline int index() const {
      return i_;
    }

    /// Pre-increment
    inline self_t& operator++() {
      if (have_dense_subiters) {
	++i_; ++data_; ++data2_; 
      } else if (is_zero_preserving) {
	++data_; ++data2_; 
	i_ = min(data_.index(),data2_.index());
      } else if (i_ == data_.index()) {
	if (i_ == data2_.index()) {
	  ++i_; ++data_; ++data2_;
	} else {
	  ++i_; ++data_;
	}
      } else if (i_ == data2_.index()) {
	++i_;  ++data2_;
      } else {
	++i_;
      }
      return (*this);
    }
    
    /// Increment by n steps
    inline self_t& operator+=(difference_type n) {
      SL_REQUIRE("Non-negative delta", i>=0);
      if (have_dense_subiters) {
	data_ += n; data_2 += n; i_ = data_.index();
      } else {
	for (int i=0; i<n; i++) {
	  ++(*this);
	}	
      }
      return (*this);
    }
 
    /// Is other pointing to an element with an index lower than this?
    inline bool operator<(const self_t& other) const {
      return i_ < other.i_;
    }

    /// Is other pointing to an element with the same index as this?
    inline bool operator==(const self_t& other) const {
      return i_ == other.i_;
    }

  };

}

#endif
