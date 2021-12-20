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

#ifndef SL_CIRCULAR_BUFFER_HPP
#define SL_CIRCULAR_BUFFER_HPP

#include <sl/cstdint.hpp>
#include <sl/serializer.hpp>
#include <sl/assert.hpp>
#include <sl/generative_types.hpp>
#include <sl/operators.hpp>

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace sl {
  
  template <typename T, typename A = std::allocator<T> >
  class circular_buffer {
  public:
    typedef circular_buffer<T,A> this_t;

    typedef T                                            value_type;
    typedef A                                            allocator_type;
    typedef typename allocator_type::size_type           size_type;
    typedef typename allocator_type::difference_type     difference_type;
    typedef typename allocator_type::reference           reference;
    typedef typename allocator_type::const_reference     const_reference;
    typedef typename allocator_type::pointer             pointer;
    typedef typename allocator_type::const_pointer       const_pointer;

  public: // iterators
    
    template <bool is_const, bool is_reverse>
    class cb_iterator {
    public:
      typedef cb_iterator<is_const, is_reverse> this_t;
    
      typedef std::random_access_iterator_tag        iterator_category; // FIXME: Currently slow
      //typedef std::bidirectional_iterator_tag      iterator_category;
      typedef T                                      value_type;
      typedef typename gen_if<is_const, const circular_buffer, circular_buffer>::type array_t;
      typedef typename gen_if<is_const, const value_type&, value_type&>::type  reference;
      typedef typename gen_if<is_const, const value_type*, value_type*>::type  pointer;
      typedef typename array_t::size_type            size_type;
      typedef typename array_t::difference_type      difference_type;

    protected:
      array_t  *cb_;
      size_type idx_;

    public:
      
      inline cb_iterator(array_t* a, size_type idx):
        cb_(a), idx_(idx) {
      }

      inline bool is_valid(const array_t* cb) const {
	return (cb_) && (cb_ == cb) && (idx_ <= cb_->size());
      }

    public: // Moving
    
      inline size_type index() const {
	return is_reverse ? cb_->size()-idx_-1 : idx_;
      }

      inline this_t& operator ++() {
	assert(cb_ && idx_<cb_->size());
	++idx_;
	return *this;
      }
      
      inline this_t& operator --() {
	assert(idx_>0);
	--idx_;
	return *this;
      }
      
      SL_OP_INCREMENTABLE(this_t);
      SL_OP_DECREMENTABLE(this_t);
      
      inline this_t& operator +=(difference_type n) {
	assert((cb_) && ((n==0) || (n>0 && idx_+n<=cb_->size()) || (n<0 && idx_+n>=0)));
	idx_ += n;
	return *this;
      }
      
      inline this_t& operator -=(difference_type n) {
	return (*this) += (-n);
      }
      
      SL_OP_ADDABLE2(this_t, difference_type);
      SL_OP_SUBTRACTABLE2(this_t, difference_type);

      inline difference_type operator -(const this_t& other) const {
	return (idx_ - other.idx_);
      }
      
    public: // deref
	  
      inline reference operator*() const {
	return (*cb_)[index()];
      }

      inline pointer operator->() const {
	return &((*cb_)[index()]);
      }
      
      inline reference operator [](difference_type n) const {
	return *(this + n);
      }
      
      inline pointer to_pointer() const {
	return (index()<cb_->size()) ? &((*cb_)[index()]) : 0;
      }

    public: // comparison
      
      inline bool operator<(const this_t& other) const {
	return idx_ < other.idx_;
      }
      
      inline bool operator==(const this_t& other) const {
	return idx_ == other.idx_;
      }
      
      SL_OP_COMPARABLE1(this_t);
      SL_OP_EQUALITY_COMPARABLE1(this_t);
    };

  public: // iterators

    typedef cb_iterator<false, false> iterator;
    typedef cb_iterator<true, false>  const_iterator;
    typedef cb_iterator<false, true>  reverse_iterator;
    typedef cb_iterator<true, true>   const_reverse_iterator;
    
  protected:

    size_type       capacity_;
    allocator_type  allocator_;
    pointer         buffer_;
    pointer         front_;
    pointer         back_; // points to next unused item

  protected: // Helpers
    
    inline pointer increment(pointer ptr, size_type by=1) {
      assert(ptr >= buffer_);
      assert(ptr <  buffer_+capacity_);
      const size_type index = ptr-buffer_;
      assert(index < capacity_);
      if (index+by >= capacity_) {
        return ptr+by-capacity_;
      } else {
        return ptr+by;
      }
    }

    inline const_pointer increment(const_pointer ptr, size_type by=1) const {
      assert(ptr >= buffer_);
      assert(ptr <  buffer_+capacity_);
      const size_type index = ptr-buffer_;
      assert(index < capacity_);
      if (index+by >= capacity_) {
        return ptr+by-capacity_;
      } else {
        return ptr+by;
      }
    }
    
    inline pointer decrement(pointer ptr, size_type by=1) {
      assert(ptr >= buffer_);
      assert(ptr <  buffer_+capacity_);
      const size_type index = ptr-buffer_;
      assert(index < capacity_);
      if (index >= by) {
        return ptr-by;
      } else {
        return ptr+(capacity_-by);
      }
    }
    
    inline const_pointer decrement(const_pointer ptr, size_type by=1) const {
      assert(ptr >= buffer_);
      assert(ptr <  buffer_+capacity_);
      const size_type index = ptr-buffer_;
      assert(index < capacity_);
      if (index >= by) {
        return ptr-by;
      } else {
        return ptr+(capacity_-by);
      }
    }
    
  public: // Create and destroy
    
    explicit circular_buffer(size_type capacity=0,
			     const allocator_type &a = allocator_type())  :
      capacity_(capacity),
      allocator_(a),
      buffer_(allocator_.allocate(capacity)),
      front_(0),
      back_(buffer_)
    {
      assert((capacity_ == 0) || (buffer_!=0));
      assert(empty());
    }
		    
    inline ~circular_buffer()  {
      clear(); // deallocates all objects
      allocator_.deallocate(buffer_, capacity_);
      buffer_ = 0;
      front_ = 0;
      back_ = 0;
    }

    allocator_type get_allocator() const {
      return allocator_;
    }

    allocator_type& get_allocator() {
      return allocator_;
    }

  public: // Copy
    
    circular_buffer(const this_t& other) :
      capacity_(other.capacity()),
      allocator_(other.get_allocator()),
      buffer_(allocator_.allocate(capacity_)),
      front_(0),
      back_(buffer_) {
      // FIXME Speed-up with direct copy
      assert((capacity_ == 0) || (buffer_!=0));
      const size_type N = other.size();
      for (std::size_t i=0; i<N; ++i) {
	push_back(other[i]);
      }
      assert(size() == other.size());
      assert(capacity() == other.capacity());
    }

    this_t &operator=(const this_t &other) {
      if (this != &other) {
  if (capacity_ != other.capacity()) {
	  allocator_.deallocate(buffer_, capacity_);
	  buffer_ = allocator_.allocate(other.capacity());
	  front_ = 0; 
	  back_ = buffer_;
	} else {
	  clear();
	}
	// FIXME Speed-up with direct copy
	const size_type N = other.size();
        for (std::size_t i=0; i<N; ++i) {
	  push_back(other[i]);
	}
	assert(size() == other.size());
	assert(capacity() == other.capacity());
      }
      return *this;
    }

  public: // swapping
    
    void swap(this_t& other) {
      std::swap(capacity_, other.capacity_);
      std::swap(allocator_, other.allocator_);
      std::swap(buffer_, other.buffer_);
      std::swap(front_, other.front_);
      std::swap(back_, other.back_);
    }
        
    friend void swap(this_t& x, this_t& y) {
      x.swap(y);
    }

  public: // serialization
      
    void store_to(output_serializer& s) const {
      s << sl::uint64_t(capacity());
      s << sl::uint64_t(size());
      for (size_type i=0; i<size(); ++i) {
	s << (*this)[i];
      }
    }
    
    void retrieve_from(input_serializer& s) {
      sl::uint64_t cap;
      sl::uint64_t sz;
      s >> cap >> sz;

      clear();
      set_capacity(cap);
      for (size_type i=0; i<sz; ++i) {
	value_type x_i; s >> x_i;
	push_back(x_i);
      }
    }

  public: // Capacity and size
    
    void set_capacity(size_type new_capacity) {
      if (new_capacity != capacity()) {

	size_type N = std::min(size(),new_capacity);
        pointer new_buffer = allocator_.allocate(new_capacity);
	
	for (size_type i=0; i<N; ++i) {
	  allocator_.construct(&(new_buffer[i]), (*this)[i]);
	}
	allocator_.deallocate(buffer_, capacity_);

	buffer_ = new_buffer;
	capacity_ = new_capacity;
	if (N) {
	  front_ = buffer_;
	  back_ = front_+N;
	} else {
	  front_ = 0; 
	  back_ = buffer_;
	}
      }
    }
    
    /// Change the size 
    void resize(size_type new_size, const_reference item = value_type()) {
      if (new_size > size()) {
	if (new_size > capacity()) {
	  set_capacity(new_size);
	}
	while (size() != new_size) {
	  push_back(item);
	}
      } else {
	while (new_size != size()) {
	  pop_back();
	}
      }
    }

    size_type size() const  {
      return !front_ ? 0
	: (back_ > front_ ? back_ : back_+capacity_) - front_;
    }
    
    size_type max_size() const  {
      return allocator_.max_size();
    }

    bool empty() const {
      return !front_;
    }

    bool full() const {
      return front_ == back_;
    }
    
    size_type capacity() const {
      return capacity_;
    }
    
  public: // Moving
    
    /// Rotate elements in the <code>circular_buffer</code>.
    void rotate(const_iterator new_begin) {
      assert(new_begin.is_valid(this)); // check for uninitialized or invalidated iterator
      assert(new_begin != end());      // check for iterator pointing to end()
      if (full()) {
	front_ = back_ = const_cast<pointer>(new_begin.to_pointer());
      } else {
	difference_type m = end() - new_begin;
	difference_type n = new_begin - begin();
	if (m < n) {
	  for (; m > 0; --m) {
	    push_front(back());
	    pop_back();
	  }
	} else {
	  for (; n > 0; --n) {
	    push_back(front());
	    pop_front();
	  }
	}
      }
    }

    /// Rotate elements in the <code>circular_buffer</code>.
    void rotate(iterator new_begin) {
      assert(new_begin.is_valid(this)); // check for uninitialized or invalidated iterator
      assert(new_begin != end());      // check for iterator pointing to end()
      if (full()) {
	front_ = back_ = const_cast<pointer>(new_begin.to_pointer());
      } else {
	difference_type m = end() - new_begin;
	difference_type n = new_begin - begin();
	if (m < n) {
	  for (; m > 0; --m) {
	    push_front(back());
	    pop_back();
	  }
	} else {
	  for (; n > 0; --n) {
	    push_back(front());
	    pop_front();
	  }
	}
      }
    }

  public: // Element access
    
    reference front()  {
      assert(front_);
      return *front_;
    }
    
    const_reference front() const  {
      assert(front_);
      return *front_;
    }
    
    reference back() {
      assert(front_);
      return *decrement(back_);
    }
    
    const_reference back() const  {
      assert(front_);
      return *decrement(back_);
    }
    
    reference operator[](size_type n)  {
      assert(front_);
      return *increment(front_,n);
    }

    const_reference operator[](size_type n) const  {
      return *increment(front_,n);
    }
    
    reference at(size_type n)  {
      if (n >= size()) throw std::out_of_range("Parameter out of range");
      return (*this)[n];
    }
    
    const_reference at(size_type n) const  {
      if (n >= size()) throw std::out_of_range("Parameter out of range");
      return (*this)[n];
    }

    iterator begin() {
      return iterator(this, 0);
    }
    
    iterator end() {
      return iterator(this, size());
    }

    const_iterator begin() const {
      return const_iterator(this, 0);
    }
    
    const_iterator end() const {
      return const_iterator(this, size());
    }

    reverse_iterator rbegin() {
      return reverse_iterator(this, 0);
    }
    
    reverse_iterator rend() {
      return reverse_iterator(this, size());
    }

    const_reverse_iterator rbegin() const {
      return const_reverse_iterator(this, 0);
    }
    
    const_reverse_iterator rend() const {
      return const_reverse_iterator(this, size());
    }
   
    void clear()  {
      if (front_) {
	do {
	  allocator_.destroy(front_);
	  front_ = increment(front_);
	} while (front_ != back_);
      }
      front_ = 0;
      back_ = buffer_;
    }

    bool push_back(const value_type& value) {
      assert(capacity()>0);
      if (full()) {
	assert(front_);
	assert(front_ == back_);
	// Overwrite
	allocator_.destroy(back_);
      }
      allocator_.construct(back_, value);

      value_type *const next = increment(back_);
      if (!front_) {
	// first entry in the buffer
	front_ = back_;
	back_ = next;
	return true;
      } else if (front_ == back_) {
	// buffer is full already, throw something away
	front_ = back_ = next;
	return false;
      } else {
	back_ = next;
	return true;
      }
    }
    
    bool push_front(const value_type& value) {
      assert(capacity()>0);
      if (empty()) {
	assert(!front_);
	allocator_.construct(back_, value);
	front_ = back_;
	back_ = increment(back_);
	return true;
      } else {
	value_type* prev = decrement(front_);
	if (front_ == back_) {
	  // buffer is full already, throw something away
	  allocator_.destroy(prev);
	  allocator_.construct(prev, value);
	  front_ = back_ = prev;
	  return false;
	} else {
	  // buffer was not full, just store
	  allocator_.construct(prev, value);
	  front_ = prev;
	  return true;
	}
      }
    }

    void pop_front()  {
      assert(!empty());

      allocator_.destroy(front_);
      value_type *const next = increment(front_);
      if (next == back_) {
	front_ = 0;
      } else {
	front_ = next;
      }
    }
    
    void pop_back()  {
      assert(!empty());
      allocator_.destroy(back_);
      value_type *const prev = decrement(back_);
      if (prev == front_) {
	front_ = 0;
      } else {
	back_ = prev;
      }
    }
    
  }; // class circular_buffer<>

} // namespace sl

#endif
