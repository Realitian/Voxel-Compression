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
#ifndef SL_ADAPTIVE_PACKED_MEMORY_ARRAY_HPP
#define SL_ADAPTIVE_PACKED_MEMORY_ARRAY_HPP

#include <sl/cstdint.hpp>
#include <sl/serializer.hpp>
#include <sl/assert.hpp>
#include <sl/generative_types.hpp>
#include <sl/operators.hpp>
#include <sl/bitops.hpp>
#include <sl/circular_buffer.hpp>
#include <vector>
#include <deque>
#include <utility>
#include <functional>
#include <cmath>
#include <cassert>

namespace sl {

  /**
   *  Adaptive PMA (packed-memory array)
   *
   *  A cache-oblivious solution for maintaining a dynamic set of elements in
   *  sorted order. It supports operations insert, delete, and scan. Because
   *  the elements are stored physically in sorted order in memory or on disk,
   *  the PMA can be used to support extremly efficient range queries.
   *
   *  The PMA is divided into $/Theta(N/logn)$ segments, and a contiguous run of
   *  segments is called a window. The PMA is viewed in terms of a tree structure,
   *  where nodes of the tree are windows. The root node is the window containing
   *  all segments and a leaf node is a window containing a single segment. The
   *  tree is implicitly rather than explicitly maintained.
   *
   *  Based on:
   *
   *    Michael A. Bender and Haodong Hu. 2006. An adaptive
   *    packed-memory array. In Proc. PODS 2006.
   *    DOI=10.1145/1142351.1142355
   *
   *  Memory complexity is O(n), search is O(log n), insert and
   *  delete are O(log^2 n). The structure is a natural choice for
   *  cache-oblivious algorithms, since consecutive elements are
   *  stored nearby.
   */
  template <class T,
	    bool is_not_allowing_duplicates = false,
	    class Compare = std::less<T>,
	    class Allocator = std::allocator<T> >
  class adaptive_packed_memory_array {
  public:
    typedef adaptive_packed_memory_array<T,is_not_allowing_duplicates,Compare,T> this_t;
    typedef std::vector<T,Allocator>              value_vector_t;
    typedef std::vector<bool>                     bool_vector_t;
    
    typedef typename value_vector_t::allocator_type     allocator_type;
    typedef typename value_vector_t::size_type          size_type;
    typedef typename value_vector_t::difference_type	difference_type;
    typedef typename value_vector_t::reference		reference;
    typedef typename value_vector_t::const_reference	const_reference;
    typedef typename value_vector_t::value_type		value_type;
    typedef typename value_vector_t::value_type		key_type;
    typedef Compare                                     key_compare;
    typedef Compare                                     value_compare;

    typedef std::pair<value_type,size_type>             predictor_pair_t;
    typedef sl::circular_buffer<predictor_pair_t>       predictor_buffer_t;

  public: // iterators
    
    template <bool is_const, bool is_reverse>
    class apma_iterator {
    public:
      typedef apma_iterator<is_const, is_reverse> this_t;
    
      // typedef std::random_access_iterator_tag        iterator_category; // FIXME: Currently slow
      typedef std::bidirectional_iterator_tag        iterator_category;
      typedef T                                      value_type;
      typedef T*                                     pointer; // FIXME
      typedef typename gen_if<is_const, const adaptive_packed_memory_array, adaptive_packed_memory_array >::type array_t;
      typedef typename gen_if<is_const, const value_type&, value_type&>::type  reference;
      typedef typename array_t::size_type            size_type;
      typedef typename array_t::difference_type      difference_type;

    protected:
      array_t  *apma_;
      size_type idx_;

    public:
      
      inline apma_iterator(array_t* a, size_type idx):
        apma_(a), idx_(idx) {
      }

    public: // Moving
    
      inline size_type index() const {
	return is_reverse ? apma_->storage_size()-idx_-1 : idx_;
      }

      inline this_t& operator ++() {
	const size_type M = apma_->storage_size();
	if (idx_<M) ++idx_;
	while ((idx_<M) && !apma_->storage_is_present(is_reverse ? M-idx_-1 : idx_)) {
	  ++idx_;
	}
	return *this;
      }
      
      inline this_t& operator --() {
	const size_type M = apma_->storage_size();
	if (idx_ > 0) --idx_;
	while ((idx_ > 0) && !apma_->storage_is_present(is_reverse ? M-idx_-1 : idx_)) {
	  --idx_;
	}
	return *this;
      }
      
      SL_OP_INCREMENTABLE(this_t);
      SL_OP_DECREMENTABLE(this_t);
      
      inline this_t& operator +=(difference_type n) {
	const size_type M = apma_->storage_size();
	if (n>0) {
	  while (n!=0 && idx_<M) {
	    if (idx_<M) ++idx_;
	    while ((idx_<M) && !apma_->storage_is_present(is_reverse ? M-idx_-1 : idx_)) {
	      ++idx_;
	    }
	    --n;
	  }
	} else if (n<0) {
	  while (n!=0 && idx_>0) {
	    if (idx_>0) ++idx_;
	    while ((idx_>0) && !apma_->storage_is_present(is_reverse ? M-idx_-1 : idx_)) {
	      --idx_;
	    }
	    ++n;
	  }
	}
	return *this;
      }
      
      inline this_t& operator -=(difference_type n) {
	return (*this) += (-n);
      }
      
      SL_OP_ADDABLE2(this_t, difference_type);
      SL_OP_SUBTRACTABLE2(this_t, difference_type);

#if 0
      // -- NO -- we should count steps...
      inline difference_type operator -(const this_t& other) const {
	return (idx_ - other.idx_);
      }
#endif       
      
    public: // deref
	  
      inline reference operator*() const {
	return apma_->storage_at(this->index());
      }
      
      inline reference operator [](difference_type n) const {
	return *(this + n);
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

    typedef apma_iterator<false, false> iterator;
    typedef apma_iterator<true, false>  const_iterator;
    typedef apma_iterator<false, true>  reverse_iterator;
    typedef apma_iterator<true, true>   const_reverse_iterator;

  protected:
    key_compare comp_;
    
    size_type implicit_tree_root_height_; // Height of root -- level count = height +1
    size_type segment_capacity_;          // Must be power of two
    size_type capacity_;                  // 2^height * segment_capacity_
    
    value_vector_t storage_;              // The actual array
    bool_vector_t  exists_;               // A bitmask indicating occupied slots
    predictor_buffer_t predictor_;        // Approximate counter of inserts falling in each segment 
    
    size_type element_count_;
    size_type first_index_;
    size_type last_index_;
    
  protected: // Helpers

    static inline bool is_pow2(size_type x) {
      return sl::bitops<size_type>::is_power2(x);
    }

    static inline size_type log2(size_type x) {
      return sl::bitops<size_type>::log2(x);
    }

    static inline size_type rounded_up_to_pow2(size_type x) {
      return sl::bitops<size_type>::next_power2(x);
    }
        
    inline bool lt(const value_type& x, const value_type&y) const {
      return comp_(x,y);
    }

    inline bool gt(const value_type& x, const value_type&y) const {
      return lt(y,x);
    }
      
    inline bool leq(const value_type& x, const value_type&y) const {
      return !gt(x,y);
    }

    inline bool geq(const value_type& x, const value_type&y) const {
      return !lt(x,y);
    }

    inline bool eq(const value_type& x, const value_type&y) const {
      return !(lt(x,y) || gt(x,y));
    }
    
  protected: // Parameters
    
    // THRESHOLDS
    //
    // The nodes at each height h have an upper density threshold t_h and a
    // lower density threshold p_h, which together determine the acceptable
    // density of keys within a window of 2^h segments. As node height
    // increases the udts decrease and ldts increase.
    // D_min = p_0 <...< p_h < t_h <...< t_0 = D_max
    // We must have 2*p_h < t_h because when we double/halve the array size, it
    // must be within bounds
    
    static inline size_type APMA_INITIAL_CAPACITY() { return 4; } // Should be a power of two
    static inline size_type APMA_SCALE_FACTOR() { return 2; }

    static inline double APMA_LEAF_LOWER_DENSITY() { return 0.20; }
    static inline double APMA_ROOT_LOWER_DENSITY() { return 0.30; } 
    
    static inline double APMA_LEAF_UPPER_DENSITY() { return 1.00; } 
    static inline double APMA_ROOT_UPPER_DENSITY() { return 0.75; } 

    inline double apma_lower_density_threshold(size_type h) const {
      double t = double(h)/double(implicit_tree_root_height_); // 1 at root, 0 at leaf
      return APMA_LEAF_LOWER_DENSITY() * (1.0-t) + APMA_ROOT_LOWER_DENSITY() * (t);
    }
    
    inline double apma_upper_density_threshold(size_type h) const {
      double t = double(h)/double(implicit_tree_root_height_); // 1 at root, 0 at leaf
      return APMA_LEAF_UPPER_DENSITY() * (1.0-t) + APMA_ROOT_UPPER_DENSITY() * (t);
    }

    inline size_type apma_window_capacity(size_type h) const {
      return segment_capacity_ << h;
    }

    inline bool apma_is_window_too_full(size_type h, size_type window_element_count) const {
      return window_element_count > apma_upper_density_threshold(h) * apma_window_capacity(h);
    }

    inline bool apma_is_too_full(size_type insert_count=0) const {
      const size_type N = element_count_ + insert_count;
      const size_type M = capacity_;
      
      return N > APMA_ROOT_UPPER_DENSITY() * M;
    }
    
    inline bool apma_is_window_too_empty(size_type h, size_type window_element_count) const {
      return window_element_count < apma_lower_density_threshold(h) * apma_window_capacity(h);
    }

    inline bool apma_is_too_empty(size_type erase_count=0) const {
      const size_type N = element_count_ - erase_count;
      const size_type M = capacity_;
      
      return N < APMA_ROOT_LOWER_DENSITY() * M;
    }

  protected: // Low-level implementation: init
         
    void apma_clear() {
      capacity_ = 0;
      implicit_tree_root_height_ = 0;
      segment_capacity_ = 0;
      
      storage_.resize(0);
      exists_.resize(0);
      predictor_.set_capacity(0);
      element_count_ = 0;
      first_index_ = 0;
      last_index_ = 0;
    }

    void apma_init_tree_shape(size_type c_desired) {
      size_type c = rounded_up_to_pow2(c_desired);
      size_type log2c = log2(c);
      segment_capacity_ = rounded_up_to_pow2(log2c);
      implicit_tree_root_height_ = log2c - log2(segment_capacity_);
      capacity_ = (segment_capacity_ << implicit_tree_root_height_);

      assert(capacity_ >= c_desired);
      assert(capacity_ == apma_window_capacity(implicit_tree_root_height_));
    }
    
    size_type apma_init(const value_type& x0) {
      assert(storage_.size() == 0);
      assert(exists_.size() == 0);
      assert(predictor_.size() == 0);

      apma_init_tree_shape(APMA_INITIAL_CAPACITY());
      
      storage_.resize(capacity_);
      exists_.resize(capacity_, false);

      apma_predictor_reset();
      
      storage_[0] = x0;
      exists_[0] = true;
      
      element_count_ = 1;
      first_index_ = 0;
      last_index_ = 0;

      assert(capacity_ == storage_.size());
      
      return first_index_;
    }

  protected: // Low-level implementation: packing

    size_type apma_packed_lower_bound(const value_type& x,
				      size_type l,
				      size_type u) const {
      assert(exists_[l]);
      assert(exists_[u]);
      assert(l<=u);
      if (lt(storage_[u],x)) {
	// x is just after u
	l = u+1; 
      } else if (lt(storage_[l], x)) {
	while (l+8<u) {
	  assert(exists_[l]);
	  assert(exists_[u]);
	    
	  size_type m = (l+u)>>1;
	  assert(exists_[m]);
	  if (lt(storage_[m],x)) {
	    // m<x -- search interval m+1..u
	    l=m+1; 
	  } else {
	    // m>=x -- search interval l..m
	    u=m;
	  }
	}
	while (l<u && lt(storage_[l], x)) {
	  ++l;
	}
      }
      assert(exists_[l]);
      assert(!lt(storage_[l],x));
      
      return l;
    }

    void apma_pack_right(size_type* packed_count,
			 size_type* packed_insert_index,
			 size_type window_start,
			 size_type window_size,
			 size_type insert_index = size_type(-1),
			 const value_type* x_insert = 0) {
      // Move all elements in window to the upper side
      *packed_count = 0;
      *packed_insert_index = size_type(-1);

      size_type window_end = window_start+window_size;
      if (!x_insert) {
	// No element to insert, just move
	difference_type next_storing_position = window_end - 1;
	for (difference_type i = next_storing_position; i >= difference_type(window_start); --i) {
	  if (exists_[i]) {
	    assert(next_storing_position>=difference_type(window_start));
	    if (size_type(i)!=size_type(next_storing_position)) {
	      storage_[next_storing_position] = storage_[i];
	      exists_[next_storing_position] = true;
	      exists_[i] = false; 
	    }
	    --next_storing_position;
	    ++(*packed_count);
	    assert(*packed_count<=window_size);
	  }
	}

	assert(*packed_count>0);
	assert(*packed_count<=element_count_);
	assert((window_size!=capacity_) || ((*packed_count) == element_count_));
      } else {
	assert(insert_index!= size_type(-1));
	assert(insert_index>=window_start);
	assert(insert_index<=window_end);
	
	// Move all elements after insert position
	size_type next_storing_position = window_end-1;
	for (difference_type i = window_end-1; i >= difference_type(insert_index); --i) {
	  if (exists_[i]) {
	    assert(difference_type(next_storing_position)>=difference_type(window_start));
	    if (size_type(i)!=next_storing_position) {
	      storage_[next_storing_position] = storage_[i];
	      exists_[next_storing_position] = true;
	      exists_[i] = false; 
	    }
	    --next_storing_position;
	    ++(*packed_count);
	    assert(*packed_count<=window_size);
	  }
	}
	// Store inserted element, eventually moving down stored ones
	*packed_insert_index = next_storing_position;
	if (!exists_[next_storing_position]) {
	  // Just store
	  assert(difference_type(next_storing_position)>=difference_type(window_start));
	  storage_[next_storing_position] = *x_insert;
	  exists_[next_storing_position] = true;
	  --next_storing_position;
	  ++(*packed_count);
	  assert(*packed_count<=window_size);
	} else {
	  // Move down!
	  value_type x_prev = *x_insert;
	  while (exists_[next_storing_position]) {
	    assert(difference_type(next_storing_position)>=difference_type(window_start));
	    std::swap(storage_[next_storing_position], x_prev);
	    --next_storing_position;
	    ++(*packed_count);
	    assert(*packed_count<=window_size);
	  }
	  // Stored last kicked-out item
	  assert(!exists_[next_storing_position]);
	  assert(difference_type(next_storing_position)>=difference_type(window_start));
	  storage_[next_storing_position] = x_prev;
	  exists_[next_storing_position] = true;
	  --next_storing_position;
	  ++(*packed_count);
	  assert(*packed_count<=window_size);
	}
	// Store elements before inserted element
	for (difference_type i = next_storing_position; i >= difference_type(window_start); --i) {
	  if (exists_[i]) {
	    assert(difference_type(next_storing_position)>=difference_type(window_start));
	    if (size_type(i)!=next_storing_position) {
	      storage_[next_storing_position] = storage_[i];
	      exists_[next_storing_position] = true;
	      exists_[i] = false; 
	    }
	    --next_storing_position;
	    ++(*packed_count);
	    assert(*packed_count<=window_size);
	  }
	}

	++element_count_;

	assert(*packed_count>0);
	assert(*packed_count<=element_count_);
	assert((window_size!=capacity_) || ((*packed_count) == element_count_));
      }

      // Update bounds
      if (first_index_ >= window_start && first_index_<window_end-*packed_count) {
	first_index_ = window_end-*packed_count;
      }
      if (last_index_ < window_end) {
	last_index_ = window_end-1;
      }

      assert(insert_index == size_type(-1) || (*packed_insert_index<capacity_));
      assert(insert_index == size_type(-1) || (exists_[*packed_insert_index]));
      assert(insert_index == size_type(-1) || (eq(*x_insert, storage_[*packed_insert_index])));
      
      assert(!*packed_count || exists_[first_index_]);
      assert(!*packed_count || exists_[last_index_]);
    }

  protected: // Low-level implementation: even redistribution

    void apma_even_redistribute_packed(size_type* output_insert_index,
				       size_type window_start,
				       size_type window_size,
				       size_type packed_start,
				       size_type packed_count,
				       size_type packed_insert_index) {
      // Now move at the right place, including elements to be inserted
      assert(packed_count<=window_size);

      bool is_moving_inserted = ((packed_insert_index>=packed_start) && (packed_insert_index<(packed_start+packed_count)));
#ifndef NDEBUG
      value_type x_old = value_type();
      if (is_moving_inserted) {
	assert(packed_insert_index<capacity_);
	assert(exists_[packed_insert_index]);
	x_old = storage_[packed_insert_index];
      }
#endif
      *output_insert_index = (is_moving_inserted?packed_insert_index:size_type(-1));

      if ((packed_count == 0) || ((packed_start==window_start) && (packed_count==window_size))) {
	// Leave as is
      } else if (packed_count==window_size) {
	// Move completely full
	size_type out_index = window_start;
	size_type in_index = packed_start;
	for(size_type i=0; i<packed_count; ++i) {
	  assert(out_index<capacity_);
	  assert(out_index>=window_start);
	  assert(out_index<window_start+window_size);
	  assert(in_index<capacity_);
	  assert(in_index>=packed_start);
	  assert(in_index<packed_start+packed_count);
	  assert(out_index<in_index);
	  storage_[out_index] = storage_[in_index];
	  exists_[out_index] = true;
	  exists_[in_index] = false;
	  if (in_index == packed_insert_index) {
	    *output_insert_index = out_index;
	  }
	  ++out_index;
	  ++in_index;
	}	  
	// Update bounds
	if (first_index_>=window_start) first_index_ = window_start;
	if (last_index_<packed_start+packed_count) last_index_ = window_start + window_size - 1;

	assert(exists_[first_index_]);
	assert(exists_[last_index_]);
      } else {
	// Copy evenly spaced
	double deltak = double(window_size)/double(packed_count);
	size_type out_count = 0;
	size_type in_index  = packed_start;
	
	// Copy elements before inserted one
	for(size_type i=0; i<packed_count; ++i) {
	  size_type out_index = window_start + size_type(double(out_count) * deltak);
	
	  assert(out_index<capacity_);
	  assert(out_index>=window_start);
	  assert(out_index<window_start+window_size);
	  assert(in_index<capacity_);
	  assert(in_index>=packed_start);
	  assert(in_index<packed_start+packed_count);
	  assert(out_index<=in_index);

	  if (out_index != in_index) {
	    storage_[out_index] = storage_[in_index];
	    exists_[out_index] = true;
	    exists_[in_index] = false;
	    if (in_index == packed_insert_index) {
	      *output_insert_index = out_index;
	    }
	  }
	  ++out_count;
	  ++in_index;
	}

	// Update bounds
	if (first_index_>=window_start) first_index_ = window_start;
	if (last_index_<packed_start+packed_count) last_index_ = window_start + size_type(double(out_count-1)*deltak);

	assert(exists_[first_index_]);
	assert(exists_[last_index_]);
      }

#ifndef NDEBUG
      if (is_moving_inserted)  {
	assert(*output_insert_index<capacity_);
	assert(*output_insert_index<=packed_insert_index);
	assert(exists_[*output_insert_index]);
	assert(eq(x_old, storage_[*output_insert_index]));
      }
#endif

    }
    
  protected: // Uneven redistribution
    
    static inline size_type apma_best_left_count(size_type min_left_count,
					    size_type max_left_count,
					    size_type child_capacity,
					    size_type x_count,
					    size_type inserts_left,
					    size_type inserts_right) {
      double il = double(inserts_left);
      double ir = double(inserts_right);
      double c = double(child_capacity);
      double n = double(x_count);
      
      return sl::median(min_left_count,
			max_left_count,
			size_type(std::max(0.0, (ir*c + il*(n-c))/(il+ir))));
    }

    void apma_uneven_redistribute_packed(size_type* output_insert_index,
					 size_type window_start,
					 size_type window_size,
					 size_type packed_start,
					 size_type packed_count,
					 const std::vector< std::pair<size_type,size_type> >& active_predictors,
					 size_type packed_insert_index) {
      size_type window_segment_count = window_size/segment_capacity_;
      size_type window_level = log2(window_segment_count);
      size_type window_capacity = segment_capacity_ * window_segment_count; 
      size_type active_predictor_count = active_predictors.size();
      
      assert(packed_count<=window_capacity);
      assert(window_segment_count>0);

      if ((window_segment_count<2) ||
	  (active_predictor_count==0) ||
	  (packed_count==window_capacity) ||
	  (packed_count==0)) {

	// Reached bottom -- evenly assign all elements to this window
	apma_even_redistribute_packed(output_insert_index,
				      window_start, window_size,
				      packed_start, packed_count,
				      packed_insert_index);
      } else {
	// Unevenly redistribute among two childrens
	size_type child_capacity = window_capacity/2;
	size_type left_count_min = std::min(size_type(std::max(apma_lower_density_threshold(window_level)*child_capacity,
							       double(packed_count)-apma_upper_density_threshold(window_level)*child_capacity)),
					    packed_count-1);
	size_type left_count_max = std::min(size_type(std::min(apma_upper_density_threshold(window_level)*child_capacity,
							       double(packed_count)-apma_lower_density_threshold(window_level)*child_capacity)),
					    packed_count-1);
	assert(left_count_min<=left_count_max);
	assert(left_count_min<packed_count);
	assert(left_count_max<packed_count);

	// Find constant regions
	std::vector< std::pair<size_type,size_type> > inserts_left;
	inserts_left.reserve(active_predictor_count+2);

	size_type left_count = left_count_min;
	size_type total_inserts=0;
	{
	  size_type pidx=0;
	  while (left_count<=left_count_max+1) {
	    while (pidx<active_predictor_count && active_predictors[pidx].first<left_count) {
	      total_inserts += active_predictors[pidx].second;
	      ++pidx;
	    }
	    inserts_left.push_back(std::make_pair(left_count,total_inserts));
	    if (pidx<active_predictor_count) {
	      // Split after next predictor
	      left_count=active_predictors[pidx].first+1; 
	    } else {
	      left_count=left_count_max+2; // Exit
	    }
	  }
	  if (inserts_left.back().first != left_count_max+1) {
	    // Insert last split
	    left_count = left_count_max+1;
	    while (pidx<active_predictor_count && active_predictors[pidx].first<left_count) {
	      total_inserts += active_predictors[pidx].second;
	      ++pidx;
	    }
	    inserts_left.push_back(std::make_pair(left_count,total_inserts));
	  }
	  // Count total inserts in this region
	  while (pidx<active_predictor_count) {
	    total_inserts += active_predictors[pidx].second;
	    ++pidx;
	  }
	}
	
	const size_type Ni = inserts_left.size();

	// Choose best split
	size_type best_left_count=left_count_min;
	float     opt_value = 1e30f;
	for (size_type i=0; i<Ni-1; ++i) {
	  size_type min_left_count_i = inserts_left[i].first;
	  size_type max_left_count_i = inserts_left[i+1].first-1;
	  size_type inserts_left_i   = inserts_left[i].second;
	  size_type inserts_right_i  = total_inserts-inserts_left_i;
	  
	  size_type best_left_count_i = apma_best_left_count(min_left_count_i, max_left_count_i,
							     child_capacity, packed_count,
							     inserts_left_i, inserts_right_i);

	  size_type gaps_left_i  = child_capacity-best_left_count;
	  size_type gaps_right_i = child_capacity-(packed_count-best_left_count);

	  float opt_value_i= std::abs(inserts_left_i/(gaps_left_i+1e-6f) -
				      inserts_right_i/(gaps_right_i+1e-6f));
	  if ((i==0) || (opt_value_i < opt_value)) {
	    best_left_count = best_left_count_i;
	    opt_value = opt_value_i;
	  }
	} // for each insert

	// Splitting at "best_left_count" -- redistribute predictors 

	std::vector< std::pair<size_type,size_type> > left_active_predictors;
	left_active_predictors.reserve(active_predictor_count);
	std::vector< std::pair<size_type,size_type> > right_active_predictors;
	right_active_predictors.reserve(active_predictor_count);
	
	{
	  size_type pidx=0;
	  while (pidx<active_predictor_count && active_predictors[pidx].first<best_left_count) {
	    left_active_predictors.push_back(active_predictors[pidx]);
	    ++pidx;
	  }
	  while (pidx<active_predictor_count) {
	    right_active_predictors.push_back(std::make_pair(active_predictors[pidx].first-best_left_count, active_predictors[pidx].second));
	    ++pidx;
	  }
	}

	// Recurse
	size_type output_insert_index_left;
	apma_uneven_redistribute_packed(&output_insert_index_left,
					window_start, child_capacity,
					packed_start, best_left_count,
					left_active_predictors,
					packed_insert_index);
	size_type output_insert_index_right;
	apma_uneven_redistribute_packed(&output_insert_index_right,
					window_start+child_capacity, child_capacity,
					packed_start+best_left_count, packed_count-best_left_count,
					right_active_predictors,
					(output_insert_index_left == size_type(-1)) ? packed_insert_index : size_type(-1));
	*output_insert_index = (output_insert_index_left == size_type(-1)) ? output_insert_index_right : output_insert_index_left;
      }
    }

    void apma_uneven_redistribute_packed(size_type* output_insert_index,
					 size_type window_start,
					 size_type window_size,
					 size_type packed_start,
					 size_type packed_count,
					 size_type packed_insert_index) {
      assert(packed_count<=window_size);

      *output_insert_index = packed_insert_index;

      if (packed_count == 0) {
	// Leave as is
      } else if ((window_size <= segment_capacity_) || (packed_count < 8) || (packed_count==window_size)) {
	// Too few elements -- not worth the effort to redistribute unevenly
	apma_even_redistribute_packed(output_insert_index,
				      window_start,
				      window_size,
				      packed_start,
				      packed_count,
				      packed_insert_index);
      } else { // Something to move with gaps
	assert(packed_count>0);

	size_type window_end = window_start+window_size;
	size_type packed_start = window_end-packed_count;

	// Extract and sort predictors within window range
	const value_type& x_lo = storage_[packed_start];
	const value_type& x_hi = storage_[window_end-1];
	
	std::vector< std::pair<size_type,size_type> >  active_predictors; // Relative to packed_start
	active_predictors.reserve(segment_capacity_);
	for (typename predictor_buffer_t::iterator p_it = predictor_.begin();
	     p_it != predictor_.end();
	     ++p_it) {
	  const predictor_pair_t& pp = *p_it;
	  if (leq(x_lo,pp.first) && leq(pp.first,x_hi)) {
	    // Within range -- search location and insert into sorted list
	    size_type location_k = std::min(window_end-1,
					    apma_packed_lower_bound(pp.first,
								    packed_start,
								    window_end-1));

	    size_type inserts_k = pp.second;
	    active_predictors.push_back(std::make_pair(location_k-packed_start, // Make it relative to region start
						       inserts_k));
	    size_type k=active_predictors.size()-1;
	    while (k>0 && (active_predictors[k]<active_predictors[k-1])) {
	      std::swap(active_predictors[k], active_predictors[k-1]);
	      --k;
	    }
	  }
	}

	const size_type Np = active_predictors.size();
	if (Np==0) {
	  // No predictors in this range -- redistribute evenly!
	  apma_even_redistribute_packed(output_insert_index,
					window_start,
					window_size,
					packed_start,
					packed_count,
					packed_insert_index);
	} else {
	  // Predictors -> uneven rebalancing
	  apma_uneven_redistribute_packed(output_insert_index,
					  window_start, window_size,
					  packed_start, packed_count,
					  active_predictors,
					  packed_insert_index);
	} // if predictors
      } // if totally full or totally empty
    }
    
  protected: // Low-level implementation: predictors

    // Predictor maintains the list of recent elements that have been displaced
    // by inserts, with approximate insert counts

    void apma_predictor_reset() {
      predictor_.clear();
      predictor_.set_capacity(segment_capacity_);
    }

    void apma_predictor_update_on_insert(size_type insert_index) {
      if (insert_index>last_index_) insert_index = last_index_;
      if (insert_index<first_index_) insert_index = first_index_;

      if (insert_index<capacity_ && exists_[insert_index]) {
	const value_type& x = storage_[insert_index];
	
	// Look for x in predictor buffer -- Linear scan is O(logN) since
	// predictor size is logN
	typename predictor_buffer_t::iterator x_it = predictor_.begin();
	while (x_it != predictor_.end() && (!eq(x_it->first,x))) ++x_it;

	if (x_it != predictor_.end()) {
	  // x is already a predictor -- update count
	  if (x_it->second >= segment_capacity_) {
	    typename predictor_buffer_t::iterator x_next = x_it;
	    ++x_next;
	    if (x_next!= predictor_.end()) {
	      // We are not the last element -- decrease that one
	      if (predictor_.back().second > 0) {
		--predictor_.back().second;
	      } 
	    }
	  } else {
	    // Increase this count
	    ++x_it->second;
	  }
	  // Move forward
	  if (x_it != predictor_.begin()) {
	    typename predictor_buffer_t::iterator x_prev = x_it;
	    --x_prev;
	    std::swap(*x_it,*x_prev);
	  }
	  // Pop back if needed
	  if (predictor_.back().second == 0) {
	    predictor_.pop_back();
	  }
	} else {
	  // x is not in element list -- insert as most recently updated
	  if (predictor_.full()) {
	    // No space -- decrease last
	    if (predictor_.back().second > 1) {
	      --predictor_.back().second;
	    } else {
	      predictor_.pop_back();
	      predictor_.push_front(std::make_pair(x,size_type(1)));
	    }
	  } else {
	    predictor_.push_front(std::make_pair(x,size_type(1)));
	  }
	}
      }
    }

    void apma_predictor_update_on_erase(size_type erase_index) {
      if (erase_index<capacity_ && exists_[erase_index]) {
	
	const value_type& x = storage_[erase_index];

	// Look for x in predictor buffer -- Linear scan is O(logN) since
	// predictor size is logN
	typename predictor_buffer_t::iterator x_it = predictor_.begin();
	while (x_it != predictor_.end() && (!eq(x_it->first,x))) ++x_it;
	
	if (x_it != predictor_.end()) {
	  // x is already a predictor -- Move the element at end and just remove the element...
	  // Should we instead push here a near one?
	  typename predictor_buffer_t::iterator x_next = x_it;
	  ++x_next;
	  while (x_next != predictor_.end()) {
	    std::swap(*x_it,*x_next);
	    ++x_it; ++x_next;
	  }
	  predictor_.pop_back();
	}
      }
    }
        
  protected:   

    /**
     * Return index in array of existing element not less than x
     * this means that insert position is at result
     */
    size_type apma_lower_bound(const value_type& x) const {
      assert(capacity_ == storage_.size());

      const size_type M = capacity_;

      if (element_count_ == 0) {
	return M;
      } else {
	size_type l= first_index_;
	size_type u= last_index_;
	assert(exists_[l]);
	assert(exists_[u]);
	if (lt(storage_[u],x)) {
	  // x is just after u
	  return M; // u+1; 
	} else if (!lt(storage_[l], x)) {
	  return l;
	} else {
	  while (l+8<u) {
	    assert(exists_[l]);
	    assert(exists_[u]);
	    
	    size_type m = (l+u)>>1;
	    while (!exists_[m]) ++m;
	    if (m==u) {
	      m = (l+u)>>1;
	      while (!exists_[m]) --m;
	    }
	    assert(exists_[m]);
	    if (lt(storage_[m],x)) {
	      // m<x -- search interval m+1..u
	      l=m+1;  while (!exists_[l]) ++l;
	    } else {
	      // m>=x -- search interval l..m
	      u=m;
	    }
	  }
	  while (l<u && lt(storage_[l], x)) {
	    ++l; while (!exists_[l]) ++l;
	  }
	  assert(exists_[l]);
	  assert(!lt(storage_[l],x));

	  return l;
	}
      }
    }
    
    /**
     * Return index in array of existing element greater than x
     * this means that insert position is at result
     */
    size_type apma_upper_bound(const value_type& x) const {
      assert(capacity_ == storage_.size());

      const size_type M = capacity_;

      if (element_count_ == 0) {
	return M;
      } else {
	size_type l= first_index_;
	size_type u= last_index_;
	assert(exists_[l]);
	assert(exists_[u]);

	if (!lt(x,storage_[u])) {
	  // out of bounds
	  return u+1;
	} else if (gt(x,storage_[l])) {
	  return l;
	} else {
	  while (l+8<u) {
	    assert(exists_[l]);
	    
	    size_type m = (l+u)>>1;
	    while (!exists_[m]) ++m;
	    if (m==u) {
	      m = (l+u)>>1;
	      while (!exists_[m]) --m;
	    }
	    assert(exists_[m]);
	    if (!lt(x,storage_[m])) {
	      // m<x -- search interval m+1..u
	      l=m+1;  while (!exists_[l]) ++l;
	    } else {
	      // m>=x -- search interval l..m
	      u=m;
	    }
	  }
	  while (l<u && !gt(storage_[l], x)) {
	    ++l; while (!exists_[l]) ++l;
	  }
	  assert(exists_[l]);
	  assert(gt(x,storage_[l]));

	  return l;
	}
      }
    }

    void apma_smallest_window_within_bounds_before_insert_in(size_type* window_start,
							     size_type* window_size,
							     size_type  insert_index) const {
      assert(element_count_ > 0);
      assert(capacity_ == storage_.size());

      const size_type M = capacity_;

      if (insert_index == M) {
	// We are inserting after last element -- move insert position
	// back since we are going to rebalance the array
	--insert_index;
      }
            
      // Known window
      size_type rstart  = insert_index;
      size_type rend    = insert_index;
      size_type rcount  = (exists_[insert_index] ? 1 : 0);

      bool is_balanced  = false;
      size_type level   = 0; // Next level to check -- we start at leaf

      while ((!is_balanced) && level <= (implicit_tree_root_height_)) {
	assert(level <= implicit_tree_root_height_);
	
	// Get the boundaries of the next window
	size_type rsz = apma_window_capacity(level);
	assert(rsz<=M);
	
	size_type rleft  = rstart-(rstart%rsz);
	size_type rright = rleft+rsz-1;
	
	if (level == implicit_tree_root_height_) {
	  assert(rsz == capacity_);
	  assert(rleft == 0);
	  assert(rright+1 == capacity_);
	  // Use top level stats
	  rcount = element_count_;
	} else {
	  // Update count adding new unvisited windows
	  for (size_type i = rleft; i < rstart; ++i) {
	    assert(i<M);
	    if (exists_[i]) ++rcount;
	  }
	  for (size_type i = rend + 1; i <= rright; ++i) {
	    assert(i<M);
	    if (exists_[i]) ++rcount;
	  }
	}
	rstart = rleft; rend = rright;
	
	is_balanced = !apma_is_window_too_full(level, rcount+1); 
	
	++level;
      }
      
      if (is_balanced) {
	*window_start = rstart;
	*window_size = rend+1-rstart;
      } else {
	// Needs a rebuild!
	*window_start = 0;
	*window_size = 0;
      }
      
      assert(*window_start+*window_size<=capacity_);
      assert((*window_size==0) || ((insert_index>= *window_start) && (insert_index< *window_start+*window_size)));
    }

    void apma_smallest_window_within_bounds_after_erase_in(size_type* window_start,
							   size_type* window_size,
							   size_type  erase_index_first,
							   size_type  erase_index_last) const {
      assert(capacity_ == storage_.size());

      const size_type M = capacity_; SL_USEVAR(M);

      // Known window
      size_type rstart  = erase_index_first;
      size_type rend    = erase_index_last;
      size_type rcount  = 0;

      bool is_balanced  = false;
      size_type level   = 0; // Next level to check -- we start at leaf
      // Grow until next level contains erased interval
      bool bracketed = false;
      while (!bracketed) {
	size_type rsz = apma_window_capacity(level);
	size_type rleft  = rstart-(rstart%rsz);
	size_type rright = rleft+rsz-1;
	bracketed = (rleft <= rstart && rend <= rright);
	++level;
      }
      assert(level>0);
      --level;
      assert(level<=implicit_tree_root_height_);

      while ((!is_balanced) && level <= (implicit_tree_root_height_)) {
	assert(level <= implicit_tree_root_height_);
	
	// Get the boundaries of the next window
	size_type rsz = apma_window_capacity(level);
	assert(rsz<=M);
	
	size_type rleft  = rstart-(rstart%rsz);
	size_type rright = rleft+rsz-1;
	
	if (level == implicit_tree_root_height_) {
	  assert(rsz == capacity_);
	  // Use top level stats
	  rcount = element_count_;
	} else {
	  // Update count adding new unvisited windows
	  for (size_type i = rleft; i < rstart; ++i) {
	    assert(i<M);
	    if (exists_[i]) ++rcount;
	  }
	  for (size_type i = rend + 1; i <= rright; ++i) {
	    assert(i<M);
	    if (exists_[i]) ++rcount;
	  }
	}
	rstart = rleft; rend = rright;
	
	is_balanced = !apma_is_window_too_empty(level, rcount); 
	
	++level;
      }
      
      if (is_balanced) {
	*window_start = rstart;
	*window_size = rend+1-rstart;
      } else {
	// Needs a rebuild!
	*window_start = 0;
	*window_size = 0;
      }

      assert(*window_start+*window_size<=capacity_);
      assert((*window_size==0) || ((erase_index_first>= *window_start) && (erase_index_last< *window_start+*window_size)));
    }

  protected: // Rebalance and reallocate
    
    /**
     * Redistribute elements in window evenly, inserting a
     * new element if insert_index != size_type(-1).
     */
    size_type apma_rebalance(size_type window_start,
			     size_type window_size,
			     size_type insert_index = size_type(-1),
			     const value_type* x_insert = 0) {
      size_type N_prime = element_count_;
      if (insert_index!=size_type(-1)) ++N_prime;

      size_type packed_count;
      size_type packed_insert_index;
      apma_pack_right(&packed_count, &packed_insert_index,
		      window_start, window_size,
		      insert_index, x_insert);
      size_type packed_start = window_start+window_size-packed_count;

      assert(element_count_ == N_prime);
      assert((window_size!=capacity_) || (packed_count == element_count_));
      assert(insert_index != size_type(-1) || (packed_insert_index==size_type(-1)));
      assert(insert_index == size_type(-1) || (packed_insert_index<capacity_));
      assert(insert_index == size_type(-1) || (exists_[packed_insert_index]));
      assert(insert_index == size_type(-1) || (eq(*x_insert, storage_[packed_insert_index])));

      size_type moved_insert_index;
#if 0
      apma_even_redistribute_packed(&moved_insert_index,
				    window_start, window_size,
				    packed_start, packed_count,
				    packed_insert_index);
#else
      apma_uneven_redistribute_packed(&moved_insert_index,
				      window_start, window_size,
				      packed_start, packed_count,
				      packed_insert_index);
#endif
      assert(insert_index != size_type(-1) || (moved_insert_index==size_type(-1)));
      assert(insert_index == size_type(-1) || (moved_insert_index<capacity_));
      assert(insert_index == size_type(-1) || (exists_[moved_insert_index]));
      assert(insert_index == size_type(-1) || (eq(*x_insert, storage_[moved_insert_index])));
      
      return moved_insert_index; // insert position
    }
    
    /**
     * Reallocate array to at least new capacity and redistribute elements evenly,
     * inserting a new element if insert_index != size_type(-1).
     */
    size_type apma_reallocate(size_type new_capacity,
			      size_type insert_index = size_type(-1),
			      const value_type* x_insert = 0) {
      assert((insert_index == size_type(-1)) || (x_insert != 0));
      assert(capacity_ == storage_.size());
       
      apma_init_tree_shape(new_capacity);
      assert(capacity_>=new_capacity);
      
      size_type N_prime = element_count_;
      if (insert_index!=size_type(-1)) ++N_prime;
      assert(N_prime<=capacity_);

      storage_.resize(capacity_);
      exists_.resize(capacity_, false);
      predictor_.set_capacity(segment_capacity_); // FIXME

      return apma_rebalance(0, capacity_, insert_index, x_insert);
    }
    
  protected:

    size_type apma_insert(size_type insert_index, const value_type& x) {
      assert(element_count_>0);
      
      const size_type M = capacity_;
      
      size_type result = M;

      if (insert_index == M) {
	insert_index = last_index_+1;
	assert(insert_index==M || !exists_[insert_index]);
      }

      apma_predictor_update_on_insert(insert_index);

      if ((insert_index<M) && !exists_[insert_index]) {
	// Slot is empty -- just store
	storage_[insert_index] = x;
	exists_[insert_index] = true;
	if (insert_index<first_index_) first_index_ = insert_index;
	if (insert_index>last_index_) last_index_ = insert_index;
	result = insert_index;
	++element_count_;
      } else if ((insert_index>0) && !exists_[insert_index-1]) { 
	// Empty slot just before insertion position -- just store
	size_type slot_index = insert_index-1;
	storage_[slot_index] = x;
	exists_[slot_index] = true;
	if (slot_index<first_index_) first_index_ = slot_index;
	if (slot_index>last_index_) last_index_ = slot_index;
	result = slot_index;
	++element_count_;
      } else {
	const size_type MAX_OFFSET = 8;

	// Look for empty slot in small area around desired position
	size_type offset_left = 1;
	while (insert_index>offset_left && exists_[insert_index-offset_left] && offset_left<MAX_OFFSET) {
	  ++offset_left;
	}
	bool is_valid_offset_left = insert_index>=offset_left && !exists_[insert_index-offset_left];
	
	size_type offset_right = 1;
	while (insert_index+offset_right<M && exists_[insert_index+offset_right] && offset_right<MAX_OFFSET) {
	  ++offset_right;
	}
	bool is_valid_offset_right = insert_index+offset_right<M && !exists_[insert_index+offset_right];

	bool is_best_left  = is_valid_offset_left  && (!is_valid_offset_right || (offset_left<offset_right));
	bool is_best_right = (!is_best_left) && is_valid_offset_right;

	if (is_best_left) {
	  // LOW COST insert to the left
	  assert(!exists_[insert_index-offset_left]);
	  for (size_type i = offset_left; i>1; --i) {
	    storage_[insert_index-i] = storage_[insert_index-(i-1)];
	  }
	  exists_[insert_index-offset_left] = true;
	  size_type slot_index = insert_index-1;
	  storage_[slot_index] = x;
	  exists_[slot_index] = true;
	  if (insert_index-offset_left<first_index_) first_index_ = insert_index-offset_left;
	  result = slot_index;
	  ++element_count_;
 	} else if (is_best_right) {
	  // LOW COST insert to the right
	  assert(!exists_[insert_index+offset_right]);
	  for (size_type i = offset_right; i>0; --i) {
	    storage_[insert_index+i] = storage_[insert_index+(i-1)];
	  }
	  exists_[insert_index+offset_right] = true;
	  size_type slot_index = insert_index;
	  storage_[slot_index] = x;
	  exists_[slot_index] = true;
	  if (insert_index+offset_right>last_index_) last_index_ = insert_index+offset_right;
	  result = slot_index;
	  ++element_count_;
	} else {
	  // HEAVY insert
	  // Find area with low enough density and insert into it
	  size_type window_start;
	  size_type window_size;
	  apma_smallest_window_within_bounds_before_insert_in(&window_start,
							      &window_size,
							      insert_index);
	  if (window_size== 0) {
	    // Whole array is too dense -- resize it!
	    size_type new_capacity = APMA_SCALE_FACTOR()*capacity_;
	    result = apma_reallocate(new_capacity, insert_index, &x);
	  } else {
	    result = apma_rebalance(window_start, window_size, insert_index, &x);
	  }
	}	  
      } 

      assert(result<capacity_);
      assert(exists_[result]);
      assert(eq(storage_[result],x));
      
      return result;
    }
			 
    size_type apma_insert(const value_type& x) {
      const size_type N = element_count_;

      size_type result;
      if (N == 0) {
	// Insert into empty array
	result = apma_init(x);
      } else {
	// Find insertion position into current array and insert
	result = apma_insert(apma_lower_bound(x), x);
      }
      return result;
    }

    size_type apma_erase(size_type erase_index_first, size_type erase_index_last) {
      const size_type M = capacity_;

      size_type result = 0; 
      // Remove all elements in range
      for (size_type idx=erase_index_first; idx<= erase_index_last; ++idx) {
	bool b = this->exists_[idx];
	if (b) {
	  apma_predictor_update_on_erase(idx);
	  this->exists_[idx] = false;
	  --element_count_;
	  ++result;
	}
      }
      
      const size_type N = element_count_;
      if (N == 0) {
	// Wipe out array
	apma_clear();
      } else {
	// Update first/last bounds
	if (first_index_ >= erase_index_first) {
	  first_index_ = erase_index_last+1;
	  while (!exists_[first_index_]) ++first_index_;
	}
	if (last_index_ <= erase_index_last) {
	  last_index_ = erase_index_first-1;
	  while (!exists_[last_index_]) --last_index_;
	}
	// Find area with high enough density and reblance it
	size_type window_start;
	size_type window_size;
	apma_smallest_window_within_bounds_after_erase_in(&window_start,
							  &window_size,
							  erase_index_first,
							  erase_index_last);
	if (window_size == 0) {
	  // Whole array is too empty -- resize it!
	  size_type new_capacity = M;
	  while (N <= APMA_ROOT_UPPER_DENSITY() * (new_capacity/APMA_SCALE_FACTOR())) new_capacity/= APMA_SCALE_FACTOR();
	  apma_reallocate(new_capacity);
	} else {
	  apma_rebalance(window_start, window_size);
	}
      }

      return result; // Number of erased items
    }

    size_type apma_erase(size_type erase_index_first) {
      return apma_erase(erase_index_first, erase_index_first);
    }
    
  public:

    adaptive_packed_memory_array(const Compare& pred = Compare(),const allocator_type& al = allocator_type()):
      comp_(pred), storage_(al)  {
      apma_clear();
    }

    template <class It>
    adaptive_packed_memory_array(It first, It beyond,
				 const Compare& pred = Compare(),const allocator_type& al = allocator_type()) :
      comp_(pred), storage_(al)  {
      this->assign(first, beyond);
    }

    adaptive_packed_memory_array(const this_t& other)
    : comp_(other.comp_),
      implicit_tree_root_height_(other.implicit_tree_root_height_),
      segment_capacity_(other.segment_capacity_),
      capacity_(other.capacity_),
      storage_(other.storage_),
      exists_(other.exists_),
      predictor_(other.predictor_),
      element_count_(other.element_count_),
      first_index_(other.first_index_),
      last_index_(other.last_index_)
    {
    }

    ~adaptive_packed_memory_array() {
      apma_clear();
    }

    this_t& operator=(const this_t& other) {
      (this->comp_) = (other.comp_);
      (this->implicit_tree_root_height_) = (other.implicit_tree_root_height_);
      (this->segment_capacity_) = (other.segment_capacity_);
      (this->capacity_) = (other.capacity_);
      (this->storage_) =(other.storage_);
      (this->exists_) = (other.exists_);
      (this->predictor_) = (other.predictor_);
      (this->element_count_) = (other.element_count_);
      (this->first_index_) = (other->first_index_);
      (this->last_index_) = (other->last_index_);
      return *this;
    }
    
    void reserve(size_type n) {
      // FIXME --- possibly reallocate storage instead
      (this->storage_).reserve(n);
      (this->exists_).reserve(n);
      (this->predictor_).reserve(log2(n)); // FIXME??
    }
    
    allocator_type get_allocator() const {
      return (this->storage_).get_allocator();
    }
    

  public: // swapping
    
    void swap(this_t& other) {
      std::swap((this->comp_),other.comp_);
      std::swap(this->implicit_tree_root_height_, other.implicit_tree_root_height_);
      std::swap(this->segment_capacity_, other.segment_capacity_);
      std::swap(this->capacity_, other.capacity_);
      (this->storage_).swap(other.storage_);
      (this->exists_).swap(other.exists_);
      std::swap(this->element_count_, other.element_count_);
      std::swap(this->first_index_, other.first_index_);
      std::swap(this->last_index_, other.last_index_);
    }
    
    friend void swap(this_t& x, this_t& y) {
      x.swap(y);
    }

  public: // Comparison
    
    key_compare key_comp() const {
      return (this->comp_);
    }
    
    value_compare value_comp() const {
      return (this->key_comp());
    }

  public: // Iterators
    
    inline iterator begin() {
      return this->capacity_ ? iterator(this, first_index_) : end();
    }
    
    inline const_iterator begin() const {
      return this->capacity_ ? const_iterator(this, first_index_) : end();
    }
    
    inline iterator end() {
      return iterator(this, this->capacity_);
    }
    
    inline const_iterator end() const {
      return const_iterator(this, this->capacity_);
    }
    
    inline reverse_iterator rbegin() {
      return this->capacity_ ? reverse_iterator(this, capacity_-last_index_-1) : rend();
    }
    
    inline const_reverse_iterator rbegin() const {
      return this->capacity_ ? const_reverse_iterator(this, capacity_-last_index_-1) : rend();
    }
    
    inline reverse_iterator rend() {
      return reverse_iterator(this, this->capacity_);
    }
    
    inline const_reverse_iterator rend() const {
      return const_reverse_iterator(this, this->capacity_);
    }

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      store_container(s, *this);
    }
    
    void retrieve_from(input_serializer& s) {
      apma_clear();
      retrieve_container(s, *this);
    }
    
  public: // Size
    
    inline size_type size() const {
      return (this->element_count_);
    }
    
    inline size_type max_size() const {
      return (this->storage_).max_size();
    }
    
    inline bool empty() const {
      return (this->element_count_) == 0;
    }

  public: // Low-level element access

    inline bool storage_is_valid_index(size_type idx) const {
      return idx<(this->capacity_);
    }
    
    inline size_type storage_size() const {
      return (this->capacity_);
    }

    inline size_type storage_first_index() const {
      return this->first_index_;
    }

    inline size_type storage_last_index() const {
      return this->last_index_;
    }

    inline size_type storage_mid_index(size_type l,
				       size_type u) const {
      assert(this->is_valid_index(l));
      assert(this->is_valid_index(u));
      assert(this->storage_exists(l));
      assert(this->storage_exists(u));

      size_type m = (l+u)>>1;
      while (!exists_[m]) ++m;
      if (m==u) {
	m = (l+u)>>1;
	while (!exists_[m]) --m;
      }
      return m;
    }
    
    inline bool storage_is_present(size_type idx) const {
      assert(storage_is_valid_index(idx));
      return (this->exists_)[idx];
    }
        
    inline const_reference storage_at(size_type idx) const {
      assert(storage_is_valid_index(idx));
      return (this->storage_)[idx];
    }
            
    inline reference storage_at(size_type idx) {
      assert(storage_is_valid_index(idx));
      return (this->storage_)[idx];
    }
            
    inline const_reference front() const {
      assert(!empty());
      return this->storage_[this->storage_first_index()];
    }

    inline const_reference back() const {
      assert(!empty());
      return this->storage_[this->storage_last_index()];
    }

  public: // Front/back removal

    void pop_front() {
      assert(!this->empty());
      this->apma_erase(this->storage_first_index());
    }
    
    void pop_back() {
      assert(!this->empty());
      this->apma_erase(this->storage_last_index());
    }
    
  public: // Search
        
    inline iterator lower_bound(const value_type& x) {
      return iterator(this, this->apma_lower_bound(x));
    }
    
    inline const_iterator lower_bound(const value_type& x) const {
      return const_iterator(this, this->apma_lower_bound(x));
    }

    inline iterator find(const value_type& x) {
      iterator it = this->lower_bound(x);
      return
	((it==this->end())||(this->lt(x, *it))) ? this->end() : it;
    }
    
    inline const_iterator find(const value_type& x) const {
      const_iterator it = this->lower_bound(x);
      return
	((it==this->end())||(this->lt(x, *it))) ? this->end() : it;
    }
    
    inline iterator upper_bound(const value_type& x) {
      return iterator(this, apma_upper_bound(x));
    }
    
    inline const_iterator upper_bound(const value_type& x) const {
      return const_iterator(this, apma_upper_bound(x));
    }
        
    inline std::pair<iterator, iterator> equal_range(const value_type& x) {
      iterator lb = this->lower_bound(x);
      iterator ub = lb; if (lb!=this->end()) { while (ub != this->end() && !this->gt(*ub,x)) ++ub; }
      return std::make_pair(lb,ub);
    }
    
    inline std::pair<const_iterator, const_iterator> equal_range(const value_type& x) const {
      const_iterator lb = this->lower_bound(x);
      const_iterator ub = lb; if (lb!=this->end()) { while (ub != this->end() && !this->gt(*ub,x)) ++ub; }
      return std::make_pair(lb,ub);
    }
    
    inline size_type count(const value_type& x) const {
      size_type result = 0;
      const_iterator it = find(x);
      while (it != end() && !this->gt(*it,x)) {
	++result; ++it;
      }
      return result;
    }

  public: // Assignment
    
    template <class It>
    void assign(It first, It beyond) {
      this->apma_clear();
      this->insert(first, beyond);
    }
    
    void assign(size_type n, const value_type& x = value_type()) {
      this->apma_clear();
      if (is_not_allowing_duplicates) {
	if (n) this->insert(x);
      } else {
	for (size_type i=0; i<n; ++i) {
	  this->insert(x);
	}
      } 
    }
    
  public: // Element insertion
    
    std::pair<iterator,bool> insert(const value_type& x) {
      if (element_count_ == 0) {
	// Insert into empty array
	apma_init(x);
	return std::make_pair(begin(), true);
      } else {
	iterator it = lower_bound(x);
	bool do_insert =
	  (!is_not_allowing_duplicates) ||
	  (it == end()) ||
	  (this->lt(x,*it));
	if (do_insert) {
	  return std::make_pair(iterator(this, apma_insert(it.index(), x)), true);
	} else {
	  return std::make_pair(it, false);
	}
      }
    }

    iterator insert(iterator hint_it, const value_type& x) {
      if (element_count_ == 0) {
	// Insert into empty array
	apma_init(x);
	return std::make_pair(begin(), true);
      } else {
	if (hint_it==end()) {
	  hint_it = iterator(this, storage_last_index());
	}
	assert(hint_it!=end());
	
	iterator lb_it = end();
	if (!lt(*hint_it,x)) {
	  // Might be the lower bound -- check
	  if (hint_it==begin()) {
	    // The first one
	    lb_it = hint_it;
	  } else {
	    iterator prev_hint_it = hint_it; --prev_hint_it;
	    if (lt(*prev_hint_it,x)) {
	      lb_it = hint_it;
	    } else {
	      lb_it = lower_bound(x);
	    }
	  }
	} else {
	  // The next one might be the lower bound -- check
	  iterator next_hint_it = hint_it; ++next_hint_it;
	  if (next_hint_it == end() || lt(*next_hint_it,x)) {
	    lb_it = next_hint_it;
	  } else {
	    lb_it = lower_bound(x);
	  }
	}
	bool do_insert =
	  (!is_not_allowing_duplicates) ||
	  (lb_it == end()) ||
	  (this->lt(x,*lb_it)); // == !equal
	if (do_insert) {
	  return iterator(this, apma_insert(lb_it.index(), x));
	} else {
	  return lb_it;
	}
      }
    }
	    
    template<class It>
    void insert(It first, It beyond) {
      // reserve(size() + std::distance(first, beyond));
      for (It it=first; it!= beyond; ++it) {
	insert(*it);
      }
    }

  public: // Element removal
    
    void erase(iterator it) {
      if (it != this->end()) {
	this->apma_erase(it.index());
      }
    }

    void erase(reverse_iterator it) {
      if (it != this->rend()) {
	this->apma_erase(it.index());
      }
    }
    
    size_type erase(iterator first, iterator beyond) {
      size_type result = 0;      
      if (first != beyond) {
	size_type first_index = first.index();
	size_type last_index = std::max(beyond.index(),size_type(1))-size_type(1);
	result = this->apma_erase(std::max(first_index,first_index_), std::min(last_index,last_index_));
      }
      return result;
    }

    size_type erase(reverse_iterator first, reverse_iterator beyond) {
      size_type result = 0;      
      if (first != beyond) {
	size_type last_index = first.index();
	size_type first_index = (beyond == this->rend()) ? size_type(0) : std::max(beyond.index(),size_type(1))-size_type(1);
	result = this->apma_erase(std::max(first_index,first_index_), std::min(last_index,last_index_));
      }
      return result;
    }

    size_type erase(const value_type& x) {
      if (is_not_allowing_duplicates) {
	iterator it = this->find(x);
	if (it == this->end()) {
	  return 0;
	} else {
	  return apma_erase(it.index(),it.index());
	}
      } else {
	std::pair<iterator,iterator> beg_end = this->equal_range(x);
	return this->erase(beg_end.first,
			   beg_end.second);
      }
    }
    
    void clear() {
      apma_clear();
    }

  public: // comparison
    
    bool is_equal(const this_t& other) const {
      return (size() == other.size()
              && std::equal(begin(), end(), other.begin()));
    }
    
    bool is_less_than(const this_t& other) const {
      return (std::lexicographical_compare(begin(), end(),
                                           other.begin(), other.end()));
    }
    
  };

}

template<class K,bool is_not_allowing_duplicates,class Compare, class A> 
inline bool operator==(const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& x,
		       const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& y) {
  return x.is_equal(y);
}

template<class K,bool is_not_allowing_duplicates,class Compare, class A> 
inline bool operator!=(const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& x,
		       const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& y) {
  return !(x == y);
}

template<class K,bool is_not_allowing_duplicates,class Compare, class A> 
inline bool operator<(const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& x,
		      const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& y) {
  return x.is_less_than(y);
}

template<class K,bool is_not_allowing_duplicates,class Compare,class A> 
inline bool operator>(const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& x,
		      const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& y) {
  return y < x;
}

template<class K,bool is_not_allowing_duplicates,class Compare, class A>
inline bool operator<=(const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& x,
		       const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& y) {
  return !(y < x);
}

template<class K, bool is_not_allowing_duplicates,class Compare,class A>
bool operator>=(const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& x,
		const sl::adaptive_packed_memory_array<K, is_not_allowing_duplicates,Compare,A>& y) {
  return (!(x < y));
}


// =============================================================================================
// STL Compliant implementations -- set 
// =============================================================================================

namespace sl {
  
  /**
   *  A std::set implemented as an adaptive packed-memory array
   */
  template <class T, class Compare = std::less<T>, class Allocator = std::allocator<T> >
  class adaptive_packed_memory_array_set: public adaptive_packed_memory_array<T,true,Compare,Allocator> {
  public:
    typedef adaptive_packed_memory_array_set<T,Compare,Allocator>  this_t;
    typedef adaptive_packed_memory_array<T,true,Compare,Allocator> super_t;
    
    typedef typename super_t::value_vector_t value_vector_t;
    typedef typename super_t::bool_vector_t  bool_vector_t;

    typedef typename super_t::allocator_type  allocator_type;
    typedef typename super_t::size_type       size_type;
    typedef typename super_t::difference_type difference_type;
    typedef typename super_t::reference	      reference;
    typedef typename super_t::const_reference const_reference;
    typedef typename super_t::value_type      value_type;
    typedef typename super_t::value_type      key_type;
    typedef typename super_t::key_compare     key_compare;
    typedef typename super_t::value_compare   value_compare;

    typedef typename super_t::iterator                iterator;
    typedef typename super_t::const_iterator          const_iterator;
    typedef typename super_t::reverse_iterator        reverse_iterator;
    typedef typename super_t::const_reverse_iterator  const_reverse_iterator;

  public:
    
    adaptive_packed_memory_array_set(const Compare& pred = Compare(),const allocator_type& al = allocator_type()):
      super_t(pred, al)  {
    }

    template <class It>
    adaptive_packed_memory_array_set(It first, It beyond,
				     const Compare& pred = Compare(), const allocator_type& al = allocator_type()) :
      super_t(first, beyond, pred, al)  {
    }

    adaptive_packed_memory_array_set(const this_t& other)
    : super_t(other) {
    }

    ~adaptive_packed_memory_array_set() {
      this->apma_clear();
    }

    this_t& operator=(const this_t& other) {
      super_t::operator=(other);
      return *this;
    }
  };
} // namespace sl

// =============================================================================================
// STL Compliant implementations -- multiset
// =============================================================================================

namespace sl {
  /**
   *  A std::multiset implemented as an adaptive packed-memory array
   */
  template <class T, class Compare = std::less<T>, class Allocator = std::allocator<T> >
  class adaptive_packed_memory_array_multiset: public adaptive_packed_memory_array<T,false,Compare,Allocator> {
  public:
    typedef adaptive_packed_memory_array_multiset<T,Compare,Allocator>  this_t;
    typedef adaptive_packed_memory_array<T,false,Compare,Allocator> super_t;
    
    typedef typename super_t::value_vector_t value_vector_t;
    typedef typename super_t::bool_vector_t  bool_vector_t;

    typedef typename super_t::allocator_type  allocator_type;
    typedef typename super_t::size_type       size_type;
    typedef typename super_t::difference_type difference_type;
    typedef typename super_t::reference	      reference;
    typedef typename super_t::const_reference const_reference;
    typedef typename super_t::value_type      value_type;
    typedef typename super_t::value_type      key_type;
    typedef typename super_t::key_compare     key_compare;
    typedef typename super_t::value_compare   value_compare;

    typedef typename super_t::iterator                iterator;
    typedef typename super_t::const_iterator          const_iterator;
    typedef typename super_t::reverse_iterator        reverse_iterator;
    typedef typename super_t::const_reverse_iterator  const_reverse_iterator;

  public:
    
    adaptive_packed_memory_array_multiset(const Compare& pred = Compare(),const allocator_type& al = allocator_type()):
      super_t(pred, al)  {
    }

    template <class It>
    adaptive_packed_memory_array_multiset(It first, It beyond,
					  const Compare& pred = Compare(),const allocator_type& al = allocator_type()) :
      super_t(first, beyond, pred, al)  {
    }

    adaptive_packed_memory_array_multiset(const this_t& other)
    : super_t(other) {
    }

    ~adaptive_packed_memory_array_multiset() {
      this->apma_clear();
    }

    this_t& operator=(const this_t& other) {
      super_t::operator=(other);
      return *this;
    }
  };  
} // namespace sl

// =============================================================================================
// STL Compliant implementations -- map 
// =============================================================================================

namespace sl {

  template <class Key, class T, class Compare = std::less<Key> >
  class map_value_type_compare {
  public:
    typedef Key                     key_type;
    typedef T                       mapped_type;
    typedef std::pair</*const*/Key, T> value_type;
  public:
    Compare comp_;
  public:
    inline map_value_type_compare(Compare c) : comp_(c) {}
    inline bool operator() (const value_type& x, const value_type& y) const {
      return comp_(x.first, y.first);
    }
  };
  
  template < class Key, class T,
	     class Compare = std::less<Key>,
	     class Allocator = std::allocator<std::pair</*const*/ Key,T> > >
  class adaptive_packed_memory_array_map: public adaptive_packed_memory_array<std::pair</*const*/ Key,T>,
									      true,
									      map_value_type_compare<Key,T,Compare>,
									      Allocator> {
  public:
    typedef Key                     key_type;
    typedef T                       mapped_type;
    typedef std::pair</*const*/Key, T> value_type;
    typedef map_value_type_compare<Key,T,Compare> value_compare;
  public:
    typedef adaptive_packed_memory_array_map<Key,T,Compare,Allocator>               this_t;
    typedef adaptive_packed_memory_array<value_type,true,value_compare,Allocator> super_t;
    
    typedef typename super_t::value_vector_t value_vector_t;
    typedef typename super_t::bool_vector_t  bool_vector_t;

    typedef typename super_t::allocator_type  allocator_type;
    typedef typename super_t::size_type       size_type;
    typedef typename super_t::difference_type difference_type;
    typedef typename super_t::reference	      reference;
    typedef typename super_t::const_reference const_reference;
    typedef Compare                           key_compare;

    typedef typename super_t::iterator                iterator;
    typedef typename super_t::const_iterator          const_iterator;
    typedef typename super_t::reverse_iterator        reverse_iterator;
    typedef typename super_t::const_reverse_iterator  const_reverse_iterator;

  public:
    
    adaptive_packed_memory_array_map(const Compare& pred = Compare(), const allocator_type& al = allocator_type()):
      super_t(value_compare(pred), al)  {
    }

    template <class It>
    adaptive_packed_memory_array_map(It first, It beyond,
				     const Compare& pred = Compare(),const allocator_type& al = allocator_type()) :
      super_t(first, beyond, value_compare(pred), al)  {
    }

    adaptive_packed_memory_array_map(const this_t& other)
    : super_t(other) {
    }

    ~adaptive_packed_memory_array_map() {
      this->apma_clear();
    }

    this_t& operator=(const this_t& other) {
      super_t::operator=(other);
      return *this;
    }

    key_compare key_comp() const {
      return this->comp_.comp_;
    }

    const value_compare& value_comp() const {
      return this->comp_;
    }
 
    size_type erase(const key_type& x) {
      return super_t::erase(value_type(x,mapped_type()));
    }

    std::pair<iterator,iterator> equal_range(const key_type& x) {
      return super_t::equal_range(value_type(x,mapped_type()));
    }
    
    std::pair<const_iterator,const_iterator> equal_range(const key_type& x) const {
      return super_t::equal_range(value_type(x,mapped_type()));
    }
       
    iterator find(const key_type& x) {
      return super_t::find(value_type(x,mapped_type()));
    }
    
    const_iterator find(const key_type& x) const {
      return super_t::find(value_type(x,mapped_type()));
    }

    iterator lower_bound(const key_type& x) {
      return super_t::lower_bound(value_type(x,mapped_type()));
    }
    
    const_iterator lower_bound(const value_type& x) const {
      return super_t::lower_bound(value_type(x,mapped_type()));
    }
    
    iterator upper_bound(const key_type& x) {
      return super_t::upper_bound(value_type(x,mapped_type()));
    }
    
    const_iterator upper_bound(const value_type& x) const {
      return super_t::upper_bound(value_type(x,mapped_type()));
    }

    size_type count(const key_type& x) const {
      return super_t::count(value_type(x,mapped_type()));
    }

    mapped_type& operator[](const key_type& x) {
      return (*((this->insert(std::make_pair(x,mapped_type()))).first)).second;
    }
  };

}
// =============================================================================================
// STL Compliant implementations -- multimap 
// =============================================================================================

namespace sl {

  template < class Key, class T,
	     class Compare = std::less<Key>,
	     class Allocator = std::allocator<std::pair</*const*/ Key,T> > >
  class adaptive_packed_memory_array_multimap: public adaptive_packed_memory_array<std::pair</*const*/ Key,T>,
										   false,
										   map_value_type_compare<Key,T,Compare>,
										   Allocator> {
  public:
    typedef Key                     key_type;
    typedef T                       mapped_type;
    typedef std::pair</*const*/Key, T> value_type;
    typedef map_value_type_compare<Key,T,Compare> value_compare;
  public:
    typedef adaptive_packed_memory_array_multimap<Key,T,Compare,Allocator>         this_t;
    typedef adaptive_packed_memory_array<value_type,false,value_compare,Allocator> super_t;
    
    typedef typename super_t::value_vector_t value_vector_t;
    typedef typename super_t::bool_vector_t  bool_vector_t;

    typedef typename super_t::allocator_type  allocator_type;
    typedef typename super_t::size_type       size_type;
    typedef typename super_t::difference_type difference_type;
    typedef typename super_t::reference	      reference;
    typedef typename super_t::const_reference const_reference;
    typedef Compare                           key_compare;

    typedef typename super_t::iterator                iterator;
    typedef typename super_t::const_iterator          const_iterator;
    typedef typename super_t::reverse_iterator        reverse_iterator;
    typedef typename super_t::const_reverse_iterator  const_reverse_iterator;

  public:
    
    adaptive_packed_memory_array_multimap(const Compare& pred = Compare(), const allocator_type& al = allocator_type()):
      super_t(value_compare(pred), al)  {
    }

    template <class It>
    adaptive_packed_memory_array_multimap(It first, It beyond,
				     const Compare& pred = Compare(),const allocator_type& al = allocator_type()) :
      super_t(first, beyond, value_compare(pred), al)  {
    }

    adaptive_packed_memory_array_multimap(const this_t& other)
    : super_t(other) {
    }

    ~adaptive_packed_memory_array_multimap() {
      this->apma_clear();
    }

    this_t& operator=(const this_t& other) {
      super_t::operator=(other);
      return *this;
    }

    key_compare key_comp() const {
      return this->comp_.comp_;
    }

    const value_compare& value_comp() const {
      return this->comp_;
    }
 
    size_type erase(const key_type& x) {
      return super_t::erase(value_type(x,mapped_type()));
    }

    std::pair<iterator,iterator> equal_range(const key_type& x) {
      return super_t::equal_range(value_type(x,mapped_type()));
    }
    
    std::pair<const_iterator,const_iterator> equal_range(const key_type& x) const {
      return super_t::equal_range(value_type(x,mapped_type()));
    }
       
    iterator find(const key_type& x) {
      return super_t::find(value_type(x,mapped_type()));
    }
    
    const_iterator find(const key_type& x) const {
      return super_t::find(value_type(x,mapped_type()));
    }

    iterator lower_bound(const key_type& x) {
      return super_t::lower_bound(value_type(x,mapped_type()));
    }
    
    const_iterator lower_bound(const value_type& x) const {
      return super_t::lower_bound(value_type(x,mapped_type()));
    }
    
    iterator upper_bound(const key_type& x) {
      return super_t::upper_bound(value_type(x,mapped_type()));
    }
    
    const_iterator upper_bound(const value_type& x) const {
      return super_t::upper_bound(value_type(x,mapped_type()));
    }

    size_type count(const key_type& x) const {
      return super_t::count(value_type(x,mapped_type()));
    }

  };

}

#endif
