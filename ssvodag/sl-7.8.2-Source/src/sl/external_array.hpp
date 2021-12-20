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
#ifndef SL_EXTERNAL_ARRAY_HPP
#define SL_EXTERNAL_ARRAY_HPP

# ifdef _MSC_VER
#   pragma warning(disable:4522)
# endif //_MSC_VER

#include <sl/assert.hpp>
#include <sl/cstdint.hpp>
#include <sl/generative_types.hpp>
#include <sl/operators.hpp>
#include <sl/os_file.hpp>
#include <sl/sort.hpp>
#include <cassert>

/**
 * FIXME:
 *  - construct elements at resize, current version leaves
 *    them not initialized!
 *  - simplify creation with factory methods
 *  - parameterize mmapping
 */

namespace sl {
  
  template <class T>
  class external_array1;

  template <class T>
  class external_array1_item_ref;

  template <class T>
  class external_array1_item_const_ref;
  
  /// Const reference to external array items
  template <class T>
  class external_array1_item_const_ref {
  public:
    typedef external_array1_item_const_ref<T> self_t;
    typedef external_array1<T> array_t;
    typedef T                  value_type;
    typedef uint64_t           size_type;
  protected:
    const array_t *xarray_a_;
    size_type xarray_idx_;
    
  public:
    
    inline external_array1_item_const_ref(const array_t* a, size_type idx)
        : xarray_a_(a), xarray_idx_(idx) {
    }
 
    inline size_type index() const {
      return xarray_idx_;
    }

    inline value_type value() const {
      return xarray_a_->item(xarray_idx_);
    }

    inline operator value_type() const {
      return value();
    }
  };

  /// Reference to external array items
  template <class T>
  class external_array1_item_ref {
  public:
    typedef external_array1_item_ref<T> self_t;
    typedef external_array1<T> array_t;
    typedef T                  value_type;
    typedef uint64_t           size_type;
  protected:
    array_t *xarray_a_;
    size_type xarray_idx_;
  public:
    inline external_array1_item_ref(array_t* a, size_type idx)
        : xarray_a_(a), xarray_idx_(idx) {
    }

    inline external_array1_item_ref(const external_array1_item_const_ref<T>& x) 
        : xarray_a_(x.xarray_a_), xarray_idx_(x.xarray_idx_) {
    }

    inline external_array1_item_ref(const external_array1_item_ref<T>& x) 
        : xarray_a_(x.xarray_a_), xarray_idx_(x.xarray_idx_) {
    }

    inline size_type index() const {
      return xarray_idx_;
    }

    inline value_type value() const {
      return xarray_a_->item(xarray_idx_);
    }

    inline operator value_type() const {
      return value();
    }

    inline self_t &operator=(const value_type& x) {
      //std::cerr << "x[" << xarray_idx_ << "]=" << x << std::endl;
      xarray_a_->put(x, xarray_idx_); return *this;
    }

    inline self_t &operator=(const external_array1_item_const_ref<T>& x) {
      //std::cerr << "x[" << xarray_idx_ << "]= const y[" << x.index() << "]" << std::endl;
      xarray_a_->put(x.value(), xarray_idx_); return *this;
    }

    inline self_t &operator=(const external_array1_item_ref<T>& x) {
      //std::cerr << "x[" << xarray_idx_ << "]= y[" << x.index() << "]" << std::endl;
      xarray_a_->put(x.value(), xarray_idx_); return *this;
    }
  };

  // Iterator on external arrays
  template <class T, bool is_const, bool is_reverse>
  class external_array1_item_iterator {
  public:
    typedef external_array1_item_iterator<T, is_const, is_reverse> self_t;
    typedef std::random_access_iterator_tag  iterator_category;
    typedef T         value_type;
    typedef T*        pointer; // FIXME
    typedef uint64_t  size_type;
    typedef int64_t   difference_type;
    typedef typename gen_if< is_const, const external_array1<T>, external_array1<T> >::type array_t;
    typedef typename gen_if< is_const, external_array1_item_const_ref<T>, external_array1_item_ref<T> >::type reference;
    //protected:
    array_t *xarray_a_;
    size_type xarray_idx_;
  public:
    inline external_array1_item_iterator(array_t* a, size_type idx):
        xarray_a_(a), xarray_idx_(idx) {
    }

  public: // Moving
    
    inline self_t& operator ++() {
      ++xarray_idx_; return *this;
    }
    inline self_t& operator --() {
      --xarray_idx_; return *this;
    }

    SL_OP_INCREMENTABLE(self_t);
    SL_OP_DECREMENTABLE(self_t);
    
    inline self_t& operator +=(difference_type n) {
      xarray_idx_ += n ;
      return (*this);
    }
      
    inline self_t& operator -=(difference_type n) {
      xarray_idx_ -= n ;
      return (*this);
    }
    
    SL_OP_ADDABLE2(self_t, difference_type);
    SL_OP_SUBTRACTABLE2(self_t, difference_type);
    
    inline difference_type operator -(const self_t& other) const {
      return (xarray_idx_ - other.xarray_idx_);
    }
    
  public: // deref

    inline reference operator* () const {
      if (is_reverse) {
        return reference(xarray_a_, xarray_a_->size()-xarray_idx_-1);
      } else {
        return reference(xarray_a_, xarray_idx_);
      }        
    }

    inline reference operator [](difference_type n) const {
      if (is_reverse) {
        return reference(xarray_a_, xarray_a_->size()-xarray_idx_-n-1);
      } else {
        return reference(xarray_a_, xarray_idx_+n);
      }        
    }
    
  public: // comparison

    inline bool operator<(const self_t& other) const {
      return xarray_idx_ < other.xarray_idx_;
    }

    inline bool operator==(const self_t& other) const {
      return xarray_idx_ == other.xarray_idx_;
    }
    
    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);
  };

  /// Arrays maintai ned in external memory
  template <class T>
  class external_array1 {
  public:
    
    typedef external_array1<T>                self_t;

  public: // std::vector compatibility
    
    typedef T                                   value_type;
    typedef T*                                  pointer; // FIXME
    typedef external_array1_item_ref<T>         reference;
    typedef external_array1_item_const_ref<T>   const_reference;
    typedef uint64_t                            size_type;
    typedef int64_t                             difference_type;
    
    typedef external_array1_item_iterator<value_type, false, false> iterator;
    typedef external_array1_item_iterator<value_type, true, false>  const_iterator;
    typedef external_array1_item_iterator<value_type, false, true>  reverse_iterator;
    typedef external_array1_item_iterator<value_type, true, true>   const_reverse_iterator;
    
  protected:
    
    std::string            file_name_;
    os_file::file_handle_t file_descriptor_;
    bool                   file_creatable_;
    bool                   file_writable_;
    bool                   file_resizable_;
    bool                   file_persistent_;
    uint64_t               file_capacity_;        // Number of bytes storable
    
  protected:

    enum { INCORE_CACHE_REGION_COUNT = 4 };     // Max number of memory mapped regions

    size_type        incore_cache_region_size_;    // Multiple of system page size >= element_size_
    mutable uint64_t incore_cache_region_begin_[INCORE_CACHE_REGION_COUNT];   // Current memory mapped region start file offset
    mutable uint64_t incore_cache_region_end_[INCORE_CACHE_REGION_COUNT];     // Current memory mapped region end file offset
    mutable char*    incore_cache_region_address_[INCORE_CACHE_REGION_COUNT]; // Current memory mapped region start address

    mutable size_type stat_page_in_count_;
    mutable uint64_t  stat_page_in_byte_count_;
    mutable size_type stat_page_out_count_;
    mutable uint64_t  stat_page_out_byte_count_;
    
  protected:

    size_type           size_;                 // Number of items
    
  protected:

    inline void init_cache_size(size_type suggested_cache_size) {
      size_type element_ram_size = sizeof(T);
      size_type page_size = os_file::memory_page_size();

      incore_cache_region_size_ =
        ((std::max(suggested_cache_size, element_ram_size) / page_size) +
         (std::max(suggested_cache_size, element_ram_size) % page_size == 0 ? 0 : 1)) * page_size;
    }
    
    inline uint64_t offset(size_type index) const {
      return uint64_t(index)*sizeof(value_type);
    }

    inline value_type* incore_address(size_type index) {
      // Search through current incore cache, compute address and move
      // active cache at first position (LRU policy)
      value_type* result = 0;
      uint64_t o = offset(index);
      for (int i=0; i<INCORE_CACHE_REGION_COUNT; ++i) {
	if (o>=incore_cache_region_begin_[i] && o<incore_cache_region_end_[i]) {
	  result = (value_type*)(incore_cache_region_address_[i] + (o-incore_cache_region_begin_[i]));
	  // Move to first position
	  for (int k=i; k>0; --k) {
	    std::swap(incore_cache_region_begin_[k], incore_cache_region_begin_[k-1]);
	    std::swap(incore_cache_region_end_[k], incore_cache_region_end_[k-1]);
	    std::swap(incore_cache_region_address_[k], incore_cache_region_address_[k-1]);
	  }
	}
      }
      return result;
    }

    inline const value_type* incore_address(size_type index) const {
      return const_cast<self_t*>(this)->incore_address(index);
    }

    inline void page_out(std::size_t k) const {
      assert(k<INCORE_CACHE_REGION_COUNT);
      if (incore_cache_region_address_[k]) {
	// Update stats
	stat_page_out_count_ += 1;
	stat_page_out_byte_count_ += (std::size_t)(incore_cache_region_end_[k] - incore_cache_region_begin_[k]);

	// Perform action
	os_file::memory_unmap(incore_cache_region_address_[k], 
			      (std::size_t)(incore_cache_region_end_[k] - incore_cache_region_begin_[k]));
	for (int i=int(k); i<INCORE_CACHE_REGION_COUNT-1; ++i) {
	  incore_cache_region_address_[i] = incore_cache_region_address_[i+1];
	  incore_cache_region_begin_[i] = incore_cache_region_begin_[i+1];
	  incore_cache_region_end_[i] = incore_cache_region_end_[i+1];
	}
	incore_cache_region_address_[INCORE_CACHE_REGION_COUNT-1] = 0;
	incore_cache_region_begin_[INCORE_CACHE_REGION_COUNT-1] = uint64_t(-1);
	incore_cache_region_end_[INCORE_CACHE_REGION_COUNT-1] = uint64_t(-1);

      }
    }

    inline void page_out() const {
      for (int i=0; i<INCORE_CACHE_REGION_COUNT; ++i) {
	if (incore_cache_region_address_[i]) {
	  os_file::memory_unmap(incore_cache_region_address_[i], 
		  (std::size_t)(incore_cache_region_end_[i] - incore_cache_region_begin_[i]));
	  incore_cache_region_address_[i] = 0;
	  incore_cache_region_begin_[i] = uint64_t(-1);
	  incore_cache_region_end_[i] = uint64_t(-1);
	}
      }
    }

    inline void page_in(size_type index_begin, size_type index_end) const {
      uint64_t begin = offset(index_begin);
      uint64_t end   = offset(index_end);

      bool found = false;
      for (int i=0; i<INCORE_CACHE_REGION_COUNT && !found; ++i) {
	found = (begin>=incore_cache_region_begin_[i] && end<=incore_cache_region_end_[i]);
	if (found) {
	  // Move to first position in LRU
	  for (int k=i; k>0; --k) {
	    std::swap(incore_cache_region_begin_[k], incore_cache_region_begin_[k-1]);
	    std::swap(incore_cache_region_end_[k], incore_cache_region_end_[k-1]);
	    std::swap(incore_cache_region_address_[k], incore_cache_region_address_[k-1]);
	  }
	}
      }

      if (!found) {
	// Region not mapped, map new one
        uint64_t new_region_begin = (begin/incore_cache_region_size_)*incore_cache_region_size_;
	uint64_t new_region_end   = (end/incore_cache_region_size_)*incore_cache_region_size_;
	while (new_region_end < end) new_region_end += incore_cache_region_size_;
        if (!is_writable()) {
          new_region_end = std::min(new_region_end, file_capacity_);
        }
	
	// Check overlaps with existing regions
	for (int i=0; i<INCORE_CACHE_REGION_COUNT; ++i) {
	  if (!(new_region_end-1 < incore_cache_region_begin_[i] ||
		incore_cache_region_end_[i] < new_region_begin+1)) {
	    SL_TRACE_OUT(1) << "Detected overlap at " << i << ": " <<
	      new_region_begin << ".." << new_region_end << " over " <<
	      incore_cache_region_begin_[i] << " " << incore_cache_region_end_[i] <<
	      std::endl;
	    // Overlap!
	    new_region_begin = std::min(new_region_begin, incore_cache_region_begin_[i]);
	    new_region_end   = std::max(new_region_end, incore_cache_region_end_[i]);
	    page_out(i); --i;
	  }
	}
	std::size_t new_region_size = (std::size_t)(new_region_end-new_region_begin);
	char* new_region_address = (char*)os_file::memory_map(new_region_size,
							      is_writable() ? os_file::OS_READ_WRITE : os_file::OS_READ_ONLY,
							      file_descriptor_,
							      new_region_begin);        
	if (!new_region_address) {
	  SL_TRACE_OUT(-1) << "map error: " << 
	    " sz=" << size_t(new_region_size) << 
	    " offset=" << new_region_begin << 
	    std::endl;
	} else {
	  // Insert in cache
	  SL_TRACE_OUT(1) << "Page fault: Range = " << index_begin << " " << index_end << std::endl;
	  page_out(INCORE_CACHE_REGION_COUNT-1);
	  for (int i=INCORE_CACHE_REGION_COUNT-1; i>0; --i) {
	    incore_cache_region_address_[i] = incore_cache_region_address_[i-1];
	    incore_cache_region_begin_[i] = incore_cache_region_begin_[i-1];
	    incore_cache_region_end_[i] = incore_cache_region_end_[i-1];
	  }
	  incore_cache_region_address_[0] = new_region_address;
	  incore_cache_region_begin_[0] = new_region_begin;
	  incore_cache_region_end_[0] = new_region_end;
	  
	  // Update stats
	  stat_page_in_count_ += 1;
	  stat_page_in_byte_count_ += (std::size_t)(incore_cache_region_end_[0] - incore_cache_region_begin_[0]);
	}
      }
    }

    inline void page_in(size_type index_begin) const {
      page_in(index_begin, index_begin+1);
    }

    bool file_is_open() const {
      return file_descriptor_ > 0;
    }
    
    void file_close() {
      page_out();
      if (file_is_open()) {
        if (is_resizable()) {
          os_file::file_resize(file_descriptor_, size_ * sizeof(value_type)); 
        }
        os_file::file_close(file_descriptor_);
        file_descriptor_ = -1;
        if (!is_persistent()) {
	  os_file::file_delete(file_name_.c_str());
        }
      }
    }
  
    void file_open(std::string fname,
                   bool creatable = true,
                   bool writable = true,
                   bool resizable = true,
                   bool persistent = true) {
      file_close();

      file_name_ = fname;
      file_creatable_ = creatable;
      file_writable_ = writable;
      file_resizable_ = resizable;
      file_persistent_ = persistent;
      file_descriptor_ = os_file::file_open(fname.c_str(),
                                            writable  ?  os_file::OS_READ_WRITE : os_file::OS_READ_ONLY,
                                            creatable ?  os_file::OS_OPEN_CREATE_IF_NOT_PRESENT : os_file::OS_OPEN_EXISTING);
      file_capacity_ = os_file::file_size(file_descriptor_);
      size_ = file_capacity_/sizeof(value_type);

      SL_TRACE_OUT(1) << "[" << file_name_ << "]" << "-> fd:" << file_descriptor_ << " ->capacity: " << file_capacity_ << std::endl;

      // Clear paging stats
      stat_clear();
    }
    
    inline void file_resize(size_type new_size) {
      SL_REQUIRE("Resizable", is_resizable());
      SL_TRACE_OUT(1) << "RESIZE(" << new_size << ") BEGIN" << std::endl;
      if (new_size < size_) {
	uint64_t new_capacity = (((new_size * sizeof(value_type))/incore_cache_region_size_+
                                  ((new_size * sizeof(value_type))%incore_cache_region_size_ == 0 ? 0 : 1)))*incore_cache_region_size_; // FIXME parameterize
        if (2*new_capacity < file_capacity_) { // FIXME parameterize
          SL_TRACE_OUT(1) << "SHRINK" << std::endl;
          page_out();
	  file_capacity_ = new_capacity;
          os_file::file_resize(file_descriptor_, new_capacity);
          SL_TRACE_OUT(1) << "CAPACITY Shrinked to " << file_capacity_ << std::endl;
        }
        size_ = new_size;
      } else if (new_size > size_) {
        if (new_size*sizeof(value_type) > size_type(file_capacity_)) {
          SL_TRACE_OUT(1) << "GROW" << std::endl;
          page_out();
	  file_capacity_ = (new_size * sizeof(value_type)/incore_cache_region_size_+1)*incore_cache_region_size_; // FIXME parameterize
          os_file::file_resize(file_descriptor_, file_capacity_);
        }
        size_ = new_size;
      }
      SL_TRACE_OUT(1) << "RESIZE END" << std::endl;
    }

  private:
    
    external_array1(const self_t& other); // Not implemented!
    self_t& operator=(const self_t& other); // Not implemented!
    
  public:

    external_array1(std::string fname,
                    bool creatable = true,
                    bool writable = true,
                    bool resizable = true,
                    bool persistent = false,
                    bool emptyatcreation = true,
		    size_type suggested_cache_size = 1024*1024) {
      file_name_ = "";
      file_descriptor_ = -1;
      file_creatable_ = false;
      file_writable_ = false;
      file_resizable_ = false;
      file_persistent_ = false;
      file_capacity_ = 0;
      incore_cache_region_size_ = 0;
      for (int i=0; i<INCORE_CACHE_REGION_COUNT; ++i) {
	incore_cache_region_address_[i] = 0;
	incore_cache_region_begin_[i] = uint64_t(-1);
	incore_cache_region_end_[i] = uint64_t(-1);
      }
      size_ = 0;

      init_cache_size(suggested_cache_size);
      file_open(fname, creatable, writable, resizable, persistent);
      if (emptyatcreation && is_open()) clear();
    }

    /**
     *  Mode: "r" readable, persistent
     *  Mode: "w" readable+writable, truncate file to zero length or create file for writing, persistent
     *  Mode: "a" readable+writable, do not truncate, but create file for writing, persistent
     *  Mode: "t" temporary read/write file, the file is created if it does not exist,
     *            otherwise it is truncated.
     */
    external_array1(std::string fname,
                    const char* mode,
		    size_type suggested_cache_size = 1024*1024) {
      file_name_ = "";
      file_descriptor_ = -1;
      file_creatable_ = false;
      file_writable_ = false;
      file_resizable_ = false;
      file_persistent_ = false;
      file_capacity_ = 0;
      incore_cache_region_size_ = 0;
      for (int i=0; i<INCORE_CACHE_REGION_COUNT; ++i) {
	incore_cache_region_address_[i] = 0;
	incore_cache_region_begin_[i] = uint64_t(-1);
	incore_cache_region_end_[i] = uint64_t(-1);
      }
      size_ = 0;
      bool creatable = false;
      bool writable = false;
      bool resizable = false;
      bool persistent = false;
      bool emptyatcreation = true;
      if (mode && mode[0] == 'r' && mode[1] == '\0') { creatable=false; writable=false; resizable=false; persistent=true; emptyatcreation=false; }
      if (mode && mode[0] == 'w' && mode[1] == '\0') { creatable=true;  writable=true;  resizable=true;  persistent=true; emptyatcreation=true; }
      if (mode && mode[0] == 'a' && mode[1] == '\0') { creatable=true;  writable=true;  resizable=true;  persistent=true; emptyatcreation=false; }
      if (mode && mode[0] == 't' && mode[1] == '\0') { creatable=true;  writable=true;  resizable=true;  persistent=false; emptyatcreation=true; }

      init_cache_size(suggested_cache_size);
      file_open(fname, creatable, writable, resizable, persistent);
      if (emptyatcreation) clear();
    }
    
    inline ~external_array1() {
      file_close();
    }

    /// Move everything out-of-core
    inline void minimize_footprint() const {
      page_out();
    }

    /// Move everything in-core (single mmap for entire range)
    inline void maximize_footprint() const {
      size_type sz = size();
      if (sz) {
	page_in(size_type(0), sz);
      }
    }
    
    inline bool is_open() const {
      return file_is_open();
    }
    inline bool is_creatable() const {
      return file_creatable_;
    }
    inline bool is_writable() const {
      return file_writable_;
    }
    inline bool is_resizable() const {
      return file_resizable_;
    }
    inline bool is_persistent() const {
      return file_persistent_;
    }

    const std::string& file_name() const {
      return file_name_;
    }
    
    inline void close() {
      file_close();
    }

    inline void reopen() {
      if (!is_open()) {
        file_open(file_name(), is_creatable(), is_writable(), is_resizable(), is_persistent());
      }
      // Clear stats
      stat_clear();
    }

    // Clear paging statistics - usually done at file open
    inline void stat_clear() const {
      stat_page_in_count_ = 0;
      stat_page_in_byte_count_ = 0;
      stat_page_out_count_ = 0;
      stat_page_out_byte_count_ = 0;
    }
    
    // Number of mmaps since open
    inline size_type stat_page_in_count() const { return stat_page_in_count_; }

    // Number of mmapped bytes since open
    inline uint64_t stat_page_in_byte_count() const { return stat_page_in_byte_count_; }

    // Number of munmaps since open
    inline size_type stat_page_out_count() const { return stat_page_out_count_; }

    // Number of munmapped bytes since open
    inline uint64_t stat_page_out_byte_count() const { return stat_page_out_byte_count_; }
    
    /**
     *  Ensure all data in range index_lo..index_hi is incore, return pointer to it.
     *  The pointer remains valid only as long as no other call to external array
     *  is made.
     */
    inline value_type* range_page_in(size_type index_lo, size_type index_hi) {
      page_in(index_lo, index_hi+1);
      return incore_address(index_lo);
    }

    /**
     *  Ensure all data in range index_lo..index_hi is incore, return pointer to it.
     *  The pointer remains valid only as long as no other call to external array
     *  is made.
     */
    inline const value_type* range_page_in(size_type index_lo, size_type index_hi) const {
      page_in(index_lo, index_hi+1);
      return incore_address(index_lo);
    }
    
    inline iterator begin() {
      return iterator(this, 0);
    }

    inline iterator end() {
      return iterator(this, size());
    }

    inline const_iterator begin() const {
      return const_iterator(this, 0);
    }

    inline const_iterator end() const {
      return const_iterator(this, size());
    }

    inline reverse_iterator rbegin() {
      return reverse_iterator(this, 0);
    }

    inline reverse_iterator rend() {
      return iterator(this, size());
    }

    inline const_reverse_iterator rbegin() const {
      return const_reverse_iterator(this, 0);
    }

    inline const_reverse_iterator rend() const {
      return const_reverse_iterator(this, size());
    }
    
    inline bool good_index(const size_type index) const {
      return (index<size());
    }

    inline size_type size() const {
      return size_;
    }

    inline size_type max_size() const {
      return ~size_type(0);
    }
    
    inline size_type capacity() const {
      return file_capacity_;
    }

    inline bool empty() const {
      return size() == 0;
    }

    inline void put(const value_type& x, size_type index) {
      SL_REQUIRE("Good index", good_index(index));
      SL_REQUIRE("Writable", is_writable());
      page_in(index);
      (*incore_address(index)) = x;
    }

    inline value_type item(size_type index) const {
      SL_REQUIRE("Good index", good_index(index));
      page_in(index);
      return *incore_address(index);
    }
    
    inline reference operator[](size_type index) {
      SL_REQUIRE("Good index", good_index(index));
      // No - can be called non const SL_REQUIRE("Writable", is_writable());
      return reference(this, index);
    }

    inline value_type operator[](size_type index) const {
      SL_REQUIRE("Good index", good_index(index));
      page_in(index);
      return *incore_address(index);
    }

    inline self_t& operator=(const self_t& other) const {
      SL_REQUIRE("Writable", is_writable());
      clear();
      for (const_iterator it = other.begin(); it!= other.end(); ++it) {
        push_back(*it);
      }
    }

    inline void reserve(size_type sz) {
      if (sz > capacity()) {
        size_type old_size = size(); 
        file_resize(sz);
        size_ = old_size;
      }
    }

    inline reference front() {
      SL_REQUIRE("Not empty", !empty());
      return reference(this, 0);
    }

    inline reference back() {
      SL_REQUIRE("Not empty", !empty());
      return reference(this, size()-1);
    }

    inline const_reference front() const {
      SL_REQUIRE("Not empty", !empty());
      return const_reference(this, 0);
    }

    inline const_reference back() const {
      SL_REQUIRE("Not empty", !empty());
      return const_reference(this, size()-1);
    }

    inline void push_back(const value_type& x) {
      SL_REQUIRE("Writable", is_writable());
      SL_REQUIRE("Resizable", is_resizable());
      resize(size()+1);
      this->operator[](size()-1) = x;
    }

    inline void pop_back() {
      SL_REQUIRE("Writable", is_writable());
      SL_REQUIRE("Resizable", is_resizable());
      SL_REQUIRE("Not empty", !empty());
      resize(size()-1);
    }
    
    void swap(self_t& other) {
      std::swap(file_name_, other.file_name_);
      std::swap(file_descriptor_, other.file_descriptor_);
      std::swap(file_creatable_, other.file_creatable_);
      std::swap(file_writable_, other.file_writable_);
      std::swap(file_resizable_, other.file_resizable_);
      std::swap(file_persistent_, other.file_persistent_);
      std::swap(file_capacity_, other.file_capacity_);
      std::swap(incore_cache_region_size_, other.incore_cache_region_size_);
      for (int i=0; i<INCORE_CACHE_REGION_COUNT; ++i) {
	std::swap(incore_cache_region_begin_[i], other.incore_cache_region_begin_[i]);
	std::swap(incore_cache_region_end_[i], other.incore_cache_region_end_[i]);
	std::swap(incore_cache_region_address_[i], other.incore_cache_region_address_[i]);
      }
      std::swap(size_, other.size_);
    }

    iterator insert(iterator pos, const value_type& x) {
      if (pos == end()) {
        // Insert at end
        push_back(x);
      } else {
        // Insert middle
        resize(size()+1);
        value_type x_copy = x;
        std::copy_backward(pos, end()-2, end()-1);
        *pos = x_copy;
      }
      return begin()+pos;
    }

    template <class InputIterator>
    void insert(iterator pos,
                InputIterator f, InputIterator l) {
      if (pos == end()) {
        for (InputIterator it=f; it!=l; ++it) {
          push_back(*it);
        }
      } else {
        size_type n = l-f;
        resize(size()+n);
        std::copy_backward(pos, end()-n-1, end()-n);
        for (InputIterator it=f; it!=l; ++it) {
          *pos = *it; ++pos;
        }
      }
    }

    iterator erase(iterator pos) {
      if (pos+1 != end()) {
        std::copy(pos+1, end(), pos);
      }
      resize(size()-1);
      return pos;
    }

    iterator erase(iterator first, iterator last) {
      iterator i(std::copy(last, end(), first));
      resize(size() - (last-first));
      return first;
    }
    
    inline void clear() {
      SL_REQUIRE("Writable", is_writable());
      SL_REQUIRE("Resizable", is_resizable());
      resize(0);
    }
    
    inline void resize(size_type new_size) {
      SL_REQUIRE("Writable", is_writable());
      SL_REQUIRE("Resizable", is_writable());
      file_resize(new_size);
    }

    inline bool operator==(const self_t& other) const {
      if (size() != other.size()) {
        return false;
      } else {
        bool result = true;
        for (size_type i=0; i<size() && result; ++i) {
          result = (*this)[i] == other[i];
        }
        return result;
      }
    }

    inline bool operator<(const self_t& other) {
      if (size() >= other.size()) {
        return false;
      } else {
        for (size_type i=0; i<size(); ++i) {
          if ((*this)[i] < other[i]) {
            return true;
          } else if ((*this)[i] > other[i]) {
            return false;
          }
        }
        return false;
      }
    }

    SL_OP_COMPARABLE1(self_t);
    SL_OP_EQUALITY_COMPARABLE1(self_t);

  };
}

// Redefinition of std::sort to use sl::quicksort
namespace std {

  template <typename T, bool is_const, bool is_reverse, typename BP>
  void sort(const sl::external_array1_item_iterator<T,is_const,is_reverse>& begin,
	    const sl::external_array1_item_iterator<T,is_const,is_reverse>& end,
	    const BP& is_less) {
    return sl::quicksort(begin, end, is_less);
  }

  template <typename T, bool is_const, bool is_reverse>
  void sort(const sl::external_array1_item_iterator<T,is_const,is_reverse>& begin,
	    const sl::external_array1_item_iterator<T,is_const,is_reverse>& end) {
    return sl::quicksort(begin, end);
  }

}

#endif
