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
#ifndef SL_MEMORY_POOL_HPP
#define SL_MEMORY_POOL_HPP

#include <sl/config.hpp>
#include <functional> // std::less, std::less_equal, std::greater
#include <new>        // new[], delete[], std::nothrow
#include <exception>  // std::invalid_argument
#include <set>        // std::set
#include <iostream>   // std::ostream

namespace sl {

  /**
   *  A fast memory allocator for blocks of memory of a
   *  single size.
   */
  class memory_pool {
  private: // Noncopyable
    memory_pool(memory_pool const&); 
    memory_pool& operator=(memory_pool const&);
  protected:

    /// The beginning of a free list element
    struct pool_element {
      pool_element* next_;
    };

    /// The beginning of a block of allocated pool memory
    struct pool_block {
      std::size_t   element_count_; //< the size of the block
      std::size_t   free_count_;    //< the number of elements in the free list
      pool_element* free_elements_; //< the list of free elements within the block
    };

  protected:
    
    const std::string     name_;          //< The name of the pool

    std::set<pool_block*> free_blocks_;   //< The set of blocks with free elements
    std::set<pool_block*> full_blocks_;   //< The set of blocks without free elements

    std::size_t block_header_size_;       //< The size of the header of a pool block

    std::size_t element_requested_size_;  //< The requested element size
    std::size_t element_actual_size_;     //< The allocated element size

    std::size_t chunk_first_count_;       //< The number of elements in the first block
    std::size_t chunk_grow_factor_;       //< The grow factor
    std::size_t chunk_current_count_;     //< The number of elements allocated in the next block

    std::size_t capacity_;                //< The number of free and allocated elemements in all blocks
    std::size_t allocated_count_;         //< The number of allocated elements

  protected: // Helpers

    /**
     *  Initialize the actual element size and block header size
     *  for this pool, given that the user wants elements of
     *  size element_requested_size_
     */
    void initialize() {
      SL_REQUIRE("Good firsr count", chunk_first_count_ > 0);
      chunk_current_count_ = chunk_first_count_;
      element_actual_size_ = lcm(element_requested_size_, sizeof(pool_element));
      block_header_size_   = lcm(element_actual_size_,    sizeof(pool_block));
    }

    /// Allocate b bytes from the system
    static inline void* sys_malloc(const std::size_t n) { 
#if 0
      return (void*)new (std::nothrow) char[n]; 
#else
      return (void*)new char[n]; 
#endif
    }

    /// Return the block pointed by ptr to the system
    static inline void sys_free(void* const ptr) { 
      delete [] (char*)ptr; 
    }

    /// Request a new block of memory from the system
    void grab_memory() {
      // Grab memory from system
      pool_block* new_block = static_cast<pool_block*>(sys_malloc(block_header_size_ + 
								  chunk_current_count_ * element_actual_size_));
      new_block->element_count_ = chunk_current_count_;
      new_block->free_count_    = chunk_current_count_;
      new_block->free_elements_ = NULL;

      // Populate free list
      char* new_elements = reinterpret_cast<char*>(new_block)+block_header_size_;
      for (size_t i=0; i<new_block->element_count_; ++i) {
	pool_element* elt = reinterpret_cast<pool_element*>(new_elements + i * element_actual_size_);
	elt->next_ = new_block->free_elements_;
	new_block->free_elements_ = elt;
      }

      // Insert new block in set of blocks with free elements
      free_blocks_.insert(new_block);
      capacity_ += new_block->element_count_;

      // Compute next block size
      chunk_current_count_ *= chunk_grow_factor_;
    }

    /// The pool block containing element e in the set of blocks b or b.end() if none
    std::set<pool_block*>::iterator find_block(const std::set<pool_block*>& b_param,
					       const pool_element* e) const {
      std::set<pool_block*>& b = const_cast< std::set<pool_block*>& >(b_param);
      std::set<pool_block*>::iterator result = b.lower_bound(reinterpret_cast<pool_block*>(const_cast<pool_element*>(e)));
      if (result == b.begin()) {
	// Address below first block
	result = b.end();
      } else {
	if (b.size() > 0) --result;
	if (result != b.end()) {
	  // Range check
	  const char* elt_bgn = reinterpret_cast<const char*>(*result)+block_header_size_;
	  const char* elt_end = elt_bgn + element_actual_size_ * (*result)->element_count_;
	  const char* elt_e   = reinterpret_cast<const char*>(e);

	  if (!(elt_e >= elt_bgn && elt_e < elt_end)) {
	    result = b.end();
	  }
	}
      }
      return result;
    }

  public:

    /**
     *  Create a pool for blocks of memory of size requested_size. 
     *  The first allocation will allocate chunk_first_count element,
     *  the second one chunk_grow_factor*chunk_first_count, and so on.
     *  If chunk_first_count is set to zero, a default value is
     *  automatically computed
     */
    explicit memory_pool(const std::string& pool_name,
				const std::size_t requested_size,
				const std::size_t chunk_first_count = 0,
				const std::size_t chunk_grow_factor = 1)
      :
      name_(pool_name),
      block_header_size_(0),
      element_requested_size_(requested_size),
      element_actual_size_(0),
      chunk_first_count_((chunk_first_count == 0) ? 
			 (sl::max(std::size_t(16), std::size_t(65536)/sl::max(std::size_t(1),requested_size))) :
			 (sl::max(std::size_t(1), chunk_first_count))),
      chunk_grow_factor_(chunk_grow_factor),
      chunk_current_count_(0),
      capacity_(0),
      allocated_count_(0) { 

      SL_REQUIRE("Positive request", requested_size>0);
      SL_REQUIRE("Non null grow factor", chunk_grow_factor>=1);

      initialize();
    }

    /// Delete the memory pool, returning all memory to the system
    ~memory_pool() { 
      purge(); 
    }

    /// The name of the pool
    const std::string& name() const {
      return name_;
    }

    /// The requested size of an element
    inline std::size_t element_requested_size() const {
      return element_requested_size_;
    }

    /// The actual size of an element. May be bigger than the requested size because of alignement requirements
    inline std::size_t element_actual_size() const {
      return element_actual_size_;
    }

    /// The number of elements that may be allocated without requesting memory from the system
    inline std::size_t capacity() const {
      return capacity_; 
    }

    /// The number of elements currently allocated
    inline std::size_t allocated_count() const {
      return allocated_count_;
    }

    /// The number of elements currently in the free list
    inline std::size_t free_count() const {
      return capacity_ - allocated_count_;
    }

    /// Is e an element from the current pool?
    inline bool is_pool_element(const void *e) const {
      return 
	(find_block(full_blocks_,
		    reinterpret_cast<const pool_element*>(e)) != const_cast< std::set<pool_block*>& >(full_blocks_).end()) ||
	(find_block(free_blocks_,
		    reinterpret_cast<const pool_element*>(e)) != const_cast< std::set<pool_block*>& >(free_blocks_).end());
    }

    /// A new block of memory of size element_size
    void *allocate() {
      if (free_count() == 0) {
	grab_memory();
      }
      SL_CHECK("Has free memory", free_count() > 0);
      std::set<pool_block*>::iterator b_it = free_blocks_.begin();
      SL_CHECK("Has free block", b_it != free_blocks_.end());
      pool_element *result = (*b_it)->free_elements_;
      (*b_it)->free_elements_ = result->next_;
      --(*b_it)->free_count_;
      if ((*b_it)->free_count_ == 0) {
	full_blocks_.insert(*b_it);
	free_blocks_.erase(b_it);
      }
      ++allocated_count_;
      return result;
    }

    /// Release the block of memory pointed by elt_ptr
    void release(void* elt_ptr) {
      SL_REQUIRE("Null or pool element", is_pool_element(elt_ptr));
      SL_REQUIRE("Allocated", allocated_count() != 0);
      pool_element* e = static_cast<pool_element*>(elt_ptr);

      std::set<pool_block*>::iterator b_it = find_block(full_blocks_, e);
      if (b_it == full_blocks_.end()) {
	b_it = find_block(free_blocks_, e);
	SL_CHECK("Has block", b_it != free_blocks_.end());
	// Deallocate from partially free block
	e->next_ = (*b_it)->free_elements_;
	(*b_it)->free_elements_ = e;
	++(*b_it)->free_count_;
	if ((*b_it)->free_count_ == (*b_it)->element_count_) {
	  // Release empty block
	  capacity_ -= (*b_it)->element_count_;
	  sys_free(*b_it);
	  free_blocks_.erase(b_it);
	}
      } else {
	SL_CHECK("Has block", b_it != full_blocks_.end());
	// Deallocate from full block, move block to free list
	e->next_ = (*b_it)->free_elements_;
	(*b_it)->free_elements_ = e;
	++(*b_it)->free_count_;
	free_blocks_.insert(*b_it);
	full_blocks_.erase(b_it);
      }
      --allocated_count_;
    }
      
    // Releases *all* memory blocks, even if chunks are still allocated
    void purge() {
      if (allocated_count() > 0) {
	SL_TRACE_OUT(0) << "Purging memory pool with " << allocated_count() << " live elements!" << std::endl;
      }
      std::set<pool_block*>::iterator b_it = full_blocks_.begin();
      while (b_it != full_blocks_.end()) {
	capacity_ -= (*b_it)->element_count_;
	allocated_count_ -= ((*b_it)->element_count_ - (*b_it)->free_count_);
	sys_free(*b_it);
	full_blocks_.erase(b_it);
      }
      b_it = free_blocks_.begin();
      while (b_it != free_blocks_.end()) {
	capacity_ -= (*b_it)->element_count_;
	allocated_count_ -= ((*b_it)->element_count_ - (*b_it)->free_count_);
	sys_free(*b_it);
	free_blocks_.erase(b_it);
      } 
      SL_CHECK("Null capacity", capacity_ == 0);
      SL_CHECK("Null allocated count", allocated_count_ == 0);
      chunk_current_count_ = chunk_first_count_;
    }

  }; // memory_pool

} // namespace sl

inline std::ostream& operator <<(std::ostream& s, 
				 const sl::memory_pool& p) {
  s << "Memory Pool [" << p.name() << "]: " <<
    (p.allocated_count() * p.element_actual_size()) / 1024 <<
    "/" <<
    (p.capacity() * p.element_actual_size()) / 1024 <<
    " kB used" <<
    " (" <<
    p.allocated_count() <<
    " live instances)." <<
    std::endl;
  return s;
}
    
#define SL_DECLARE_CLASS_POOL_ALLOCATION_FEATURES(THIS_T_, POOL_VARIABLE_NAME_) \
  protected: \
    static sl::memory_pool* POOL_VARIABLE_NAME_; \
  public: \
    inline static sl::memory_pool& instance_memory_pool() { \
      if (!POOL_VARIABLE_NAME_) {\
	POOL_VARIABLE_NAME_ = new sl::memory_pool(std::string(#THIS_T_),sizeof(THIS_T_));\
      }\
      return *POOL_VARIABLE_NAME_;\
    }\
    inline void* operator new(std::size_t sz) { \
      SL_REQUIRE("Good size", sz == sizeof(THIS_T_));\
      if (sz); \
      return instance_memory_pool().allocate();\
    }\
    inline void operator delete(void* ptr, std::size_t sz) {\
      SL_REQUIRE("Has pool", (ptr == NULL) || (POOL_VARIABLE_NAME_));\
      SL_REQUIRE("Good ptr", (ptr == NULL) || (POOL_VARIABLE_NAME_->is_pool_element(ptr)));\
      SL_REQUIRE("Good size", (ptr == NULL) || (sz == POOL_VARIABLE_NAME_->element_requested_size()));\
      if (sz); \
      if (ptr) POOL_VARIABLE_NAME_->release(ptr);\
    }

#define SL_DECLARE_CLAS_POOL_STORAGE(THIS_T_, POOL_VARIABLE_NAME_) \
  THIS_T_ :: POOL_VARIABLE_NAME_ = NULL



#endif

