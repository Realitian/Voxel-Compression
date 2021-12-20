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
#ifndef ARENA_ALLOCATOR_HPP
#define ARENA_ALLOCATOR_HPP

#include <sl/config.hpp>
#include <cstddef>
#include <new> // std::bad_alloc
#include <limits>
#include <cassert>
  
namespace sl {
  
  /**
   *  A simple memory arena allocating memory with new/delete
   */
  class default_arena {
  public:
    /// Allocate bytes bytes
    void* allocate(std::size_t bytes) {
      return ::operator new(bytes);
    }

    /// Deallocate memory previously allocated from this arena
    void deallocate(void* ptr, std::size_t bytes) {
      ::operator delete(ptr);
    }

    /// Capacity
    std::size_t capacity() const {
      return std::numeric_limits<std::size_t>::max();
    }
  };

  /**
   * A memory arena optimized for allocations of locally used memory.
   * A local buffer of size N is used for the first allocations, and no memory
   * is released in that case. Upon overflow, memory is allocated using the heap.
   * This is useful for small temporary structures that occasionally overflow (e.g.,
   * small priority queues).
   */
  template <std::size_t N>
  class tmp_arena {
  protected:
    typedef union max_align {
      int                 i     ;
      long                l     ;
      long long           ll    ;
      long double         ld    ;
      double              d     ;
      void*               p     ;
      void (*             pf)() ;
      max_align*          ps    ;
    } max_align_t;
    
    static const std::size_t alignment = sizeof(max_align); //align??
    
    typedef union {
      max_align_t align_;
      char buf_[N];
    } aligned_buf_t;
  
    aligned_buf_t aligned_;
    char* ptr_;
  
    static inline std::size_t align_up(std::size_t n) {
      return (n + (alignment-1)) & (~(alignment-1));
    }
  
    inline bool is_pointer_in_buffer(void* p) const {
      return aligned_.buf_ <= p && p <= aligned_.buf_ + N;
    }
    
  private: // No copy
    tmp_arena(const tmp_arena&);
    tmp_arena& operator=(const tmp_arena&);
  public:
    inline tmp_arena() : ptr_(aligned_.buf_) {}
    inline ~tmp_arena() { ptr_ = 0; }
    
    void* allocate(std::size_t n) {
      assert(is_pointer_in_buffer(ptr_));
      n = align_up(n);
      if (aligned_.buf_ + N >= ptr_ + n) {
	// We're allocating within the arena, return ptr
	char* r = ptr_;
	ptr_ += n;
	return static_cast<void*>(r);
      } else {
	// revert to standard allocation
	return static_cast<void*>(::operator new(n));
      }
    }
  
    inline void deallocate(void* p, std::size_t n) {
      if (is_pointer_in_buffer(p)) {
	// Do nothing -- we could actually check for stack deallocation!
      } else {
	::operator delete(p);
      }
    }
  
    inline std::size_t capacity() const {
      return std::numeric_limits<std::size_t>::max();
    }

    static inline std::size_t local_memory_size() { return N; }
  
    inline std::size_t used_local_memory_size () const {
      return static_cast<std::size_t>(ptr_ - aligned_.buf_);
    }
    
    inline void reset() {
      ptr_ = aligned_.buf_;
    }
  };

} // namespace sl

namespace sl {  

  /**
   *  An instance stateful allocator 
   */
  template <typename ArenaType, typename T>
  class arena_allocator;
    
  template <typename ArenaType>
  class arena_allocator<ArenaType, void> {
  public:
    // arena_allocator specific interface
    typedef ArenaType arena_type;
    
  public:
    explicit arena_allocator(arena_type& arena) throw ()
      :   arena_(arena)
    {}

    const arena_type& arena() const { return arena_; }
    arena_type& arena() { return arena_; }
    
  public:
    // ISO/IEC 14882:2003 20.1.5 [lib.allocator.requirements]
    typedef void value_type;
    typedef void* pointer;
    typedef const void* const_pointer;
    
    template<class Other>
    struct rebind
    {
      typedef arena_allocator<arena_type, Other> other;
    };
    
    template <class Other>
    arena_allocator(
		    const arena_allocator<arena_type, Other>& other) throw ()
      :   arena_(other.arena())
    {
    }
    
  private:
    arena_type& arena_;
  };


  template <typename ArenaType, typename T>
  class arena_allocator {
  public:
    // arena_allocator specific interface
    typedef ArenaType arena_type;

  protected:
    arena_type& arena_;

  public:
    explicit arena_allocator(arena_type& arena) throw () : arena_(arena) {
    }

    const arena_type& arena() const { return arena_; }
    arena_type& arena() { return arena_; }
    
  public:
    // ISO/IEC 14882:2003 20.1.5 [lib.allocator.requirements]
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef T *pointer;
    typedef const T *const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
    
    template<class Other>
    struct rebind
    {
      typedef arena_allocator<arena_type, Other> other;
    };
    
    template <class Other>
    arena_allocator(
		    const arena_allocator<arena_type, Other>& other) throw ()
      : arena_(const_cast<arena_allocator<arena_type, Other>& >(other).arena()) {
    }

    pointer address(reference val) const throw () {
      return addressof(val);
    }
    
    const_pointer address(const_reference val) const throw () {
      return addressof(val);
    }

    void construct(pointer ptr, const_reference val) {
      new (ptr) value_type(val);
    }

    void destroy(pointer ptr) {
      ptr->~value_type();
    }
            
    pointer allocate(size_type n, void* hint = 0) {
      return reinterpret_cast<pointer>(arena_.allocate(n * sizeof(value_type)));
    }
            
    void deallocate(pointer ptr, size_type n) throw () {
      arena_.deallocate(ptr, n * sizeof(value_type));
    }
            
    size_type max_size() const throw() {
      return arena_.capacity() / sizeof(value_type);
    }
  };

  // ISO/IEC 14882:2003 20.1.5 [lib.allocator.requirements]
  template <
    typename ArenaType,
    typename T,
    typename U>
  bool operator==(
		  const arena_allocator<ArenaType, T>& lhs,
		  const arena_allocator<ArenaType, U>& rhs) {
    return
      addressof(lhs.arena()) == addressof(rhs.arena());
  }
  
  template <
    typename ArenaType,
    typename T,
    typename U>
  bool operator!=(
		  const arena_allocator<ArenaType, T>& lhs,
		  const arena_allocator<ArenaType, U>& rhs)
  {
    return !(lhs == rhs);
  }

} // namespace sl

#endif
