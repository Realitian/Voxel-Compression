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
/*===========================================================================
  This library is released under the MIT license. See FSBAllocator.html
  for further information and documentation.

  Copyright (c) 2008-2011 Juha Nieminen

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
  =============================================================================*/

#ifndef _SL_FSB_ALLOCATOR_HPP
#define _SL_FSB_ALLOCATOR_HPP

#include <new>
#include <cassert>
#include <vector>

namespace sl {

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_BOOST
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OPENMP
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_PTHREAD
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC
#define FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
#include <boost/thread.hpp>
typedef boost::mutex fsb_allocator_mutex;
#endif

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OPENMP
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_BOOST
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_PTHREAD
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC
#define FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
#include <omp.h>

  class fsb_allocator_mutex
  {
    omp_lock_t mutex;

  public:
    fsb_allocator_mutex() { omp_init_lock(&mutex); }
    ~fsb_allocator_mutex() { omp_destroy_lock(&mutex); }
    void lock() { omp_set_lock(&mutex); }
    void unlock() { omp_unset_lock(&mutex); }
  };
#endif

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_PTHREAD
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_BOOST
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OPENMP
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC
#define FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
#include <pthread.h>

  class fsb_allocator_mutex
  {
    pthread_mutex_t mutex;

  public:
    fsb_allocator_mutex() { pthread_mutex_init(&mutex, NULL); }
    ~fsb_allocator_mutex() { pthread_mutex_destroy(&mutex); }
    void lock() { pthread_mutex_lock(&mutex); }
    void unlock() { pthread_mutex_unlock(&mutex); }
  };
#endif

#if defined(FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC) || defined(FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC_WITH_SCHED)
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_BOOST
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OPENMP
#undef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_PTHREAD
#define FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC_WITH_SCHED
#include <sched.h>
#endif
  class fsb_allocator_mutex
  {
    volatile int lockFlag;

  public:
    fsb_allocator_mutex(): lockFlag(0) {}
    void lock()
    {
      while(!__sync_bool_compare_and_swap(&lockFlag, 0, 1))
        {
#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_GCC_WITH_SCHED
	  sched_yield();
#endif
        }
    }
    void unlock() { lockFlag = 0; }
  };
#endif

  template<unsigned ElemSize>
  class fsb_allocator_elem_allocator
  {
    typedef std::size_t data_t;
    static const data_t BlockElements = 512;

    static const data_t DSize = sizeof(data_t);
    static const data_t ElemSizeInDSize = (ElemSize + (DSize-1)) / DSize;
    static const data_t UnitSizeInDSize = ElemSizeInDSize + 1;
    static const data_t BlockSize = BlockElements*UnitSizeInDSize;

    class MemBlock
    {
      data_t* block;
      data_t firstFreeUnitIndex, allocatedElementsAmount, endIndex;

    public:
      MemBlock():
	block(0),
	firstFreeUnitIndex(data_t(-1)),
	allocatedElementsAmount(0)
      {}

      bool isFull() const
      {
	return allocatedElementsAmount == BlockElements;
      }

      void clear()
      {
	delete[] block;
	block = 0;
	firstFreeUnitIndex = data_t(-1);
      }

      void* allocate(data_t vectorIndex)
      {
	if(firstFreeUnitIndex == data_t(-1))
	  {
	    if(!block)
	      {
		block = new data_t[BlockSize];
		if(!block) return 0;
		endIndex = 0;
	      }

	    data_t* retval = block + endIndex;
	    endIndex += UnitSizeInDSize;
	    retval[ElemSizeInDSize] = vectorIndex;
	    ++allocatedElementsAmount;
	    return retval;
	  }
	else
	  {
	    data_t* retval = block + firstFreeUnitIndex;
	    firstFreeUnitIndex = *retval;
	    ++allocatedElementsAmount;
	    return retval;
	  }
      }

      void deallocate(data_t* ptr)
      {
	*ptr = firstFreeUnitIndex;
	firstFreeUnitIndex = ptr - block;

	if(--allocatedElementsAmount == 0)
	  clear();
      }
    };

    struct BlocksVector
    {
      std::vector<MemBlock> data;

      BlocksVector() { data.reserve(1024); }

      ~BlocksVector()
      {
	for(std::size_t i = 0; i < data.size(); ++i)
	  data[i].clear();
      }
    };

    static BlocksVector blocksVector;
    static std::vector<data_t> blocksWithFree;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
    static fsb_allocator_mutex mutex;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_BOOST
    struct Lock: boost::mutex::scoped_lock
    {
      Lock(): boost::mutex::scoped_lock(mutex) {}
    };
#else
    struct Lock
    {
      Lock() { mutex.lock(); }
      ~Lock() { mutex.unlock(); }
    };
#endif
#endif

  public:
    static void* allocate()
    {
#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
      Lock lock;
#endif

      if(blocksWithFree.empty())
        {
	  blocksWithFree.push_back(blocksVector.data.size());
	  blocksVector.data.push_back(MemBlock());
        }

      const data_t index = blocksWithFree.back();
      MemBlock& block = blocksVector.data[index];
      void* retval = block.allocate(index);

      if(block.isFull())
	blocksWithFree.pop_back();

      return retval;
    }

    static void deallocate(void* ptr)
    {
      if(!ptr) return;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
      Lock lock;
#endif

      data_t* unitPtr = (data_t*)ptr;
      const data_t blockIndex = unitPtr[ElemSizeInDSize];
      MemBlock& block = blocksVector.data[blockIndex];

      if(block.isFull())
	blocksWithFree.push_back(blockIndex);
      block.deallocate(unitPtr);
    }
  };

  template<unsigned ElemSize>
  typename fsb_allocator_elem_allocator<ElemSize>::BlocksVector
  fsb_allocator_elem_allocator<ElemSize>::blocksVector;

  template<unsigned ElemSize>
  std::vector<typename fsb_allocator_elem_allocator<ElemSize>::data_t>
  fsb_allocator_elem_allocator<ElemSize>::blocksWithFree;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
  template<unsigned ElemSize>
  fsb_allocator_mutex fsb_allocator_elem_allocator<ElemSize>::mutex;
#endif


  template<unsigned ElemSize>
  class fsb_allocator2_elem_allocator
  {
    static const std::size_t BlockElements = 1024;

    static const std::size_t DSize = sizeof(std::size_t);
    static const std::size_t ElemSizeInDSize = (ElemSize + (DSize-1)) / DSize;
    static const std::size_t BlockSize = BlockElements*ElemSizeInDSize;

    struct Blocks
    {
      std::vector<std::size_t*> ptrs;

      Blocks()
      {
	ptrs.reserve(256);
	ptrs.push_back(new std::size_t[BlockSize]);
      }

      ~Blocks()
      {
	for(std::size_t i = 0; i < ptrs.size(); ++i)
	  delete[] ptrs[i];
      }
    };

    static Blocks blocks;
    static std::size_t headIndex;
    static std::size_t* freeList;
    static std::size_t allocatedElementsAmount;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
    static fsb_allocator_mutex mutex;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_BOOST
    struct Lock: boost::mutex::scoped_lock
    {
      Lock(): boost::mutex::scoped_lock(mutex) {}
    };
#else
    struct Lock
    {
      Lock() { mutex.lock(); }
      ~Lock() { mutex.unlock(); }
    };
#endif
#endif

    static void freeAll()
    {
      for(std::size_t i = 1; i < blocks.ptrs.size(); ++i)
	delete[] blocks.ptrs[i];
      blocks.ptrs.resize(1);
      headIndex = 0;
      freeList = 0;
    }

  public:
    static void* allocate()
    {
#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
      Lock lock;
#endif

      ++allocatedElementsAmount;

      if(freeList)
        {
	  std::size_t* retval = freeList;
	  freeList = reinterpret_cast<std::size_t*>(*freeList);
	  return retval;
        }

      if(headIndex == BlockSize)
        {
	  blocks.ptrs.push_back(new std::size_t[BlockSize]);
	  headIndex = 0;
        }

      std::size_t* retval = &(blocks.ptrs.back()[headIndex]);
      headIndex += ElemSizeInDSize;
      return retval;
    }

    static void deallocate(void* ptr)
    {
      if(ptr)
        {
#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
	  Lock lock;
#endif

	  std::size_t* sPtr = (std::size_t*)ptr;
	  *sPtr = reinterpret_cast<std::size_t>(freeList);
	  freeList = sPtr;

	  if(--allocatedElementsAmount == 0)
	    freeAll();
        }
    }

    static void cleanSweep(std::size_t unusedValue = std::size_t(-1))
    {
#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
      Lock lock;
#endif

      while(freeList)
        {
	  std::size_t* current = freeList;
	  freeList = reinterpret_cast<std::size_t*>(*freeList);
	  *current = unusedValue;
        }

      for(std::size_t i = headIndex; i < BlockSize; i += ElemSizeInDSize)
	blocks.ptrs.back()[i] = unusedValue;

      for(std::size_t blockInd = 1; blockInd < blocks.ptrs.size();)
        {
	  std::size_t* block = blocks.ptrs[blockInd];
	  std::size_t freeAmount = 0;
	  for(std::size_t i = 0; i < BlockSize; i += ElemSizeInDSize)
	    if(block[i] == unusedValue)
	      ++freeAmount;

	  if(freeAmount == BlockElements)
            {
	      delete[] block;
	      blocks.ptrs[blockInd] = blocks.ptrs.back();
	      blocks.ptrs.pop_back();
            }
	  else ++blockInd;
        }

      const std::size_t* lastBlock = blocks.ptrs.back();
      for(headIndex = BlockSize; headIndex > 0; headIndex -= ElemSizeInDSize)
	if(lastBlock[headIndex-ElemSizeInDSize] != unusedValue)
	  break;

      const std::size_t lastBlockIndex = blocks.ptrs.size() - 1;
      for(std::size_t blockInd = 0; blockInd <= lastBlockIndex; ++blockInd)
        {
	  std::size_t* block = blocks.ptrs[blockInd];
	  for(std::size_t i = 0; i < BlockSize; i += ElemSizeInDSize)
            {
	      if(blockInd == lastBlockIndex && i == headIndex)
		break;

	      if(block[i] == unusedValue)
		deallocate(block + i);
            }
        }
    }
  };

  template<unsigned ElemSize>
  typename fsb_allocator2_elem_allocator<ElemSize>::Blocks
  fsb_allocator2_elem_allocator<ElemSize>::blocks;

  template<unsigned ElemSize>
  std::size_t fsb_allocator2_elem_allocator<ElemSize>::headIndex = 0;

  template<unsigned ElemSize>
  std::size_t* fsb_allocator2_elem_allocator<ElemSize>::freeList = 0;

  template<unsigned ElemSize>
  std::size_t fsb_allocator2_elem_allocator<ElemSize>::allocatedElementsAmount = 0;

#ifdef FSBALLOCATOR_USE_THREAD_SAFE_LOCKING_OBJECT
  template<unsigned ElemSize>
  fsb_allocator_mutex fsb_allocator2_elem_allocator<ElemSize>::mutex;
#endif


  template<typename Ty>
  class fsb_allocator
  {
  public:
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef Ty *pointer;
    typedef const Ty *const_pointer;
    typedef Ty& reference;
    typedef const Ty& const_reference;
    typedef Ty value_type;

    pointer address(reference val) const { return &val; }
    const_pointer address(const_reference val) const { return &val; }

    template<class Other>
    struct rebind
    {
      typedef fsb_allocator<Other> other;
    };

    fsb_allocator() throw() {}

    template<class Other>
    fsb_allocator(const fsb_allocator<Other>&) throw() {}

    template<class Other>
    fsb_allocator& operator=(const fsb_allocator<Other>&) { return *this; }

    pointer allocate(size_type count, const void* = 0)
    {
      assert(count == 1);
      return static_cast<pointer>
	(fsb_allocator_elem_allocator<sizeof(Ty)>::allocate());
    }

    void deallocate(pointer ptr, size_type)
    {
      fsb_allocator_elem_allocator<sizeof(Ty)>::deallocate(ptr);
    }

    void construct(pointer ptr, const Ty& val)
    {
      new ((void *)ptr) Ty(val);
    }

    void destroy(pointer ptr)
    {
      ptr->Ty::~Ty();
    }

    size_type max_size() const throw() { return 1; }
  };


  template<typename Ty>
  class fsb_allocator2
  {
  public:
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef Ty *pointer;
    typedef const Ty *const_pointer;
    typedef Ty& reference;
    typedef const Ty& const_reference;
    typedef Ty value_type;

    pointer address(reference val) const { return &val; }
    const_pointer address(const_reference val) const { return &val; }

    template<class Other>
    struct rebind
    {
      typedef fsb_allocator2<Other> other;
    };

    fsb_allocator2() throw() {}

    template<class Other>
    fsb_allocator2(const fsb_allocator2<Other>&) throw() {}

    template<class Other>
    fsb_allocator2& operator=(const fsb_allocator2<Other>&) { return *this; }

    pointer allocate(size_type count, const void* = 0)
    {
      assert(count == 1);
      return static_cast<pointer>
	(fsb_allocator2_elem_allocator<sizeof(Ty)>::allocate());
    }

    void deallocate(pointer ptr, size_type)
    {
      fsb_allocator2_elem_allocator<sizeof(Ty)>::deallocate(ptr);
    }

    void construct(pointer ptr, const Ty& val)
    {
      new ((void *)ptr) Ty(val);
    }

    void destroy(pointer ptr)
    {
      ptr->Ty::~Ty();
    }

    size_type max_size() const throw() { return 1; }

    void cleanSweep(std::size_t unusedValue = std::size_t(-1))
    {
      fsb_allocator2_elem_allocator<sizeof(Ty)>::cleanSweep(unusedValue);
    }
  };

  typedef fsb_allocator2<std::size_t> FSBRefCountAllocator;

} // namespace sl

#endif
