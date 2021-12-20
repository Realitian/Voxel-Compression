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
#ifndef SL_SORT_HPP
#define SL_SORT_HPP

#include <algorithm>
#include <functional>
#include <vector>
#include <cstdlib>

namespace sl {

  /** 
   *  Copy range to incore std::vector, sort, and copy back.
   *  We assume that IT points to a structure slow at random
   *  acces but fast with streaming (as for external memory
   *  arrays)
   */
  template<typename IT, typename BP> 
  void incore_sort(IT begin, IT end, const BP& is_less) {
    typedef typename std::iterator_traits<IT>::value_type value_t;

    std::vector<value_t>* incore_copy = new std::vector<value_t>();
	incore_copy->reserve((std::size_t)(end-begin));
	for (IT it=begin; it!= end; ++it) {
      incore_copy->push_back(*it);
	}
	std::sort(incore_copy->begin(), incore_copy->end(), is_less);
    std::copy(incore_copy->begin(), incore_copy->end(), begin);
    delete incore_copy; incore_copy = 0;
  }

  /** 
   *  Return the iterator containing the median
   *  value, selected among first, last, and mid
   */
  template<typename IT, typename BP>
  IT pivot_median(IT begin, IT end, const BP& is_less) {
    typedef typename std::iterator_traits<IT>::value_type value_t;
    
    IT pivot(begin+(end-begin)/2);
    IT last(end); --last;
    value_t begin_value = *begin;
    value_t pivot_value = *pivot;
    value_t last_value  = *last;
    if (is_less(begin_value, pivot_value) && 
	is_less(last_value, begin_value)) {
      // L<B<P
      pivot=begin;
    } else if (is_less(begin_value, last_value) && 
	       is_less(pivot_value, begin_value)) {
      // P<B<L
      pivot=begin;
    } else if (is_less(last_value, pivot_value) && 
	       is_less(begin_value, last_value)) {
      // B<L<P
      pivot=last;
    } else if (is_less(last_value, begin_value) && 
	       is_less(pivot_value, last_value)) {
      // P<L<B
      pivot=last;
    }
    return pivot;
  }

  template <typename IT>
  void iter_swap(IT& it0, IT& it1) {
	typedef typename std::iterator_traits<IT>::value_type value_t;
  	value_t tmp = *it0;
	*it0 = *it1;
	*it1 = tmp;
  }

  template <typename IT, typename BP> 
  IT quicksort_partition(IT begin, IT end, IT pivot, const BP& is_less) {
    typedef typename std::iterator_traits<IT>::value_type value_t;
    IT first(begin);
    IT last(end); --last;
    IT pivot_store(last);
    value_t pivot_value(*pivot);

    // Move pivot out of the way
    iter_swap(pivot, pivot_store);
    
    for(;;) {
      while (last!=first && is_less(*first, pivot_value)) {
	++first;
      }
      while (last!=first && !is_less(*last, pivot_value)) {
	--last;
      }
      if (last==first) break;
      iter_swap(first, last);
    }
    // Move pivot into final position and return final
    iter_swap(first, pivot_store);
    return first;
  }


  /**
   * Quick-and-dirty implementation of in-place quicksort
   * Quicksort is so inherently sequential that it's practical to run it 
   * on external memory arrays.
   * The implementation is based on ideas from Sedgewick. 
   * In-place partitioning is used. After partitioning, 
   * the partition with the fewest elements is (recursively) 
   * sorted first, requiring at most O(logn) space. 
   * Then the other partition is sorted using iteration. 
   */
  template<typename IT, typename BP>
  void quicksort(IT begin, IT end, const BP& is_less) {
    typedef typename std::iterator_traits<IT>::value_type value_t;

	const typename IT::difference_type INCORE_TMP_MEMORY= 20*1024*1024;
	const typename IT::difference_type INCORE_ITEM_COUNT=INCORE_TMP_MEMORY/sizeof(value_t);

	typename IT::difference_type size(end-begin);
    while (size>1) {
      // Use incore sort if below threshold, otherwise
      // proceed with quicksort partitioning
      if (size<INCORE_ITEM_COUNT) {
	// Copy to incore vector and sort
	incore_sort(begin, end, is_less);
	begin=end; // Nothing remains to be sorted
      } else {
	// Chose a pivot and partition
	IT pivot= pivot_median(begin, end, is_less);
	pivot=quicksort_partition(begin, end, pivot, is_less);

	// Recursively sort smaller partition and
	// iteratively continue with larger one
	// (this trick, due to Sedgewick, reduces 
	// stack depth)
	typename IT::difference_type size_lo=(pivot-begin);
	typename IT::difference_type size_hi=size-size_lo;
	IT pivot_plus1 = pivot; ++pivot_plus1;
	if (size_lo<size_hi) {
	  quicksort(begin, pivot, is_less);
	  begin = pivot_plus1;
	} else {
	  quicksort(pivot_plus1, end, is_less);
	  end = pivot;
	}
      }
      size = (end-begin);
    }
  }

  template<typename IT>
  void quicksort(IT begin, IT end) {	 
    typedef typename std::iterator_traits<IT>::value_type value_t;
    quicksort(begin, end, std::less<value_t>());
  }

  template <class IT>
  bool is_sorted(IT iter1, IT iter2)
  {
    IT prev = iter1;
    for (++iter1;  iter1 != iter2;  ++iter1, ++prev) {
      if (*iter1 < *prev) {
	return false;
      }
    }
    return true;
  }
 
 
  template <class IT, class BP>
  bool is_sorted(IT iter1, IT iter2, BP pred)
  {
    IT prev = iter1;
    for (++iter1;  iter1 != iter2;  ++iter1, ++prev) {
      if (pred(*iter1, *prev)) {
	return false;
      }
    }
    return true;
  }
}

#endif
