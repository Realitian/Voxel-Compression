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
#ifndef SL_KEYED_HEAP_HPP
#define SL_KEYED_HEAP_HPP

#include <sl/config.hpp>
#include <set>
#include <map>
#include <cassert>

namespace sl {

  /**
   *  A priority queue with elements associated to keys. In addition to
   *  providing insertion of elements, and inspection and removal of the top element,
   *  it provides queries and removal by key. Keys must be hashable, while
   *  (data, key) pairs must be ordered. 
   */
  template <
    class T_KEY, class T_DATA,
    class T_DATA_KEY_COMPARE = std::less< std::pair<T_DATA, T_KEY> >,
    class T_KEY_COMPARE = std::less< T_KEY >
  >
  class keyed_heap {
  public:
    typedef T_KEY               key_t;
    typedef T_DATA              data_t;
    typedef T_DATA_KEY_COMPARE  data_key_compare_t;
    typedef T_KEY_COMPARE       key_compare_t;

    typedef std::pair<data_t, key_t> data_key_pair_t;

    typedef std::set<data_key_pair_t, data_key_compare_t>                     data_key_set_t;
    typedef std::map<key_t, typename data_key_set_t::iterator, key_compare_t> key_iter_map_t;

    typedef typename data_key_set_t::const_reverse_iterator                   const_data_key_iterator_t;
    
    typedef std::size_t         size_type;

  protected:

    // Declaring the following as pointers provides a
    // fast std::swap without the need for redefining it...

    data_key_set_t  *sorted_data_key_pairs_;
    key_iter_map_t  *key_to_sorted_data_key_pairs_iter_;

  public:

    /// Initialize this to empty
    inline keyed_heap() {
      sorted_data_key_pairs_              = new data_key_set_t();
      key_to_sorted_data_key_pairs_iter_  = new key_iter_map_t();
    }

    /// Destruct this
    inline ~keyed_heap() {
      delete sorted_data_key_pairs_; sorted_data_key_pairs_ = 0;
      delete key_to_sorted_data_key_pairs_iter_; key_to_sorted_data_key_pairs_iter_ = 0;
    }

    /// Set the contents of this equal to the contents of rhs
    const keyed_heap& operator=(const keyed_heap& rhs) {
      sorted_data_key_pairs_->clear();
      key_to_sorted_data_key_pairs_iter_->clear();
      typename data_key_set_t::const_iterator ci;
      for (ci = rhs.sorted_data_key_pairs_->begin(); ci != rhs.sorted_data_key_pairs_->end(); ++ci) {
        push(*ci);
      }
      return *this;
    }

    /// Construct this from a copy
    inline keyed_heap(const keyed_heap& rhs) {
      sorted_data_key_pairs_              = new data_key_set_t();
      key_to_sorted_data_key_pairs_iter_  = new key_iter_map_t();
      *this = rhs;
    }

    inline void clear() {
      sorted_data_key_pairs_->clear();
      key_to_sorted_data_key_pairs_iter_->clear();
    }

    /// The number of elements in this
    inline size_type size() const {
      return sorted_data_key_pairs_->size();
    }

    /// True iff this is empty
    inline bool empty() const {
      return sorted_data_key_pairs_->empty();
    }

    /// The element with highest priority
    inline const data_key_pair_t& top() const {
      assert(!empty());
      return *(sorted_data_key_pairs_->rbegin());
    }

    /// True iff key x is present
    inline bool has(const key_t& x) const {
      return key_to_sorted_data_key_pairs_iter_->find(x) != key_to_sorted_data_key_pairs_iter_->end();
    }

    /// Insert x into priority queue, eventually removing oldest value
    inline void push(const data_key_pair_t& x) {
      erase(x.second);
      std::pair<typename data_key_set_t::iterator, bool> it_ok = sorted_data_key_pairs_->insert(x);
      key_to_sorted_data_key_pairs_iter_->insert(std::make_pair(x.second, it_ok.first));
    }
    
    /// Insert x into priority queue. Requires a new x not currently present in the data structure
    inline void push_new(const data_key_pair_t& x) {
      assert(!has(x.second));
      std::pair<typename data_key_set_t::iterator, bool> it_ok = sorted_data_key_pairs_->insert(x);
      key_to_sorted_data_key_pairs_iter_->insert(std::make_pair(x.second, it_ok.first));
    }

    /// Remove highest priority element from queue
    inline void pop() {
      assert(!empty());
      erase(top().second);
    }

    /// Returns an iterator pointing to the beginning of the queue
    inline const_data_key_iterator_t begin() const {
      return sorted_data_key_pairs_->rbegin();
    }

    /// Returns an iterator pointing to the end of the queue
    inline const_data_key_iterator_t end() const {
      return sorted_data_key_pairs_->rend();
    }

    /// Erase the element associated to key x
    inline void erase(const key_t& x) {
      typename key_iter_map_t::iterator it = key_to_sorted_data_key_pairs_iter_->find(x);
      if (it != key_to_sorted_data_key_pairs_iter_->end()) {
        sorted_data_key_pairs_->erase(it->second);
        key_to_sorted_data_key_pairs_iter_->erase(it);
      }
    }

    
  };
}

#endif
