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
#ifndef SL_LRU_CACHE_HPP
#define SL_LRU_CACHE_HPP

#include <iostream>
#include <list>
#include <map>
#include <cassert>

namespace sl {
      
  /**
   * Cache maintened using Least Recently Used eviction policy.
   */
  template <class KEY,class DATA>
  class lru_cache {
	
  public:
    typedef KEY				       	      key_t;
    typedef DATA				      data_t;
    typedef lru_cache<KEY,DATA>                       this_t;
    typedef typename std::pair<key_t, data_t>         id_data_pair_t;
    typedef typename std::list<id_data_pair_t>        list_t;
    typedef typename list_t::iterator		      list_iterator_t;
    typedef typename list_t::const_iterator	      const_list_iterator_t;
    typedef typename std::map<key_t,list_iterator_t > map_t;
    typedef typename map_t::iterator		      map_iterator_t;
    typedef typename map_t::const_iterator	      const_map_iterator_t;
	
    typedef list_iterator_t iterator_t;
    typedef const_list_iterator_t const_iterator_t;

  protected:
    std::size_t capacity_;

    list_t      lru_sorted_key_data_list_;
    map_t	id_iterator_map_;
    std::size_t	list_size_;
	
  public:

    /// Create a cache with the given capacity
    lru_cache(std::size_t max_size);

    /// Destruct cache
    ~lru_cache();

    /// Erase all elements
    void clear();

    /// Resize cache to new capacity
    void resize(std::size_t new_capacity);
		
    /// Current number of elements
    std::size_t size() const;

    /// Max number of elements
    std::size_t capacity() const;

    /// Iterator to end of cache.
    iterator_t end();

    /// Iterator to end of cache.
    const_iterator_t end() const;

    /// Iterator to most recently used element
    iterator_t begin();

    /// Iterator to most recently used element
    const_iterator_t begin() const;

    // The oldest item
    const id_data_pair_t& back() const;
    
    // The youngest item
    const id_data_pair_t& front() const;

    /// Find element associated to key. end() if not present.
    iterator_t find(const key_t& key);

    /// Find element associated to key. end() if not present.
    const_iterator_t find(const key_t& key) const;

    /// True iff key is present
    bool has(const key_t& key) const;

    /// Insert new element. Only possible if !has(key)
    void insert(const key_t& key, const data_t& data);

    /// Erase element pointed by iterator
    void erase(const iterator_t& it);

    /// Erase element associated to key if present
    void erase(const key_t& key);
  };     
      
      
} // namespace sl
#endif // SL_LRU_CACHE_HPP

#ifndef SL_LRU_CACHE_IPP
#define SL_LRU_CACHE_IPP

namespace sl {

  template < class KEY,class DATA >
  inline lru_cache<KEY, DATA>::lru_cache(std::size_t max_size) {
    capacity_ = max_size;
    list_size_ = 0;
  }
    
  template < class KEY,class DATA >
  inline lru_cache<KEY, DATA>::~lru_cache() {
    this->clear();
  }

    
  template < class KEY,class DATA >
  inline void lru_cache<KEY, DATA>::clear() {
    assert(list_size_ == lru_sorted_key_data_list_.size());

    lru_sorted_key_data_list_.clear();
    list_size_ = 0;
    id_iterator_map_.clear();
  }

  /**
   * Erase data so that current size is <= to new size.
   */
  template < class KEY,class DATA >
  inline void lru_cache<KEY, DATA>::resize(std::size_t new_size) {
    while (this->size()> new_size) {
      iterator_t it = --(this->end());
      this->erase(it); 
    }
    this->capacity_ = new_size;
  }
    
  template < class KEY,class DATA >
  inline std::size_t lru_cache<KEY, DATA>::size() const { 
    assert(list_size_ == lru_sorted_key_data_list_.size());

    //  return this->lru_sorted_key_data_list_.size();
    return list_size_;
  }
    
  template < class KEY,class DATA >
  inline std::size_t lru_cache<KEY, DATA>::capacity() const {
    return this->capacity_;
  }


  template <class KEY,class DATA >
  typename lru_cache<KEY, DATA>::iterator_t lru_cache<KEY, DATA>::end() {
    return this->lru_sorted_key_data_list_.end();
  }

  template <class KEY,class DATA >
  typename lru_cache<KEY, DATA>::const_iterator_t lru_cache<KEY, DATA>::end() const {
    return this->lru_sorted_key_data_list_.end();
  }
	
  template <class KEY,class DATA >
  typename lru_cache<KEY, DATA>::iterator_t lru_cache<KEY, DATA>::begin() {
    return this->lru_sorted_key_data_list_.begin();
  }
	
  template <class KEY,class DATA >
  typename lru_cache<KEY, DATA>::const_iterator_t lru_cache<KEY, DATA>::begin() const {
    return this->lru_sorted_key_data_list_.begin();
  }
	
  template <class KEY,class DATA >
  const typename lru_cache<KEY, DATA>::id_data_pair_t& lru_cache<KEY, DATA>::back() const {
    return *(--(this->lru_sorted_key_data_list_.end()));
  }

  template <class KEY,class DATA >
  const typename lru_cache<KEY, DATA>::id_data_pair_t& lru_cache<KEY, DATA>::front() const {
    return *(this->lru_sorted_key_data_list_.begin());
  }
   
  template <class KEY,class DATA >
  typename lru_cache<KEY, DATA>::iterator_t lru_cache<KEY, DATA>::find(const key_t& key) {
    assert(list_size_ == lru_sorted_key_data_list_.size());

    map_iterator_t map_it = this->id_iterator_map_.find(key); 
    if (map_it == this->id_iterator_map_.end()) {
      return this->end();
    } else {
      list_iterator_t list_it = map_it->second;
	  	  
      // bring to top
      this->lru_sorted_key_data_list_.push_front((*list_it));
      // update map 
      map_it->second = this->lru_sorted_key_data_list_.begin();
      this->lru_sorted_key_data_list_.erase(list_it);

      return this->begin();
    }

    assert(list_size_ == lru_sorted_key_data_list_.size());
  }
	
  template <class KEY,class DATA >
  typename lru_cache<KEY, DATA>::const_iterator_t lru_cache<KEY, DATA>::find(const key_t& key) const {
    assert(list_size_ == lru_sorted_key_data_list_.size());

    map_iterator_t map_it = const_cast<this_t*>(this)->id_iterator_map_.find(key); 
    if (map_it == const_cast<this_t*>(this)->id_iterator_map_.end()) {
      return this->end();
    } else {
      list_iterator_t list_it = map_it->second;
	  	  
      // bring to top
      const_cast<this_t*>(this)->lru_sorted_key_data_list_.push_front((*list_it));
      // update map 
      const_cast<this_t*>(this)->it->second = lru_sorted_key_data_list_.begin();
      const_cast<this_t*>(this)->lru_sorted_key_data_list_.erase(list_it);

      return this->begin();
    }
    assert(list_size_ == lru_sorted_key_data_list_.size());
  }
	
  template <class KEY,class DATA >
  bool lru_cache<KEY, DATA>::has(const key_t& key) const {
    return this->id_iterator_map_.find(key) != this->id_iterator_map_.end();
  }

  template <class KEY,class DATA >
  void lru_cache<KEY, DATA>::insert(const key_t& key, const data_t& data) {
    assert(!has(key));

    // Erase old elements
    if (this->size()+1 > this->capacity()) {
      iterator_t it = --(this->end());
      this->erase(it); 
    }       

    // Insert at first position
    if (this->capacity() > 0) {
      this->lru_sorted_key_data_list_.push_front(std::make_pair(key,data));     
      this->id_iterator_map_[key]= this->lru_sorted_key_data_list_.begin();
      ++list_size_;
    }
  }

  template <class KEY,class DATA >
  inline void lru_cache<KEY, DATA>::erase(const key_t& key) {
    iterator_t it = this->find(key);
    if (it != this->end()) {
      this->erase(it);
    }
  }

  template <class KEY,class DATA >
  inline void lru_cache<KEY, DATA>::erase(const iterator_t& it) {
    assert(it != this->end());
    assert(list_size_ == lru_sorted_key_data_list_.size());
	
    map_iterator_t map_it = this->id_iterator_map_.find(it->first);
    assert(map_it != this->id_iterator_map_.end());

    this->id_iterator_map_.erase(map_it);
    this->lru_sorted_key_data_list_.erase(it);
    --list_size_;

    assert(list_size_ == lru_sorted_key_data_list_.size());
  }

} // namespace sl

#endif // SL_LRU_CACHE_IPP
