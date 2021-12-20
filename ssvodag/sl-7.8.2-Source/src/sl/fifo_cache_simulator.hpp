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
#ifndef SL_FIFO_CACHE_SIMULATOR_HPP
#define SL_FIFO_CACHE_SIMULATOR_HPP

#include <sl/config.hpp>
#include <deque>

namespace sl {

  template <class T>
  class fifo_cache_simulator {
  public:
    typedef T value_t;
    typedef typename std::deque<value_t>::const_iterator const_iterator;
    
  protected:
    std::deque<value_t> cache_;
    std::size_t   capacity_;
    std::size_t   hit_count_;
    std::size_t   miss_count_;
    
  public:
    
    inline fifo_cache_simulator(std::size_t n = 0) {
      capacity_ = n;
      hit_count_ = 0;
      miss_count_ = 0;
    }

    inline void clear() {
      cache_.clear();
    }

    inline std::size_t capacity() const {
      return capacity_;
    }
    
    inline void resize(std::size_t n) {
      capacity_ = n;
      while (cache_.size()>capacity_) cache_.pop_front();
    }

    inline bool has(const value_t& x) const {
      for (const_iterator it = begin();
           it != end();
           ++it) {
        if ((*it)==x) return true;
      }
      return false;
    }

    inline std::size_t hit_count(const value_t& x) const {
      return has(x) ? std::size_t(1) : std::size_t(0);
    }
    
    inline void insert(const value_t& x) {
      if (has(x)) {
        ++hit_count_;
      } else {
        ++miss_count_;
        cache_.push_back(x);
        if (cache_.size()>capacity_) cache_.pop_front();
      }
    }

    inline void reset_counts() {
      hit_count_ = 0;
      miss_count_ = 0;
    }

    inline std::size_t hit_count() const {
      return hit_count_;
    }

    inline std::size_t miss_count() const {
      return miss_count_;
    }

    inline const_iterator begin() const {
      return cache_.begin();
    }

    inline const_iterator end() const {
      return cache_.end();
    }
    
  };
  
}

#endif





