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
#ifndef SL_RESERVABLE_PRIORITY_QUEUE_HPP
#define SL_RESERVABLE_PRIORITY_QUEUE_HPP

#include <sl/serializer.hpp>
#include <queue>
#include <vector>

namespace sl {
  
  /**
   *  A vector-based stl priority queue with the ability of preallocating space
   */
  template <class T, class Compare=std::less<T> >
  class reservable_priority_queue: public std::priority_queue<T, std::vector<T>, Compare> {
  public:
    typedef T value_type; 
    typedef typename std::priority_queue<T>::container_type container_type;
    typedef typename std::priority_queue<T>::size_type size_type;

  public:
    reservable_priority_queue(size_type capacity = 0) { reserve(capacity); };
    void reserve(size_type capacity) { this->c.reserve(capacity); } 
    size_type capacity() const { return this->c.capacity(); } 
  };

  // FIXME: Serialization not yet implemented
}

#endif

