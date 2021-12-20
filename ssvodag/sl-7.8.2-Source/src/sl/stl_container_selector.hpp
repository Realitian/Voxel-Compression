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

#ifndef SL_STL_CONTAINER_SELECTOR_HPP
#define SL_STL_CONTAINER_SELECTOR_HPP

#include <sl/hash.hpp>
#include <map>
#include <set>
#if HAVE_STLEXT_UNORDERED_CONTAINERS
#include <sl/stlext_unordered_containers.hpp>
#endif

namespace sl {

  /**
   * A selector of set implementations. If stl provides unordered containers,
   * hash-table based containers are select as unordered and fast implementations,
   * otherwise, all versions revert to the standard ordered versions.
   */
  template <class Key,
	    class Hash = sl::hash<Key>,
	    class Equal = std::equal_to<Key>,
	    class Less = std::less<Key>,
	    class Allocator = std::allocator<Key> >
  struct stl_set_selector {
    typedef std::set<Key,Less,Allocator> ordered_t;
    typedef std::set<Key,Less,Allocator> multi_ordered_t;
#if HAVE_STLEXT_UNORDERED_CONTAINERS
    typedef stlext::unordered_set<Key,Hash,Equal,Allocator> unordered_t;
    typedef stlext::unordered_multiset<Key,Hash,Equal,Allocator> multi_unordered_t;
#else
    typedef ordered_t unordered_t;
    typedef multi_ordered_t multi_unordered_t;
#endif
    typedef unordered_t fast_t;
    typedef multi_unordered_t multi_fast_t;
  };
  
  /**
   * A selector of map implementations. If stl provides unordered containers,
   * hash-table based containers are select as unordered and fast implementations,
   * otherwise, all versions revert to the standard ordered versions.
   */
  template <class Key,
	    class Item,
	    class Hash = sl::hash<Key>,
	    class Equal = std::equal_to<Key>,
	    class Less = std::less<Key>,
	    class Allocator = std::allocator<std::pair<const Key, Item> > >
  struct stl_map_selector {
    typedef std::map<Key,Item,Less,Allocator> ordered_t;
    typedef std::map<Key,Item,Less,Allocator> multi_ordered_t;
#if HAVE_STLEXT_UNORDERED_CONTAINERS
    typedef stlext::unordered_map<Key,Item,Hash,Equal,Allocator> unordered_t;
    typedef stlext::unordered_multimap<Key,Item,Hash,Equal,Allocator> multi_unordered_t;
#else
    typedef ordered_t unordered_t;
    typedef multi_ordered_t multi_unordered_t;
#endif
    typedef unordered_t fast_t;
    typedef multi_unordered_t multi_fast_t;
  };

} // namespace sl
  
#endif
