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

/**
 * Definition of namespace stlext, and inclusion of unorderer set/map
 * STL extensions in a hopefully compiler-independent way.
 */
#ifndef STL_UNORDERED_CONTAINERS_HPP
#define STL_UNORDERED_CONTAINERS_HPP

#include <sl/config.hpp>
#include <sl/hash.hpp>

#if HAVE_STL_UNORDERED_CONTAINERS
#include <unordered_set>
#include <unordered_map>
#define STLEXT std

#elif HAVE_STL_TR1_UNORDERED_CONTAINERS
#include <tr1/unordered_set>
#include <tr1/unordered_map>

#define STLEXT std::tr1

#elif HAVE_STL_TR1_UNORDERED_CONTAINERS_MSVC
#include <unordered_set>
#include <unordered_map>

#define STLEXT std::tr1

#else

#error "No unordered containers found!"

#endif


namespace stlext {
  using STLEXT::unordered_map;
  using STLEXT::unordered_set;
  using STLEXT::unordered_multimap;
  using STLEXT::unordered_multiset;
}

#if 0

/// The following is really nice but can break some software...

// Redirect to sl hasher...
template <class T> inline std::size_t STLEXT::hash<T>::operator()(T v) const {
  sl::hash<T> hasher;
  return hasher(v);
}

#endif

#endif
