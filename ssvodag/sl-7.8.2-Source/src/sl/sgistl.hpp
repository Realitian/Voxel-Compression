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
 * Definition of namespace sgistl, and inclusion of SGI-specific
 * STL extensions in a hopefully compiler-independent way.
 */
#ifndef SGISTL_HPP
#define SGISTL_HPP

#include <sl/config.hpp>

#ifdef __GNUC__
#if __GNUC__ < 3
#include <hash_map.h>
#include <hash_set.h>
namespace sgistl { using ::hash_map; using ::hash_set; }; // inherit globals
#else
#include <ext/hash_map>
#include <ext/hash_set>
#if __GNUC_MINOR__ == 0
//namespace sgistl = std;               // GCC 3.0
#define sgistl std 
#else
//namespace sgistl = __gnu_cxx;       // GCC 3.1 and later
#define sgistl  __gnu_cxx
#endif
#endif
#else      // ...  there are other compilers, right?
#include <hash_map>
#include <hash_set>
//namespace sgistl = std;
#endif

//----------- Serialization support stuff, move elsewhere
#include <sl/serializer.hpp>

namespace sl {
  
// FIXME should check sgistl version instead
#ifdef _WIN32 
#define sgistl stdext
#define SL_HASHCOMPARE_DECL class __HashCmp__
#define SL_HASHCOMPARE_DECL2 __HashCmp__
#else
#define SL_HASHCOMPARE_DECL class __HashFn__, class __EqCmp__
#define SL_HASHCOMPARE_DECL2 __HashFn__, __EqCmp__
#endif    

  /// Serialization traits for sgistl::hash_set container
  template <class Key, SL_HASHCOMPARE_DECL, class Allocator>
  struct serialization_traits< sgistl::hash_set<Key,SL_HASHCOMPARE_DECL2,Allocator> > {
    typedef sgistl::hash_set<Key,SL_HASHCOMPARE_DECL2,Allocator> serializable_t;
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

  /// Serialization traits for sgistl::hash_set container
  template <class Key, SL_HASHCOMPARE_DECL, class Allocator>
  struct serialization_traits< sgistl::hash_multiset<Key,SL_HASHCOMPARE_DECL2,Allocator> > {
    typedef sgistl::hash_multiset<Key,SL_HASHCOMPARE_DECL2,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

  /// Serialization traits for sgistl::hash_set container
  template <class Key, class T, SL_HASHCOMPARE_DECL, class Allocator>
  struct serialization_traits< sgistl::hash_map<Key,T,SL_HASHCOMPARE_DECL2,Allocator> > {
    typedef sgistl::hash_map<Key,T,SL_HASHCOMPARE_DECL2,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      sl::uint32_t n; s >> n;
      x=serializable_t();
      Key k;
      T v;
      while (n--) {
        s >> k >> v;
        x[k]=v;
      }
      return s;
    }
  };

  /// Serialization traits for sgistl::hash_set container
  template <class Key, class T, SL_HASHCOMPARE_DECL, class Allocator>
  struct serialization_traits< sgistl::hash_multimap<Key,T,SL_HASHCOMPARE_DECL2,Allocator> > {
    typedef sgistl::hash_multimap<Key,T,SL_HASHCOMPARE_DECL2,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };
}

#endif
