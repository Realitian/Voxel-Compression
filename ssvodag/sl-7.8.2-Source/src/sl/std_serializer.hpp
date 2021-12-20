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
#ifndef SL_STD_SERIALIZER_HPP
#define SL_STD_SERIALIZER_HPP

#include <sl/serializer.hpp>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <map>
#if HAVE_STLEXT_UNORDERED_CONTAINERS
#include <sl/stlext_unordered_containers.hpp>
#endif

namespace sl {
  
  /// Serialization traits for std::pair container
  template <class T1, class T2>
  struct serialization_traits< std::pair<T1, T2> > {
    typedef std::pair<T1, T2> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return s << x.first << x.second;
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return s >> x.first >> x.second;
    }
  };

  /// Serialization traits for std::list container
  template <class T, class Allocator>
  struct serialization_traits< std::list<T,Allocator> > {
    typedef std::list<T,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_list_container(s,x);
    }
  };

  /// Serialization traits for std::vector container
  template <class T, class Allocator>
  struct serialization_traits< std::vector<T,Allocator> > {
    typedef std::vector<T,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::vector_container_serialize_helper<serializable_t,sl::serialization_traits<T>::is_array_base_type>::store(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::vector_container_serialize_helper<serializable_t,sl::serialization_traits<T>::is_array_base_type>::retrieve(s,x);
    }
  };

  /// Serialization traits for std::string container
  template <class charT, class traits, class Allocator>
  struct serialization_traits< std::basic_string<charT,traits,Allocator> > {
    typedef std::basic_string<charT,traits,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      s << sl::uint32_t(x.size());
      if(x.size()) {
        s.write_string(x.size()+1,x.c_str());
      }
      return s;
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      sl::uint32_t sz; s >> sz;
      if(sz) {
        charT* t = new charT[sz+1];
        s.read_string(sz+1,t);
        if(t[sz]!=charT(0)) SL_FAIL("string on serializer not terminating with '\\0'");
        x=t;
        delete [] t;
        if(x.length()!=sz) SL_FAIL("string on serializer has incorrect length");
      } else {
        x="";
      }
      return s;
    }
  };

  /// Serialization traits for std::set container
  template <class Key, class Compare, class Allocator>
  struct serialization_traits< std::set<Key,Compare,Allocator> > {
    typedef std::set<Key,Compare,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

  /// Serialization traits for std::set container
  template <class Key, class Compare, class Allocator>
  struct serialization_traits< std::multiset<Key,Compare,Allocator> > {
    typedef std::multiset<Key,Compare,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

  /// Serialization traits for std::set container
  template <class Key, class T, class Compare, class Allocator>
  struct serialization_traits< std::map<Key,T,Compare,Allocator> > {
    typedef std::map<Key,T,Compare,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      sl::uint32_t n; s >> n;
      x=std::map<Key,T,Compare,Allocator>();
      Key k;
      T v;
      while (n--) {
        s >> k >> v;
        x[k]=v;
      }
      return s;
    }
  };

  /// Serialization traits for std::set container
  template <class Key, class T, class Compare, class Allocator>
  struct serialization_traits< std::multimap<Key,T,Compare,Allocator> > {
    typedef std::multimap<Key,T,Compare,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

} // namespace sl


#if HAVE_STLEXT_UNORDERED_CONTAINERS

namespace sl {

  /// Serialization traits for stlext::unordered_set container
  template <class Key, class __HashFn__, class __EqCmp__, class Allocator>
  struct serialization_traits< stlext::unordered_set<Key,__HashFn__, __EqCmp__,Allocator> > {
    typedef stlext::unordered_set<Key,__HashFn__, __EqCmp__,Allocator> serializable_t;
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

  /// Serialization traits for stlext::unordered_set container
  template <class Key, class __HashFn__, class __EqCmp__, class Allocator>
  struct serialization_traits< stlext::unordered_multiset<Key,__HashFn__, __EqCmp__,Allocator> > {
    typedef stlext::unordered_multiset<Key,__HashFn__, __EqCmp__,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };

  /// Serialization traits for stlext::unordered_set container
  template <class Key, class T, class __HashFn__, class __EqCmp__, class Allocator>
  struct serialization_traits< stlext::unordered_map<Key,T,__HashFn__, __EqCmp__,Allocator> > {
    typedef stlext::unordered_map<Key,T,__HashFn__, __EqCmp__,Allocator> serializable_t;
    
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

  /// Serialization traits for stlext::unordered_set container
  template <class Key, class T, class __HashFn__, class __EqCmp__, class Allocator>
  struct serialization_traits< stlext::unordered_multimap<Key,T,__HashFn__, __EqCmp__,Allocator> > {
    typedef stlext::unordered_multimap<Key,T,__HashFn__, __EqCmp__,Allocator> serializable_t;
    
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const serializable_t& x) {
      return sl::store_container(s,x);
    }
    static input_serializer& retrieve(input_serializer& s, serializable_t& x) {
      return sl::retrieve_set_container(s,x);
    }
  };
}
#endif // HAVE_STLEXT_UNORDERED_CONTAINERS

#endif

