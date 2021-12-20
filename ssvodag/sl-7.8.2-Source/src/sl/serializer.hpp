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
#ifndef SL_SERIALIZER_HPP
#define SL_SERIALIZER_HPP

#include <sl/cstdint.hpp>
#include <map>
#include <algorithm>

namespace sl {

  class input_serializer;
  class output_serializer;

  /**
   *  Serialization traits for objects of type T. By default,
   *  we expect that store_to and retrieve_from are defined
   *  on objects of type T. Redefinition must be made in
   *  case this is not the case.
   */
  template <class T>
  struct serialization_traits {
    static const bool is_array_base_type = false;
    
    static output_serializer& store(output_serializer& s, const T& x) {
      x.store_to(s); return s;
    }
    static input_serializer& retrieve(input_serializer& s, T& x) {
      x.retrieve_from(s); return s;
    }
  };
  
  /// Objects that can be serialized
  class serializable {
  public:
    inline serializable() {}
    inline virtual ~serializable() {};

    virtual void store_to(output_serializer& s) const = 0;
    virtual void retrieve_from(input_serializer& s) = 0;
  };
  
  /**
     The abstract base class for serializing an object.
     It is similar to an ostream.
     The same << notation as for ostreams can be used. It is implemented
     for the basic data types and for classes derived from Osiris. For the
     latter the store function of the object is called. 
  */
  class output_serializer {
  protected: 
    uint32_t version_;
    std::map<void*,uint32_t> pointer_id_;
  public:
    output_serializer(uint32_t v=0);
    virtual ~output_serializer() {
    }
    // returns the version of the dump data format.
    uint32_t version() const { return version_;} 

#if (HAVE_LONG_LONG_BITS >= 64)
    /**
       write an array of signed 64-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const int64_t* p);
#endif

    /**
       write an array of signed 32-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const int32_t* p);

    /**
       write an array of signed 16-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const int16_t* p);

    /**
       write an array of signed 8-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const int8_t* p);

#if (HAVE_LONG_LONG_BITS >= 64)
    /**
       write an array of unsigned 64-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const uint64_t* p);
#endif

    /**
       write an array of unsigned 32-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const uint32_t* p);

    /**
       write an array of unsigned 16-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array */
    virtual void write_array(std::size_t n, const uint16_t* p);
    
    /**
       write an array of unsigned 8-bit integers.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const uint8_t* p);
    
    /**
       write an array of float.
       The default is unsigned the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array */
    virtual void write_array(std::size_t n, const float* p);

    /**
       write an array of double.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const double* p);
    
    /**
       write an array of long double.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const long double* p);

    /**
       write an array of booleans.
       The default is calling the write_simple for each element in the array.
       This function is used by the implementation of dumping the array classes 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param p a pointer to the array
    */
    virtual void write_array(std::size_t n, const bool* p);

    /**
       write a string.
       Perform character set translations when possible.
       The default is calling the unsigned 8-bit integer version of write_array.
       This function is used by the implementation of dumping the std::string class 
       and usually does not need to be called directly. 
       @param n the number of elements
       @param s a pointer to the C-style string
    */

    virtual void write_string(std::size_t n, const char* s);

    /**
       register an object to prepare serializing a pointer to it.
       after writing an object it has to be registered with the dump to
       allow writing a pointer to the object.
    */

    void register_object_address(void* p);

    /**
       serialize a pointer.
       after registering an object's address a pointer to it can be serialized.
       This is done by writing an integer number associated with the object
       when its address is registered.
    */
    void write_pointer(void* p);
  
  public:

#if (HAVE_LONG_LONG_BITS >= 64)
    // implements writing a signed 64-bit integer. The default is using the 32-bit version.
    virtual void write_simple(const int64_t x) ;
#endif
    // implements writing a signed 32-bit integer. This has to be implemented by all classes derived from output_serializer.
    virtual void write_simple(const int32_t x) = 0;
    // implements writing a signed 16-bit integer. The default is using the 32-bit version.
    virtual void write_simple(const int16_t x);
    // implements writing a signed 8-bit integer. The default is using the 16-bit version.
    virtual void write_simple(const int8_t x);
#if (HAVE_LONG_LONG_BITS >= 64)
    // implements writing an unsigned 64-bit integer. The default is using the signed version.
    virtual void write_simple(const uint64_t x);
#endif
    // implements writing an unsigned 32-bit integer. The default is using the signed version.
    virtual void write_simple(const uint32_t x);
    // implements writing an unsigned 16-bit integer. The default is using the signed version.
    virtual void write_simple(const uint16_t x);
    // implements writing an unsigned 8-bit integer. The default is using the signed version.
    virtual void write_simple(const uint8_t x);
    // implements writing a float. The default is using the double version.
    virtual void write_simple(const float x);
    // implements writing a double. This has to be implemented by all classes derived from output_serializer.
    virtual void write_simple(const double x) = 0;
    // implements writing a long double. The default is using the double version.
    virtual void write_simple(const long double x) ;
    // implements writing a bool. The default is using the signed 32-bit integer version.
    virtual void write_simple(const bool x);
  };

} // end namespace sl

/// Serialize an object x to serializer s
template <class T>
sl::output_serializer& operator<<(sl::output_serializer& s, const T& x) {
  sl::serialization_traits<T>::store(s, x);
  return s;
}

namespace sl {
   
  /**
     The abstract base class for de-serializing an object.
     It is similar to an istream.
     The same >> notation as for istreams can be used for basic data types
     Deserialization of Sl objects from a dump is usually implemented
     through a constructor which takes an input_serializer as argument.
     Alternatively the >> operator calls the retrieve function by default.
  */
  class input_serializer {
  protected:
    uint32_t version_;
    std::vector<void*> pointers_;
  public:
    /** The constructor. It takes a version number as
        argument. Different version numbers allow to distinguish
        different data formats, athough this is not used yet.*/
    input_serializer(uint32_t v=0);
    virtual ~input_serializer() {}

    // returns the version of the data format.
    uint32_t version() const { return version_;} 

#if (HAVE_LONG_LONG_BITS >= 64)
    /** read an array of signed 64-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, int64_t* p);
#endif

    /** read an array of signed 32-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, int32_t* p);
    /** read an array of signed 16-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, int16_t* p);
    /** read an array of signed 8-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, int8_t* p);

#if (HAVE_LONG_LONG_BITS >= 64)
    /** read an array of unsigned 64-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, uint64_t* p);
#endif

    /** read an array of unsigned 32-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, uint32_t* p);
    /** read an array of unsigned 16-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, uint16_t* p);
    /** read an array of unsigned 8-bit integers.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, uint8_t* p);
    /** read an array of floats.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, float* p);
    /** read an array of doubles.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, double* p);
    /** read an array of long doubles.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, long double* p);
    /** read an array of booleans.
        The default is calling the read_simple for each element in the array. 
        This function is used by the implementation of deserializing the array classes 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param p a pointer to the array */
    virtual void read_array(std::size_t n, bool* p);

    /** read a C-style string.
        Perform character set translations when possible.
        The default is calling the unsigned 8-bit integer version of read_array.
        This function is used by the implementation of dumping the std::string class 
        and usually does not need to be called directly. 
        @param n the number of elements
        @param s a pointer to memory allocated for a C-style string. At least n+1 bytes must
        be allocated. */
    virtual void read_string(std::size_t n, char* s);
    
    /** register an object to prepare deserializing a pointer to it.
        after reading an object it has to be registered with the dump to
        allow reading a pointer to the object. */
    void register_object_address(void* p);

    /** deserialize a pointer.
        After registering an object's address a pointer to it can be deserialized.
        This is done by reading an integer number associated with the object
        when its address is registered. */
    void* read_pointer();
 
  public:

#if (HAVE_LONG_LONG_BITS >= 64)
    // implements reading a signed 64-bit integer. This has to be implemented by all classes derived from input_serializer.
    virtual void read_simple(int64_t& x) = 0;
#endif
    // implements reading a signed 32-bit integer. This has to be implemented by all classes derived from input_serializer.
    virtual void read_simple(int32_t& x) = 0;
    // implements reading a signed 16-bit integer. The default is using the 32-bit version.
    virtual void read_simple(int16_t& x);
    // implements reading a signed 8-bit integer. The default is using the 16-bit version.
    virtual void read_simple(int8_t& x);
#if (HAVE_LONG_LONG_BITS >= 64)
    // implements reading an unsigned 64-bit integer. The default is using the signed version.
    virtual void read_simple(uint64_t& x);
#endif
    // implements reading an unsigned 32-bit integer. The default is using the signed version.
    virtual void read_simple(uint32_t& x);
    // implements reading an unsigned 16-bit integer. The default is using the signed version.
    virtual void read_simple(uint16_t& x);
    // implements reading an unsigned 8-bit integer. The default is using the signed version.
    virtual void read_simple(uint8_t& x);
    // implements reading a float. The default is using the double version.
    virtual void read_simple(float& x);
    // implements reading a double. This has to be implemented by all classes derived from output_serializer.
    virtual void read_simple(double& x) = 0;
    // implements reading a long double. The default is using the double version.
    virtual void read_simple(long double& x);
    // implements reading a boolean. The default is using the signed 32-bit integer version.
    virtual void read_simple(bool& x);
  };

} // namespace sl

/// Deserialize an object x from serializer s
template <class T>
sl::input_serializer& operator>>(sl::input_serializer& s, T& x) {
  sl::serialization_traits<T>::retrieve(s, x);
  return s;
}

// ======================================
// Helpers
namespace sl {

#define SL_DECLARE_BASE_SERIALIZER(T) \
  template<> struct serialization_traits<T> {\
  static const bool is_array_base_type=true;\
    static output_serializer& store(output_serializer& s, const T& x) {\
      s.write_simple(x); return s;\
    }\
    static input_serializer& retrieve(input_serializer& s, T& x) {\
      s.read_simple(x); return s;\
    }\
  }

#if (HAVE_LONG_LONG_BITS >= 64)
  SL_DECLARE_BASE_SERIALIZER(int64_t);
  SL_DECLARE_BASE_SERIALIZER(uint64_t);
#endif
  
  SL_DECLARE_BASE_SERIALIZER(int32_t);
  SL_DECLARE_BASE_SERIALIZER(uint32_t);
  SL_DECLARE_BASE_SERIALIZER(int16_t);
  SL_DECLARE_BASE_SERIALIZER(uint16_t);
  SL_DECLARE_BASE_SERIALIZER(int8_t);
  SL_DECLARE_BASE_SERIALIZER(uint8_t);
  SL_DECLARE_BASE_SERIALIZER(float);
  SL_DECLARE_BASE_SERIALIZER(double);
  SL_DECLARE_BASE_SERIALIZER(long double);
  SL_DECLARE_BASE_SERIALIZER(bool);

#undef SL_DECLARE_BASE_SERIALIZER

  template <class iterator_t>
  sl::output_serializer& store_container_range(sl::output_serializer& s,
                                                      iterator_t it_bgn,
                                                      iterator_t it_end) {
    for (iterator_t it = it_bgn; it!= it_end; ++it) {
      s << *it;
    }
    return s;
  }
   
  template <class iterator_t>
  sl::input_serializer& retrieve_container_range(sl::input_serializer& s,
                                                        iterator_t it_bgn,
                                                        iterator_t it_end) {
    for (iterator_t it = it_bgn; it!= it_end; ++it) {
      s >> *it;
    }
    return s;
  }
    
  /// serialize a std::-like container
  template<class C>
  inline sl::output_serializer& store_container(sl::output_serializer& s, const C& x) {
    s << uint32_t(x.size());
    return store_container_range(s,x.begin(),x.end());
  }

  /// serialize a fixed size std::-like container
  template<class C>
  inline sl::output_serializer& store_fixed_size_container(sl::output_serializer& s, const C& x) {
    return store_container_range(s,x.begin(),x.end());
  }

  
  /// deserialize a std::list-like container
  template<class C, class value_t>
  sl::input_serializer& retrieve_list_container(sl::input_serializer& s, C& x, const value_t*) {
    uint32_t n; s >> n;
    x=C();
    value_t val;
    while (n--) {
      s >> val;
      x.push_back(val);
    }
    return s;
  }

  /// deserialize a std::list-like container
  template<class C>
  inline sl::input_serializer& retrieve_list_container(sl::input_serializer& s, C& x) {
    return retrieve_list_container(s,x,(0 ? (&(*x.begin())) : 0));
  }
  
  /// deserialize a std::vector-like container
  template<class C>
  inline sl::input_serializer& retrieve_vector_container(sl::input_serializer& s, C& x) {
    uint32_t n; s >> n;
    x.resize(n);
    return retrieve_container_range(s, x.begin(), x.end());
  }

  /// deserialize a fixed-size vector-like container
  template<class C>
  inline sl::input_serializer& retrieve_fixed_size_vector_container(sl::input_serializer& s, C& x) {
    return retrieve_container_range(s, x.begin(), x.end());
  }

  template  <class C, bool OPTIMIZED>
  struct vector_container_serialize_helper {
    static inline sl::input_serializer& retrieve(sl::input_serializer& s, C& x) {
      return retrieve_vector_container(s, x);
    }
    static inline sl::output_serializer& store(sl::output_serializer& s, const C& x) {
      return store_container(s, x);
    }
  };

  template  <class C>
  struct vector_container_serialize_helper<C, true> {
    static inline sl::input_serializer& retrieve(sl::input_serializer& s, C& x) {
      uint32_t n; s >> n;
      x.resize(n);
      if (n) s.read_array(n,&(x[0]));
      return s;
    }
    static inline sl::output_serializer& store(sl::output_serializer& s, const C& x) {
      uint32_t n = uint32_t(x.size());
      s << n;
      if (n) s.write_array(n,&(x[0]));
      return s;
    }
  };
  
  template  <class iterator_t, class const_iterator_t, bool OPTIMIZED>
  struct fixed_size_vector_container_serialize_helper {
    static inline sl::input_serializer& retrieve(sl::input_serializer& s, iterator_t bgn, iterator_t end) {
      return retrieve_container_range(s, bgn, end);
    }
    static inline sl::output_serializer& store(sl::output_serializer& s, const_iterator_t bgn, const_iterator_t end) {
      return store_container_range(s, bgn, end);
    }
  };

  template  <class iterator_t, class const_iterator_t>
  struct fixed_size_vector_container_serialize_helper<iterator_t, const_iterator_t, true> {
    static inline sl::input_serializer& retrieve(sl::input_serializer& s, iterator_t bgn, iterator_t end) {
      uint32_t n = uint32_t(end-bgn);
      if (n) s.read_array(n,&(*(bgn)));
      return s;
    }
    static inline sl::output_serializer& store(sl::output_serializer& s, const_iterator_t bgn, const_iterator_t end) {
      uint32_t n = uint32_t(end-bgn);
      if (n) s.write_array(n,&(*(bgn)));
      return s;
    }
  };

  /// deserialize a std::stack-like container
  template<class C, class value_t>
  sl::input_serializer& retrieve_stack_container(sl::input_serializer& s, C& x, const value_t*) {
    uint32_t n; s >> n;
    x=C(); // empty stack
    value_t elem;
    while (n--) {
      s >> elem;
      x.push(elem);
    }
    return s;
  }

  /// deserialize a std::stack-like container
  template<class C>
  inline sl::input_serializer& retrieve_stack_container(sl::input_serializer& s, C& x) {
    return retrieve_stack_container(s,x,(0 ? (&(*x.begin())) : 0));
  }

  /// deserialize a std::set-like container
  template<class C, class value_t>
  sl::input_serializer& retrieve_set_container(sl::input_serializer& s, C& x, const value_t*) {
    uint32_t n; s >> n;
    x=C(); // empty set
    value_t elem;
    while (n--) {
      s >> elem;
      x.insert(elem);
    }
    return s;
  }

  /// deserialize a std::set-like container
  template<class C>
  inline sl::input_serializer& retrieve_set_container(sl::input_serializer& s, C& x) {
    return retrieve_set_container(s,x,(0 ? (&(*x.begin())) : 0));
  }

} // namespace sl

#endif
