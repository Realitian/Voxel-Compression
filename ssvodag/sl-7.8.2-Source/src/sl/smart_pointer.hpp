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
#ifndef SL_SMART_POINTER_HPP
#define SL_SMART_POINTER_HPP

#include <sl/utility.hpp>      // for generic superclass interface
#include <sl/type_traits.hpp>  // for is_convertible
#include <algorithm>           // for std::swap
#include <functional>          // for std::less

namespace sl {

  /**
   * sized_raw_array_pointer mimics a built-in pointer except that it also
   * contains a size, used to validate access to elements.
   */
  template <class T> 
  class sized_raw_array_pointer {
  protected:
    T*     pointer_;
    size_t count_;
  public:
    typedef T value_t;

    explicit inline sized_raw_array_pointer( T* p=0,  size_t sz=0) 
      : 
      pointer_(p), count_(sz) 
    {
      SL_REQUIRE("Null implies zero size", (p!=0) || (sz==0));
    }  // never throws

    /// The number of elements in the array
    inline size_t count() const {
      return count_;
    }
      
    /// The first element of the array
    inline T& operator*() const { 
      SL_REQUIRE("Not empty", pointer_);
      return *pointer_; 
    }

    /// Pointer to the first element of the array
    inline T* operator->() const { 
      return pointer_; 
    }

    /// Pointer to the first element of the array
    inline T* raw_pointer() const        { 
      return pointer_; 
    } 

    /// The i-th element of the array
    inline T& operator[](std::size_t i) const { 
      SL_REQUIRE("Good index", i<count());
      return pointer_[i]; 
    }

    /// Is the pointer null?
    inline operator bool() const { 
      return pointer_; 
    }
    
  };  // sized_raw_array_pointer

  /**
   *  A const pointer from a non-const one
   */
  template <class T>
  static inline sl::sized_raw_array_pointer<const T> to_constant_pointer(const sl::sized_raw_array_pointer<T>& ptr) {
    return sl::sized_raw_array_pointer<const T>(ptr.raw_pointer(), ptr.size());
  }
  
    
} // namespace sl

namespace sl {

  /**
   * scoped_pointer mimics a built-in pointer except that it guarantees deletion
   * of the object pointed to, either on destruction of the scoped_pointer or via
   * an explicit reset(). 
   */
  template <class T> 
  class scoped_pointer {
  private:
    scoped_pointer(const scoped_pointer&); // non copyable
  protected:
    T* pointer_;
  public:
    typedef T value_t;

    explicit inline scoped_pointer( T* p=0 ) : pointer_(p) {}  // never throws
    inline ~scoped_pointer()                 { delete pointer_; }

    inline void reset( T* p=0 )          { if ( pointer_ != p ) { delete pointer_; pointer_ = p; } }
    inline T& operator*() const          { return *pointer_; }  // never throws
    inline T* operator->() const         { return pointer_; }   // never throws
    inline T* raw_pointer() const        { return pointer_; }   // never throws

    inline operator bool() const { return pointer_; }
  };  // scoped_pointer

} // namespace sl

namespace sl {

  /**
   * scoped_raw_array_pointer extends scoped_pointer to arrays. Deletion of the array pointed to
   * is guaranteed, either on destruction of the scoped_raw_array_pointer or via an explicit
   * reset(). 
   */
  template<class T> 
  class scoped_raw_array_pointer {
  private:
    scoped_raw_array_pointer(const scoped_raw_array_pointer&); // non copyable
  protected:
    T* pointer_;

  public:
    typedef T value_t;
    
    explicit scoped_raw_array_pointer( T* p=0 ) : pointer_(p) {}  // never throws
    ~scoped_raw_array_pointer()                    { delete [] pointer_; }
    
    void reset( T* p=0 )               { if ( pointer_ != p ) {delete [] pointer_; pointer_=p;} }
    
    T* raw_pointer() const                     { return pointer_; }  // never throws
    T& operator[](std::size_t i) const { return pointer_[i]; }  // never throws

    inline operator bool() const { return pointer_; }
  };  // scoped_raw_array_pointer

} // namespace sl

namespace sl {
  
  /**
   *  Base class for implementing reference counted objects
   */
  class reference_counted {
  protected:
    long  reference_count_;
  public:

    /// Initialize a reference counted object
    inline reference_counted() :
      reference_count_(0) {
    }

    /// The number of pointers referencing this object
    inline long use_count() const        { 
      return reference_count_; 
    }

    /// Reference this object
    inline void ref() {
      ++reference_count_;
    }

    /// Dereference this object
    inline void deref() {
      SL_REQUIRE("Was referenced", use_count() > 0);
      --reference_count_;
    }

    inline long* refcount_pointer() const {
      return (long*)(&reference_count_);
    }

  };

  template <class T>
  struct is_reference_counted {
    static const bool value = is_convertible<T*,reference_counted*>::value;
  };

} // namespace sl

namespace sl {
 
  namespace detail {

    /**
     *  Base class for implementing pointers with reference counted copy
     *  semantics. Specialized versions implement intrusive and non
     *  intrusive reference counts
     */
    template<class I, class T, bool B_T_HAS_EMBEDDED_REFCOUNT> 
    class shared_pointer_base {
    public:
      typedef I derived_t;
      typedef T value_t;
    };

    /**
     *  Base class for implementing pointers with reference counted copy
     *  semantics. Non-intrusive version, allocating a reference count
     *  on the heap.
     */
    template<class I, class T> 
    class shared_pointer_base<I,T,false> {
    public:
      static const bool is_intrusive = false;

      typedef shared_pointer_base<I,T,false> this_t;
      typedef I derived_t;
      typedef T value_t;
      
      /// Features to access this as a derived_t pointer
      SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);

    protected:

      T*     pointer_;             // pointer to reference counted object
      long*  refcountpointer_;     // pointer to reference counter

#if 0
    protected:
#else
    public: // MIPSPro bug workaround
#endif

      inline void delete_pointer() {
	derived_ref().delete_pointer();
      }

      inline void deref() {
	SL_REQUIRE("Not void", pointer_);
	SL_CHECK("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	--*refcountpointer_;
	if (*refcountpointer_ == 0) { 
	  delete_pointer(); pointer_ = NULL;
	  delete refcountpointer_; refcountpointer_ = NULL;
	}
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
      }
      
      inline void ref() {
	SL_REQUIRE("Not void", pointer_);
	SL_CHECK("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	++*(refcountpointer_);
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
      }
      
      inline void share(T* rpointer, long* rrefcountpointer) {
	SL_CHECK("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	if (refcountpointer_ != rrefcountpointer) {
	  if (pointer_) deref();
	  pointer_ = rpointer;
	  refcountpointer_ = rrefcountpointer;
	  if (pointer_) ref();
	}
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
      }
  
    public:
      
      explicit inline shared_pointer_base(T* p, long* r) : 
	pointer_(p), 
	refcountpointer_(r ? r : (p ? new long(0) : NULL)) {
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	if (pointer_) ref();
      }

      inline ~shared_pointer_base() {
	SL_CHECK("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	if (pointer_) deref();      
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
      }

      inline void reset(T* p=0) {
	SL_CHECK("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	if ( pointer_ != p) {
	  if (pointer_) {
	    SL_CHECK("Has pointer", pointer_);
	    SL_CHECK("Has refcount", refcountpointer_);
	    --*refcountpointer_;
	    if (*refcountpointer_ == 0) {
	      delete_pointer();
	      if (p) {
		*refcountpointer_ = 1;
		pointer_ = p;
	      } else {
		delete refcountpointer_; refcountpointer_ = NULL;
		pointer_ = NULL;
	      }
	    } else if (p) {
	      refcountpointer_ = new long(1);
	      pointer_ = p;
	    } else {
	      refcountpointer_ = NULL;
	      pointer_ = NULL;
	    }
	  } else {
	    SL_CHECK("Null pointer", !pointer_);
	    SL_CHECK("No refcount", !refcountpointer_);
	    pointer_ = p;
	    if (pointer_) {
	      refcountpointer_ = new long(1);
	    }
	  }
	}
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
      }

      /// The embedded raw pointer
      inline T* raw_pointer() const { return pointer_; }

      /// Is the pointer non void?
      inline operator bool() const { return pointer_; }

      inline long* refcount_pointer() const { 
	return refcountpointer_; 
      }

      /// The number of references to the pointed object
      inline long use_count() const { 
	return refcount_pointer() ? *refcount_pointer() : 0; 
      }

      /// Is the pointed object referenced only by one pointer()?
      inline bool is_unique() const { return use_count() == 1; }

      inline void swap(this_t& other) {  // never throws
	SL_CHECK("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
	std::swap(pointer_,other.pointer_); 
	std::swap(refcountpointer_,other.refcountpointer_); 
	SL_ENSURE("Pointer and refcount consistent", (!pointer_) == (!refcountpointer_));
      } 
    };  // shared_pointer_base

    /**
     *  Base class for implementing pointers with reference counted copy
     *  semantics. Intrusive version, specialized for reference_counted
     *  objects with an embedded reference count.
     */
    template<class I, class T> 
    class shared_pointer_base<I,T,true> {
    public:
      static const bool is_intrusive = true;

      typedef shared_pointer_base<I,T,true> this_t;
      typedef I derived_t;
      typedef T value_t;
      
      /// Features to access this as a derived_t pointer
      SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);

    protected:

      T*     pointer_;             // pointer to reference counted object

#if 0
    protected:
#else
    public: // MIPSPro bug workaround
#endif

      inline void delete_pointer() {
	derived_ref().delete_pointer();
      }


      inline void deref() {
	SL_REQUIRE("Not null", pointer_);
	pointer_->deref();
	if (use_count() == 0) {
	  delete_pointer();
	  pointer_ = NULL;
	}
      }
      
      inline void ref() {
	SL_REQUIRE("Not null", pointer_);
	pointer_->ref();
      }
      
      inline void share(T* rpointer, long*) {
	if (pointer_ != rpointer) {
	  if (pointer_) deref();
	  pointer_ = rpointer;
	  if (pointer_) ref();
	}
      }
  
    public:
      
      explicit inline shared_pointer_base(T* p, long*): pointer_(p) {
	if (pointer_) ref();
      }

      inline ~shared_pointer_base() {
	if (pointer_) deref();      
      }

      inline void reset(T* p=0) {
	if ( pointer_ != p ) {
	  if (pointer_) deref();
	  pointer_ = p;
	  if (pointer_) ref();
	}
      }

      /// The embedded raw pointer
      inline T* raw_pointer() const { return pointer_; }

      /// Is the pointer non void?
      inline operator bool() const { return pointer_; }

      /// The number of references to the pointed object
      inline long use_count() const { 
	return refcount_pointer() ? *refcount_pointer() : 0; 
      }

      /// Is the pointed object referenced only by one pointer()?
      inline bool is_unique() const { return use_count() == 1; }

      inline long* refcount_pointer() const { 
	return pointer_ ? pointer_->refcount_pointer() : NULL;
      }

      inline void swap(this_t& other) {  // never throws
	std::swap(pointer_,other.pointer_); 
      } 
    };  // shared_pointer_base

  } // namespace detail

} // namespace sl

namespace sl {

  /**
   * An enhanced relative of scoped_pointer with reference counted copy semantics.
   * The object pointed to is deleted when the last shared_pointer pointing to it
   * is destroyed or reset.
   */
  template<class T> 
  class shared_pointer: public detail::shared_pointer_base<shared_pointer<T>, T, is_reference_counted<T>::value > {
  public:
    typedef detail::shared_pointer_base<shared_pointer<T>, T, is_reference_counted<T>::value> super_t;
    typedef T value_t;
#if 0
  protected:
    friend class super_t;
#else
  public: // MIPSPro bug workaround
#endif

    inline void delete_pointer() {
      delete this->pointer_;
    }
  public:
    explicit inline shared_pointer(T* p =0) : super_t(p, NULL) {
    }
    
    inline shared_pointer(const shared_pointer& r): super_t(r.raw_pointer(), r.refcount_pointer()) { 
    }  // never throws
  
    inline shared_pointer& operator=(const shared_pointer& r) {
      this->share(r.raw_pointer(),r.refcount_pointer());
      return *this;
    }
    
    template<class Y>
    inline shared_pointer(const shared_pointer<Y>& r) : super_t(r.raw_pointer(),r.refcount_pointer()) {  // never throws 
    }

    template<class Y>
    inline shared_pointer& operator=(const shared_pointer<Y>& r) { 
      this->share(r.raw_pointer(),r.refcount_pointer());
      return *this;
    }

    inline T& operator*() const          { return *(this->pointer_); }  // never throws
    inline T* operator->() const         { return (this->pointer_); }  // never throws
  };  // shared_pointer


  /**
   *  A const pointer from a non-const one
   */
  template <class T>
  static inline sl::shared_pointer<const T> to_constant_pointer(const sl::shared_pointer<T>& ptr) {
    return sl::shared_pointer<const T>(ptr);
  }

} // namespace sl
  
template<class T, typename U>
inline bool operator==(const sl::shared_pointer<T>& a, const sl::shared_pointer<U>& b) { 
  return a.raw_pointer() == b.raw_pointer(); 
}

template<class T, typename U>
inline bool operator!=(const sl::shared_pointer<T>& a, const sl::shared_pointer<U>& b) { 
  return a.raw_pointer() != b.raw_pointer(); 
}

namespace sl {

  /**
   * shared_raw_array_pointer extends shared_pointer to arrays.
   * The array pointed to is deleted when the last shared_raw_array_pointer pointing to it
   * is destroyed or reset.
   */
  template<class T> 
  class shared_raw_array_pointer: public detail::shared_pointer_base<shared_pointer<T>, T, false> {
  public:
    typedef detail::shared_pointer_base<shared_pointer<T>, T, is_reference_counted<T>::value> super_t;
    typedef T value_t;

#if 0
  protected:
    friend class super_t;
#else
  public: // MIPSPro bug workaround
#endif

    inline void delete_pointer() {
      delete[] this->pointer_;
    }

  public:
    explicit inline shared_raw_array_pointer(T* p =0) : super_t(p, NULL) {
    }
    
    inline shared_raw_array_pointer(const shared_raw_array_pointer& r): super_t(r.raw_pointer(), r.refcount_pointer()) { 
    }  // never throws
  
    inline shared_raw_array_pointer& operator=(const shared_raw_array_pointer& r) {
      this->share(r.raw_pointer(),r.refcount_pointer());
      return *this;
    }

    inline T& operator[](std::size_t i) const { return (this->pointer_)[i]; }  // never throws

  };  // shared_raw_array_pointer

  /**
   *  A const pointer from a non-const one
   */
  template <class T>
  static inline sl::shared_raw_array_pointer<const T> to_constant_pointer(const sl::shared_raw_array_pointer<T>& ptr) {
    return sl::shared_raw_array_pointer<const T>(ptr);
  }

} // namespace sl
  
template<class T, typename U>
inline bool operator==(const sl::shared_raw_array_pointer<T>& a, const sl::shared_raw_array_pointer<U>& b) { 
  return a.raw_pointer() == b.raw_pointer(); 
}

template<class T, typename U>
inline bool operator!=(const sl::shared_raw_array_pointer<T>& a, const sl::shared_raw_array_pointer<U>& b) { 
  return a.raw_pointer() != b.raw_pointer(); 
}

//  specializations for things in namespace std  -----------------------------//

namespace std {

  // Specialize std::swap to use the fast, non-throwing swap that's provided
  // as a member function instead of using the default algorithm which creates
  // a temporary and uses assignment.

  template<class T>
  inline void swap(sl::shared_pointer<T>& a, sl::shared_pointer<T>& b)
    { a.swap(b); }

  template<class T>
  inline void swap(sl::shared_raw_array_pointer<T>& a, sl::shared_raw_array_pointer<T>& b)
    { a.swap(b); }

  // Specialize std::less so we can use shared pointers and arrays as keys in
  // associative collections.

  // It's still a controversial question whether this is better than supplying
  // a full range of comparison operators (<, >, <=, >=).
  
  template<class T>
  struct less< sl::shared_pointer<T> >
    : binary_function<sl::shared_pointer<T>, sl::shared_pointer<T>, bool>
  {
    bool operator()(const sl::shared_pointer<T>& a,
        const sl::shared_pointer<T>& b) const
      { return less<T*>()(a.raw_pointer(),b.raw_pointer()); }
  };

  template<class T>
  struct less< sl::shared_raw_array_pointer<T> >
    : binary_function<sl::shared_raw_array_pointer<T>, sl::shared_raw_array_pointer<T>, bool>
  {
    bool operator()(const sl::shared_raw_array_pointer<T>& a,
        const sl::shared_raw_array_pointer<T>& b) const
      { return less<T*>()(a.raw_pointer(),b.raw_pointer()); }
  };

} // namespace std

#endif  // SL_SMART_POINTER_HPP


