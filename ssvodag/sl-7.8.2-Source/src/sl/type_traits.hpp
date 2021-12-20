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
#ifndef SL_TYPE_TRAITS_HPP
#define SL_TYPE_TRAITS_HPP

#include <sl/config.hpp>

namespace sl {

  namespace detail {

    template <class T>
    struct ptr_traits{};

    template <class T>
    struct ptr_traits<T*> {
      static const bool is_const = false;
      static const bool is_volatile = false;
      typedef T non_const_type;
      typedef T non_volatile_type;
      typedef T unqualified_type;
      static const char* what() { return ""; }
    };

    template <class T>
    struct ptr_traits<const T*> {
      static const bool is_const = true;
      static const bool is_volatile = false;
      typedef T non_const_type;
      typedef const T non_volatile_type;
      typedef T unqualified_type;
      static const char* what() { return "const"; }
    };

    template <class T>
    struct ptr_traits<volatile T*> {
      static const bool is_const = true;
      static const bool is_volatile = false;
      typedef volatile T non_const_type;
      typedef T non_volatile_type;
      typedef T unqualified_type;
      static const char* what() { return "volatile"; }
    };

    template <class T>
    struct ptr_traits<const volatile T*> {
      static const bool is_const = true;
      static const bool is_volatile = true;
      typedef volatile T non_const_type;
      typedef const T non_volatile_type;
      typedef T unqualified_type;
      static const char* what() { return "const volatile"; }
    };

  } // namespace detail

  /**
   *  is a type T  declared const? - is_const<T>
   */
  template <class T>
  struct is_const {
    static const bool value = detail::ptr_traits<T*>::is_const;
  };

  /**
   *  is a type T  declared volatile - is_volatile<T>
   */
  template <class T>
  struct is_volatile {
    static const bool value = detail::ptr_traits<T*>::is_const;
  };

  /** 
   *  is a type T the same as type U? - is_same<T,U>
   */
  template <class T, class U> struct is_same { static const bool value = false; };
  template <class T> struct is_same<T, T> { static const bool value = true; };

  /**
   * is a type T void? - is_void<T>
   */
  template <class T> struct is_void { static const bool value = false; };
  template <> struct is_void<void> { static const bool value = true; };

  /** 
   *  is a pointer T* declared restrict? -- is_restrict<T*>
   */
  template <class T> struct is_restrict { static const bool value = false; };
#if HAVE_RESTRICT_ENABLED
  template <class T> struct is_restrict<T*restrict> { static const bool value = true; };
#endif  

  namespace detail{
    struct any_conversion {
      template <class T>
      any_conversion(const T&);
      template <class T>
      any_conversion(T&);
    };

    template <class T>
    struct checker {
      static double  ok_if_char(any_conversion);
      static char    ok_if_char(T);
    };
  } // namespace detail

  /**
   *  Is the type T_From convertible to the type T_To?
   */
  template <class T_From, class T_To>
  struct is_convertible {
  private:
    static T_From from_;
  public:
    static double ok_if_char(...);
    static char ok_if_char(T_To);
#if defined(__GNUC__)
    static const bool value = ( sizeof(detail::checker<T_To>::ok_if_char(from_)) == sizeof(char) );
#else 
    static const bool value = ( sizeof(ok_if_char(from_)) == sizeof(char) );
#endif
  };

  template <class T_From>
  struct is_convertible<T_From, void> { static const bool value = false; };

  template <class T_To>
  struct is_convertible<void, T_To> { static const bool value = false; };

  template <>
  struct is_convertible<void, void> { static const bool value = true; };

};

#endif /* SL_TYPE_TRAITS_HPP */
