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
#ifndef SL_INTEGER_HPP
#define SL_INTEGER_HPP

#include <sl/assert.hpp>

namespace sl {

  namespace detail {
    ///  fast integers from least integers
    template< class LeastInt >
    struct int_fast_t { typedef LeastInt fast; }; // imps may specialize
    
    ///  Convert signed integer category to type 
    template< int Category > struct int_least_helper {}; // default is empty
    ///  Convert unsigned integer category to type 
    template< int Category > struct uint_least_helper {}; // default is empty
    
#if HAVE_LONG_LONG
    template<> struct int_least_helper<1> { typedef long long least; };
    template<> struct int_least_helper<2> { typedef long least; };
    template<> struct int_least_helper<3> { typedef int least; };
    template<> struct int_least_helper<4> { typedef short least; };
    template<> struct int_least_helper<5> { typedef signed char least; };

    template<> struct uint_least_helper<1> { typedef unsigned long long least; };
    template<> struct uint_least_helper<2> { typedef unsigned long least; };
    template<> struct uint_least_helper<3> { typedef unsigned int least; };
    template<> struct uint_least_helper<4> { typedef unsigned short least; };
    template<> struct uint_least_helper<5> { typedef unsigned char least; };
#elif HAVE_MSVC_INT64
    template<> struct int_least_helper<1> { typedef __int64 least; };
    template<> struct int_least_helper<2> { typedef long least; };
    template<> struct int_least_helper<3> { typedef int least; };
    template<> struct int_least_helper<4> { typedef short least; };
    template<> struct int_least_helper<5> { typedef signed char least; };

    template<> struct uint_least_helper<1> { typedef unsigned __int64 least; };
    template<> struct uint_least_helper<2> { typedef unsigned long least; };
    template<> struct uint_least_helper<3> { typedef unsigned int least; };
    template<> struct uint_least_helper<4> { typedef unsigned short least; };
    template<> struct uint_least_helper<5> { typedef unsigned char least; };
#else
    template<> struct int_least_helper<1> { typedef long least; };
    template<> struct int_least_helper<2> { typedef int least; };
    template<> struct int_least_helper<3> { typedef short least; };
    template<> struct int_least_helper<4> { typedef signed char least; };
    
    template<> struct uint_least_helper<1> { typedef unsigned long least; };
    template<> struct uint_least_helper<2> { typedef unsigned int least; };
    template<> struct uint_least_helper<3> { typedef unsigned short least; };
    template<> struct uint_least_helper<4> { typedef unsigned char least; };
#endif

  }; // namespace detail


  /// A signed integer type specified by the number of bits (including sign)
  template< int Bits > 
  struct int_t {
    typedef typename detail::int_least_helper
    <
#if HAVE_LONG_LONG
    (Bits-1 <= 8 * sizeof(long long)) +
#elif HAVE_MSVC_INT64
    (Bits-1 <= 8 * sizeof(__int64)) +
#endif
    (Bits-1 <= 8 * sizeof(long)) +
    (Bits-1 <= 8 * sizeof(int)) +
    (Bits-1 <= 8 * sizeof(short)) +
    (Bits-1 <= 8 * sizeof(signed char))
      >::least  least;
    typedef typename detail::int_fast_t<least>::fast  fast;
  };

  /// An unsigned integer type specified by the number of bits 
  template< int Bits > 
  struct uint_t {
    typedef typename detail::uint_least_helper
    < 
#if HAVE_LONG_LONG
    (Bits-1 <= 8 * sizeof(unsigned long long)) +
#elif HAVE_MSVC_INT64
    (Bits-1 <= 8 * sizeof(unsigned __int64)) +
#endif
    (Bits-1 <= 8 * sizeof(unsigned long)) +
    (Bits-1 <= 8 * sizeof(unsigned int)) +
    (Bits-1 <= 8 * sizeof(unsigned short)) +
    (Bits-1 <= 8 * sizeof(unsigned char))
      >::least  least;
    typedef typename detail::int_fast_t<least>::fast  fast;
  };

} // namespace sl

#endif  // SL_INTEGER_HPP

