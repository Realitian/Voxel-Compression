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
#ifndef SL_CAST_HPP
#define SL_CAST_HPP

#include <sl/assert.hpp>
#include <typeinfo>

#include <limits>

#ifdef _WIN32
#undef min
#undef max
#endif

namespace sl {

  /// dynamic cast of x to type T_TARGET, throwing bad_cast if impossible
  template <class T_TARGET, class T_SOURCE>
  inline T_TARGET polymorphic_cast(T_SOURCE* x, T_TARGET* dummy = 0) {
    SL_REQUIRE("Dummy is dummy", dummy == 0);
    if (dummy) {} // avoid unused parameter warning;
    T_TARGET result = dynamic_cast<T_TARGET>(x);
    if ( result == 0 ) throw std::bad_cast();
    return result;
  }

  /// static cast of x to type T_TARGET, with a check of dynamic tupe
  template <class T_TARGET, class T_SOURCE>
  inline T_TARGET polymorphic_downcast(T_SOURCE* x, T_TARGET* dummy = 0) {
    SL_REQUIRE("Dummy is dummy", dummy == 0);
    SL_REQUIRE("Correct type", dynamic_cast<T_TARGET>(x) == x );  
    if (dummy) {} // avoid unused parameter warning;
    return static_cast<T_TARGET>(x);
  }

  /// numeric cast of arg to type T_TARGET, throwing bad_cast when overflowing or undeflowing
  template<class T_TARGET, class T_SOURCE>
  inline T_TARGET numeric_cast(T_SOURCE arg, T_TARGET* dummy = 0) {
    SL_REQUIRE("Dummy is dummy", dummy == 0);
    if (dummy) {} // avoid unused parameter warning;
    // typedefs abbreviating respective trait classes
    typedef std::numeric_limits<T_SOURCE> arg_traits;
    typedef std::numeric_limits<T_TARGET> result_traits;
 
    SL_COMPILE_TIME_CHECK("argument must be numeric", arg_traits::is_specialized);
    SL_COMPILE_TIME_CHECK("result must be numeric", result_traits::is_specialized);
 
    if( (arg < 0 && !result_traits::is_signed) ||  // loss of negative range
	(arg_traits::is_signed &&
	 arg < result_traits::min()) ||        // underflow
	arg > result_traits::max() )            // overflow
      throw std::bad_cast();
    return static_cast<T_TARGET>(arg);
  }                                                                                                                     

}; // namespace sl


#endif


