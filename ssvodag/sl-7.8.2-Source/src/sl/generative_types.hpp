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
#ifndef SL_GENERATIVE_TYPES_HPP
#define SL_GENERATIVE_TYPES_HPP

#include <sl/config.hpp>

namespace sl {

  /**
   *  An empty struct representing the value true, for use
   *  in template specializations.
   */
  struct true_t { enum { value = true }; }; 
  
  /**
   *  An empty struct representing the value false, for use
   *  in template specializations.
   */
  struct false_t { enum { value = false }; };
    
  /**
   *  A class representing the boolean function not
   */
  template<bool b>
  struct gen_not {
    enum { value = true };
  };

  template<>
  struct gen_not<true> {
    enum { value = false };
  };

  /**
   *  A class representing the boolean function and
   */
  template<bool b1, bool b2>
  struct gen_and {
    enum { value = false };
  };

  template<>
  struct gen_and<true, true> {
    enum { value = true };
  };

  /**
   *  A class representing the boolean function or
   */
  template<bool b1, bool b2>
  struct gen_or {
    enum { value = true };
  };

  template<>
  struct gen_or<false, false> {
    enum { value = false };
  };
  
    
  /**
   *  A class representing the boolean condition n1==n2
   */
  template<bool n1, bool n2>
  struct gen_equal {
    enum { value = (n1==n2) };
  };

  /**
   *  A class representing the boolean condition n1!=n2
   */
  template<bool n1, bool n2>
  struct gen_not_equal {
    enum { value = (n1!=n2) };
  };

  /**
   *  Type selection based on boolean condition. Based on ideas 
   *  popularized by 
   *  the gcl library by Tobias Neubert, Krzysztof Czarnecki, and 
   *  Ulrich Eisenecker (1998).
   */
  template <bool cond, class T_iftrue, class T_iffalse>
  struct gen_if { };

  template <class T_iftrue, class T_iffalse>
  struct gen_if<true, T_iftrue, T_iffalse> { typedef T_iftrue type; };

  template <class T_iftrue, class T_iffalse>
  struct gen_if<false, T_iftrue, T_iffalse> { typedef T_iffalse type; };

}; // namespace sl
 
#endif /* SL_GENERATIVE_TYPES_HPP */
