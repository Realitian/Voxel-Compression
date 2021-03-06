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
#ifndef SL_CONV_TO_HPP 
#define SL_CONV_TO_HPP

#include <sl/config.hpp>

namespace sl {

  /**
   *  Generic conversion -- template specializations
   *  will define ways to create objects of one type
   *  from objects of another using the "from" function.
   */
  template <typename OUT_T>
  class conv_to {
  public:
    typedef OUT_T result_t;
  };
  
} // namespace sl

#endif
