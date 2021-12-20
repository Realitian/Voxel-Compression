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
#ifndef SL_XPR_HPP
#define SL_XPR_HPP

#include <sl/generative_types.hpp>
#include <sl/utility.hpp>

namespace sl {

  /// A marker for identifying dynamic bounds
  enum { XPR_DYNAMIC_BOUND = 4294967295U };

  /// True iff sz is a dynamic bound specification
  inline static bool xpr_is_dynamic_bound(size_t sz) {
    return sz == XPR_DYNAMIC_BOUND;
  }

}

#endif
