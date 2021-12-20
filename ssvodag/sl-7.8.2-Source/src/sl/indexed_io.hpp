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
#ifndef SL_INDEXED_IO_HPP
#define SL_INDEXED_IO_HPP

#include <sl/indexed.hpp>
#include <iostream>
#include <iomanip>

/// Write rhs to stream s
template <
  class  G_numtype, 
  size_t G_rank, 
  class  G_derived, 
  class  G_userdefined,
  class  G_discriminant
>
std::ostream& operator<<(std::ostream& s, 
			 const sl::indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>& rhs) {
  for (typename sl::indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>::const_iterator it= rhs.begin();
       !it.off();
       ++it) {
    s << it.index() 
      << " : " 
      << std::setw(6) << std::setprecision(3) << it.value() 
      << std::endl;
  }
  return s;
}

#endif
