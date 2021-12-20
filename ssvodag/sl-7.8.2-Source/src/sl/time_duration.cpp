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
#include <sl/time_duration.hpp>

std::ostream& operator<<(std::ostream& os, const sl::time_duration& d) {
  sl::int64_t v, h, m, s, us;
  v = d.as_microseconds();
 
  us = v>0 ? v : -v;
  s  = (us / 1000000);
  m  = (s / 60);
  h  =  m / 60;

  char c_sign = (v < 0 ? '-' : '+');
  return os << c_sign
	    << (int)(h)
	    << (int)(m % 60)
	    << (int)(s % 60)
	    << (int)(us % 1000000);
}  
