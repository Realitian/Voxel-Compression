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
#include <sl/assert.hpp>
#include <iostream>
#include <cstdlib>

void sl::detail::failure(const char* tp,
			 const char* tag,
			 const char* expr,
			 const char* file,
			 int line,
			 const char* fname,
			 bool do_abort) {
  std::cerr << std::endl;
  if (fname) {
    std::cerr << file << ":" << line << ": " << "In function `" << fname << "':" << std::endl;
  }
  if (expr) {
    std::cerr << file << ":" << line << ":   " << tp << " check failed: " << std::endl;
    std::cerr << file << ":" << line << ":   " << "    " << tag << ": " << expr << std::endl;
  } else {
    std::cerr << file << ":" << line << ":   failure:" << std::endl;
    std::cerr << file << ":" << line << ":   " << "    " << tag << std::endl;
  }    
  if (do_abort) {
    std::cerr << "Abort." << std::endl;
    ::abort();
  } else {
    std::cerr << file << ":" << line << ":   " << "(continuing...)" << std::endl;
  }
}






