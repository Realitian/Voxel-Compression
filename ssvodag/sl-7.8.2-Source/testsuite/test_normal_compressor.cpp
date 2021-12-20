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
#include <sl/normal_compressor.hpp>
#include <iostream>

template < int Bits > 
class normal_compressor_tester {
public:
  static void do_it() {
    sl::normal_compressor<Bits> compressor;
    std::cerr << "Testing normal compression to " << Bits << " bits" << std::endl;
    std::cerr << "  Actual bits:       " << compressor.N << std::endl;
    std::cerr << "    Octant bits:     " << compressor.N_octant << std::endl;
    std::cerr << "    Mantissa X bits: " << compressor.N_mantissa_x << std::endl;
    std::cerr << "    Mantissa Y bits: " << compressor.N_mantissa_y << std::endl;
    double accuracy = compressor.accuracy();
    std::cerr << "  Accuracy = " << accuracy << " rad = " << accuracy*180.0/3.14159265 << " degrees" << std::endl;
  }
};

int main() {
  normal_compressor_tester<8>::do_it();
  normal_compressor_tester<16>::do_it();
  normal_compressor_tester<32>::do_it();
  return 0;
}
