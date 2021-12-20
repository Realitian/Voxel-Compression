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
#include <sl/tester.hpp>
#include <sl/external_array.hpp>
#include <sl/sort.hpp>
#include <iostream>

static std::size_t failed_test_count = 0;

std::string to_string(const sl::external_array1<float>& xarray) {
  std::ostringstream os;
  os << "XARRAY[SZ=" << xarray.size() << "]:= {";
  for (std::size_t i=0; i<xarray.size(); ++i) {
    os << " " << xarray[i];
    if (i==5 && (xarray.size()-5>i)) {
      os << " ... ";
      i = (std::size_t)(xarray.size())-5;
    }
  }
  os << " }";
  return os.str();
}

void test_external_array() {

  sl::tester tester("external_array");
  
  sl::external_array1<float>* xarray = new sl::external_array1<float>("test.arr","w");
  xarray->clear();

  tester.test("Creation", to_string(*xarray), std::string("XARRAY[SZ=0]:= { }"));
  xarray->push_back(1.0);
  xarray->push_back(2.0);
  xarray->push_back(3.0);
  tester.test("Push 1/2/3", to_string(*xarray), std::string("XARRAY[SZ=3]:= { 1 2 3 }"));
  delete xarray;
  xarray = new sl::external_array1<float>("test.arr","a");

  tester.test("Destroy + reopen", to_string(*xarray), std::string("XARRAY[SZ=3]:= { 1 2 3 }"));
  (*xarray)[0] = -1.0;
  tester.test("Write at 0", to_string(*xarray), std::string("XARRAY[SZ=3]:= { -1 2 3 }"));
  xarray->clear();
  tester.test("Clear", to_string(*xarray), std::string("XARRAY[SZ=0]:= { }"));

  for (std::size_t i=1; i<=10000; ++i) {
    xarray->push_back((float)(i)*-1.0f);
  }

  tester.test("10000 push", to_string(*xarray), std::string("XARRAY[SZ=10000]:= { -1 -2 -3 -4 -5 -6 ...  -9997 -9998 -9999 -10000 }"));

  std::sort(xarray->begin(), xarray->end());
  tester.test("Standard sorting", to_string(*xarray), std::string("XARRAY[SZ=10000]:= { -10000 -9999 -9998 -9997 -9996 -9995 ...  -4 -3 -2 -1 }"));

  xarray->clear();
  tester.test("Clear", to_string(*xarray), std::string("XARRAY[SZ=0]:= { }"));

  for (std::size_t i=1; i<=10000; ++i) {
    xarray->push_back((float)(i)*-1.0f);
  }
  tester.test("10000 push (2)", to_string(*xarray), std::string("XARRAY[SZ=10000]:= { -1 -2 -3 -4 -5 -6 ...  -9997 -9998 -9999 -10000 }"));

  std::sort(xarray->begin(), xarray->end());
  tester.test("Inplace sorting", to_string(*xarray), std::string("XARRAY[SZ=10000]:= { -10000 -9999 -9998 -9997 -9996 -9995 ...  -4 -3 -2 -1 }"));

  delete xarray;
  failed_test_count += tester.failed_test_count();
}

int main() {
  test_external_array();
  
  return failed_test_count;
}




