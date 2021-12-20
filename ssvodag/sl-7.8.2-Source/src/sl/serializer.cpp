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
#include <sl/serializer.hpp> 
#include <string>

sl::output_serializer::output_serializer(sl::uint32_t v) : version_(v) {}

//-----------------------------------------------------------------------
// The default write operators
//-----------------------------------------------------------------------

#define SL_IMPLEMENT_AS(A,B) \
void sl::output_serializer::write_simple(const A x) \
{ write_simple(B(x)); }

#if (HAVE_LONG_LONG_BITS >= 64)
SL_IMPLEMENT_AS(sl::int64_t,sl::int32_t)
#endif
SL_IMPLEMENT_AS(sl::int16_t,sl::int32_t)
SL_IMPLEMENT_AS(sl::int8_t,sl::int16_t)
#if (HAVE_LONG_LONG_BITS >= 64)
SL_IMPLEMENT_AS(sl::uint64_t,sl::int64_t)
#endif
SL_IMPLEMENT_AS(sl::uint32_t,sl::int32_t)
SL_IMPLEMENT_AS(sl::uint16_t,sl::int16_t)
SL_IMPLEMENT_AS(sl::uint8_t,sl::int8_t)
SL_IMPLEMENT_AS(float,double)
SL_IMPLEMENT_AS(long double,double)
SL_IMPLEMENT_AS(bool,sl::int32_t)

#undef SL_IMPLEMENT_AS

//-----------------------------------------------------------------------
// write_array
// 
// simply writes each element.
// should be implemented in an optimized way by the derived classes
//-----------------------------------------------------------------------

#define SL_SERIALIZER_STORE_ARRAY(T) \
void sl::output_serializer::write_array(std::size_t n, const T* p)  \
{ \
  for (std::size_t i=0;i<n;++i) \
    write_simple(p[i]); \
}

#if (HAVE_LONG_LONG_BITS >= 64)
SL_SERIALIZER_STORE_ARRAY(sl::int64_t)
#endif
SL_SERIALIZER_STORE_ARRAY(sl::int32_t)
SL_SERIALIZER_STORE_ARRAY(sl::int16_t)
SL_SERIALIZER_STORE_ARRAY(sl::int8_t)
#if (HAVE_LONG_LONG_BITS >= 64)
SL_SERIALIZER_STORE_ARRAY(sl::uint64_t)
#endif
SL_SERIALIZER_STORE_ARRAY(sl::uint32_t)
SL_SERIALIZER_STORE_ARRAY(sl::uint16_t)
SL_SERIALIZER_STORE_ARRAY(sl::uint8_t)
SL_SERIALIZER_STORE_ARRAY(float)
SL_SERIALIZER_STORE_ARRAY(double)
SL_SERIALIZER_STORE_ARRAY(long double)
SL_SERIALIZER_STORE_ARRAY(bool)

#undef SERIALIZER_STORE_ARRAY

void sl::output_serializer::write_string(std::size_t n, const char* p) 
{
  for (std::size_t i=0;i<n;i++)
    write_simple(sl::uint8_t(p[i]));
}

// pointers are stored in an assocatiove array and assigned numbers
void sl::output_serializer::register_object_address(void* p) {
  std::map<void*,sl::uint32_t>::iterator it = pointer_id_.find(p);
  if (it == pointer_id_.end()) {
    std::size_t N = pointer_id_.size()+1;
    pointer_id_[p] = N;
  }
}

// instead of a pointer its number is writen.
void sl::output_serializer::write_pointer(void* p) {
  std::map<void*,sl::uint32_t>::iterator it = pointer_id_.find(p);
  SL_REQUIRE("Exists", it != pointer_id_.end());
  write_simple(it->second);
}

sl::input_serializer::input_serializer(sl::uint32_t v) : version_(v) {}

//-----------------------------------------------------------------------
// operator >> for simple data types
//-----------------------------------------------------------------------

#define SL_IMPLEMENT_AS(A,B) \
void sl::input_serializer::read_simple(A& x) \
{ B xx; read_simple(xx); x = (A)(xx); }

#if (HAVE_LONG_LONG_BITS >= 64)
SL_IMPLEMENT_AS(sl::int64_t,sl::int32_t)
#endif
SL_IMPLEMENT_AS(sl::int16_t,sl::int32_t)
SL_IMPLEMENT_AS(sl::int8_t,sl::int32_t)
#if (HAVE_LONG_LONG_BITS >= 64)
SL_IMPLEMENT_AS(sl::uint64_t,sl::int64_t)
#endif
SL_IMPLEMENT_AS(sl::uint32_t,sl::int32_t)
SL_IMPLEMENT_AS(sl::uint16_t,sl::int16_t)
SL_IMPLEMENT_AS(sl::uint8_t,sl::int8_t)
SL_IMPLEMENT_AS(float,double)
SL_IMPLEMENT_AS(long double,double)

//implementation of SL_IMPLEMENT_AS(bool,int32_t)
void sl::input_serializer::read_simple(bool& x) 
{ sl::int32_t xx; read_simple(xx); x = (xx!=0); }

#undef SL_IMPLEMENT_AS

//-----------------------------------------------------------------------
// read_array
// 
// simply reads each element.
// should be implemented in an optimized way by the derived classes
//-----------------------------------------------------------------------

#define SL_SERIALIZER_RETRIEVE_ARRAY(T) void sl::input_serializer::read_array(std::size_t n, T* p) \
{ \
  for (std::size_t i=0;i<n;i++) \
    read_simple(p[i]); \
}

#if (HAVE_LONG_LONG_BITS >= 64)
SL_SERIALIZER_RETRIEVE_ARRAY(sl::int64_t)
#endif
SL_SERIALIZER_RETRIEVE_ARRAY(sl::int32_t)
SL_SERIALIZER_RETRIEVE_ARRAY(sl::int16_t)
SL_SERIALIZER_RETRIEVE_ARRAY(sl::int8_t)
#if (HAVE_LONG_LONG_BITS >= 64)
SL_SERIALIZER_RETRIEVE_ARRAY(sl::uint64_t)
#endif
SL_SERIALIZER_RETRIEVE_ARRAY(sl::uint32_t)
SL_SERIALIZER_RETRIEVE_ARRAY(sl::uint16_t)
SL_SERIALIZER_RETRIEVE_ARRAY(sl::uint8_t)
SL_SERIALIZER_RETRIEVE_ARRAY(float)
SL_SERIALIZER_RETRIEVE_ARRAY(double)
SL_SERIALIZER_RETRIEVE_ARRAY(long double)
SL_SERIALIZER_RETRIEVE_ARRAY(bool)

#undef SERIALIZER_RETRIEVE_ARRAY


void sl::input_serializer::read_string(std::size_t n, char* p) {
  sl::uint8_t c;
  for (std::size_t i=0;i<n;i++) {
    read_simple(c);
    p[i]=char(c);
  }
}

void sl::input_serializer::register_object_address(void* p) {
  pointers_.push_back(p);
}

void* sl::input_serializer::read_pointer() {
  sl::int32_t n; read_simple(n);
  SL_REQUIRE("Exists",n<sl::int32_t(pointers_.size()));
  return pointers_[n];
}
