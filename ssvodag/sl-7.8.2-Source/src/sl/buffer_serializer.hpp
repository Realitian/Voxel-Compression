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
#ifndef SL_BUFFER_SERIALIZER_HPP
#define SL_BUFFER_SERIALIZER_HPP

#include <sl/serializer.hpp>
#include <vector>
#include <string.h> // memcpy

namespace sl {

  /**
   *  Objects that serialize other objects to a memory buffer
   *  using a standard binary encoding
   */
  class output_buffer_serializer: public output_serializer {
  protected: 
    std::vector<uint8_t> buffer_;

  protected:
    
    inline void write_raw(const void* ptr, std::size_t sz) {
      SL_REQUIRE("ptr exists", ptr);
      std::size_t N = buffer_.size();
      buffer_.resize(N+sz);
      memcpy(&(buffer_[N]), ptr, sz);
    }
    
  public:
    output_buffer_serializer(uint32_t v=0)
        : output_serializer(v) {
    }
    
    virtual ~output_buffer_serializer() {
    }
    
    void clear() {
      buffer_.clear();
    }

    const std::vector<uint8_t>& buffer() const {
      return buffer_;
    }
    const uint8_t* buffer_address() const {
      return &(buffer_[0]);
    }
    std::size_t buffer_size() const {
      return buffer_.size();
    }

    std::vector<uint8_t>& buffer() {
      return buffer_;
    }
    
    uint8_t* buffer_address() {
      return &(buffer_[0]);
    }
    
  protected:
#if (HAVE_LONG_LONG_BITS >= 64)
    inline void write_simple(const int64_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const uint64_t x) { write_raw(&x, sizeof(x)); }
#endif
    inline void write_simple(const int32_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const uint32_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const int16_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const uint16_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const int8_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const uint8_t x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const float x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const double x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const long double x) { write_raw(&x, sizeof(x)); }
    inline void write_simple(const bool x) { write_raw(&x, sizeof(x)); }

#if (HAVE_LONG_LONG_BITS >= 64)
    inline void write_array(std::size_t n, const int64_t* x)  { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const uint64_t* x) { write_raw(x, n*sizeof(*x)); }
#endif
    inline void write_array(std::size_t n, const int32_t* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const uint32_t* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const int16_t* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const uint16_t* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const int8_t* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const uint8_t* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const float* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const double* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const long double* x) { write_raw(x, n*sizeof(*x)); }
    inline void write_array(std::size_t n, const bool* x) { write_raw(x, n*sizeof(*x)); }
  };

  /**
   *  Objects that deserialize other objects from a memory buffer
   *  using a standard binary encoding
   */
  class input_buffer_serializer: public input_serializer {
  protected: 
    std::vector<uint8_t> buffer_;
    std::size_t          cursor_;
    
  protected:
    
    inline void read_raw(void* ptr, std::size_t sz) {
      SL_REQUIRE("Ptr exists", ptr);
      SL_REQUIRE("Not off", cursor_+sz <= buffer_.size());
      memcpy(ptr, &(buffer_[cursor_]), sz);
      cursor_+= sz;
    }

  public:
    input_buffer_serializer(uint32_t v=0)
        : input_serializer(v), cursor_(0) {
    }
    
    virtual ~input_buffer_serializer() {
    }
    
    void clear() {
      buffer_.clear();
      cursor_ = 0;
    }

    void reset() {
      cursor_ = 0;
    }

    bool off() const {
      return cursor_ >= buffer_.size();
    }
    
    const std::vector<uint8_t>& buffer() const {
      return buffer_;
    }
    
    const uint8_t* buffer_address() const {
      return &(buffer_[0]);
    }
    
    const uint8_t* next_read_address() const {
      return &(buffer_[cursor_]);
    }

    std::size_t buffer_size() const {
      return buffer_.size();
    }

    std::vector<uint8_t>& buffer() {
      return buffer_;
    }
    
    uint8_t* buffer_address() {
      return &(buffer_[0]);
    }

    uint8_t* next_read_address() {
      return &(buffer_[cursor_]);
    }
    
  protected:
#if (HAVE_LONG_LONG_BITS >= 64)
    inline void read_simple(int64_t& x)  { read_raw(&x, sizeof(x)); }
    inline void read_simple(uint64_t& x) { read_raw(&x, sizeof(x)); }
#endif
    inline void read_simple(int32_t& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(uint32_t& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(int16_t& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(uint16_t& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(int8_t& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(uint8_t& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(float& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(double& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(long double& x) { read_raw(&x, sizeof(x)); }
    inline void read_simple(bool& x) { read_raw(&x, sizeof(x)); }

#if (HAVE_LONG_LONG_BITS >= 64)
    inline void read_array(std::size_t n, int64_t* x)  { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, uint64_t* x) { read_raw(x, n*sizeof(*x)); }
#endif
    inline void read_array(std::size_t n, int32_t* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, uint32_t* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, int16_t* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, uint16_t* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, int8_t* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, uint8_t* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, float* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, double* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, long double* x) { read_raw(x, n*sizeof(*x)); }
    inline void read_array(std::size_t n, bool* x) { read_raw(x, n*sizeof(*x)); }
  };
  
} // namespace sl

#endif
