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
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                                                                           -
//                       ****************************                        -
//                        ARITHMETIC CODING EXAMPLES                         -
//                       ****************************                        -
//                                                                           -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                                                                           -
// Fast arithmetic coding implementation                                     -
// -> 32-bit variables, 32-bit product, periodic updates, table decoding     -
//                                                                           -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                                                                           -
// Version 1.00  -  April 25, 2004                                           -
//                                                                           -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                                                                           -
//                                  WARNING                                  -
//                                 =========                                 -
//                                                                           -
// The only purpose of this program is to demonstrate the basic principles   -
// of arithmetic coding. It is provided as is, without any express or        -
// implied warranty, without even the warranty of fitness for any particular -
// purpose, or that the implementations are correct.                         -
//                                                                           -
// Permission to copy and redistribute this code is hereby granted, provided -
// that this warning and copyright notices are not removed or altered.       -
//                                                                           -
// Copyright (c) 2004 by Amir Said (said@ieee.org) &                         -
//                       William A. Pearlman (pearlw@ecse.rpi.edu)           -
//                                                                           -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                                                                           -
// A description of the arithmetic coding method used here is available in   -
//                                                                           -
// Lossless Compression Handbook, ed. K. Sayood                              -
// Chapter 5: Arithmetic Coding (A. Said), pp. 101-152, Academic Press, 2003 -
//                                                                           -
// A. Said, Introduction to Arithetic Coding Theory and Practice             -
// HP Labs report HPL-2004-76  -  http://www.hpl.hp.com/techreports/         -
//                                                                           -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#include "sl/assert.hpp"
#include "sl/arithmetic_codec.hpp"
#include <cstdlib>

namespace sl {

  const unsigned AC__MinLength = 0x01000000U;   // threshold for renormalization
  const unsigned AC__MaxLength = 0xFFFFFFFFU;      // maximum AC interval length

  // Maximum values for binary models
  const unsigned BM__LengthShift = 13;     // length bits discarded before mult.
  const unsigned BM__MaxCount    = 1 << BM__LengthShift;  // for adaptive models

  // Maximum values for general models
  const unsigned DM__LengthShift = 15;     // length bits discarded before mult.
  const unsigned DM__MaxCount    = 1 << DM__LengthShift;  // for adaptive models

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void arithmetic_codec::propagate_carry() {
    unsigned char * p;            // carry propagation on compressed data buffer
    for (p = ac_pointer_ - 1; *p == 0xFFU; p--) {
      assert(current_encode_byte_count() < buffer_size_);
      *p = 0;
    }
    assert(current_encode_byte_count() < buffer_size_);
    ++*p;
  }
  
  void arithmetic_codec::renorm_enc_interval() {
    do {                                          // output and discard top byte
      assert(current_encode_byte_count() < buffer_size_);
      *ac_pointer_++ = (unsigned char)(base_ >> 24);
      base_ <<= 8;
    } while ((length_ <<= 8) < AC__MinLength);        // length multiplied by 256
  }
  
  void arithmetic_codec::renorm_dec_interval() {
    do {                                          // read least-significant byte
      assert(unsigned(ac_pointer_+1 - code_buffer_) < buffer_size_+4);
      value_ = (value_ << 8) | unsigned(*++ac_pointer_);
    } while ((length_ <<= 8) < AC__MinLength);        // length multiplied by 256
  }

  void arithmetic_codec::put_bit(unsigned bit) {
#ifndef NDEBUG
    if (mode_ != 1) SL_FAIL("encoder not initialized");
#endif

    length_ >>= 1;                                              // halve interval
    if (bit) {
      unsigned init_base = base_;
      base_ += length_;                                               // move base
      if (init_base > base_) propagate_carry();               // overflow = carry
    }

    if (length_ < AC__MinLength) renorm_enc_interval();        // renormalization
  }

  unsigned arithmetic_codec::get_bit() {
#ifndef NDEBUG
    if (mode_ != 2) SL_FAIL("decoder not initialized");
#endif

    length_ >>= 1;                                              // halve interval
    unsigned bit = (value_ >= length_);                              // decode bit
    if (bit) value_ -= length_;                                       // move base

    if (length_ < AC__MinLength) renorm_dec_interval();        // renormalization

    return bit;                                         // return data bit value
  }

  void arithmetic_codec::put_bits(unsigned data, unsigned bits) {
#ifndef NDEBUG
    if (mode_ != 1) SL_FAIL("encoder not initialized");
    if ((bits < 1) || (bits > 20)) SL_FAIL("invalid number of bits");
    if (data >= (1U << bits)) SL_FAIL("invalid data");
#endif

    unsigned init_base = base_;
    base_ += data * (length_ >>= bits);            // new interval base and length

    if (init_base > base_) propagate_carry();                 // overflow = carry
    if (length_ < AC__MinLength) renorm_enc_interval();        // renormalization
  }

  unsigned arithmetic_codec::get_bits(unsigned bits) {
#ifndef NDEBUG
    if (mode_ != 2) SL_FAIL("decoder not initialized");
    if ((bits < 1) || (bits > 20)) SL_FAIL("invalid number of bits");
#endif

    unsigned s = value_ / (length_ >>= bits);      // decode symbol, change length

    value_ -= length_ * s;                                      // update interval
    if (length_ < AC__MinLength) renorm_dec_interval();        // renormalization

    return s;
  }

  void arithmetic_codec::encode(unsigned bit,
                                static_bit_model & M) {
#ifndef NDEBUG
    if (mode_ != 1) SL_FAIL("encoder not initialized");
#endif

    unsigned x = M.bit_0_prob * (length_ >> BM__LengthShift);   // product l x p0
    // update interval
    if (bit == 0)
      length_  = x;
    else {
      unsigned init_base = base_;
      base_   += x;
      length_ -= x;
      if (init_base > base_) propagate_carry();               // overflow = carry
    }

    if (length_ < AC__MinLength) renorm_enc_interval();        // renormalization
  }

  unsigned arithmetic_codec::decode(static_bit_model & M) {
#ifndef NDEBUG
    if (mode_ != 2) SL_FAIL("decoder not initialized");
#endif

    unsigned x = M.bit_0_prob * (length_ >> BM__LengthShift);   // product l x p0
    unsigned bit = (value_ >= x);                                     // decision
    // update & shift interval
    if (bit == 0)
      length_  = x;
    else {
      value_  -= x;                                 // shifted interval base = 0
      length_ -= x;
    }

    if (length_ < AC__MinLength) renorm_dec_interval();        // renormalization

    return bit;                                         // return data bit value
  }

  void arithmetic_codec::encode(unsigned bit,
                                adaptive_bit_model & M) {
#ifndef NDEBUG
    if (mode_ != 1) SL_FAIL("encoder not initialized");
#endif

    unsigned x = M.bit_0_prob * (length_ >> BM__LengthShift);   // product l x p0
    // update interval
    if (bit == 0) {
      length_ = x;
      ++M.bit_0_count;
    }
    else {
      unsigned init_base = base_;
      base_   += x;
      length_ -= x;
      if (init_base > base_) propagate_carry();               // overflow = carry
    }

    if (length_ < AC__MinLength) renorm_enc_interval();        // renormalization

    if (--M.bits_until_update == 0) M.update();         // periodic model update
  }

  unsigned arithmetic_codec::decode(adaptive_bit_model & M) {
#ifndef NDEBUG
    if (mode_ != 2) SL_FAIL("decoder not initialized");
#endif

    unsigned x = M.bit_0_prob * (length_ >> BM__LengthShift);   // product l x p0
    unsigned bit = (value_ >= x);                                     // decision
    // update interval
    if (bit == 0) {
      length_ = x;
      ++M.bit_0_count;
    }
    else {
      value_  -= x;
      length_ -= x;
    }

    if (length_ < AC__MinLength) renorm_dec_interval();        // renormalization

    if (--M.bits_until_update == 0) M.update();         // periodic model update

    return bit;                                         // return data bit value
  }

  void arithmetic_codec::encode(unsigned data,
                                static_data_model & M) {
#ifndef NDEBUG
    if (mode_ != 1) SL_FAIL("encoder not initialized");
    if (data >= M.data_symbols) SL_FAIL("invalid data symbol");
#endif

    unsigned x, init_base = base_;
    // compute products
    if (data == M.last_symbol) {
      x = M.distribution[data] * (length_ >> DM__LengthShift);
      base_   += x;                                            // update interval
      length_ -= x;                                          // no product needed
    }
    else {
      x = M.distribution[data] * (length_ >>= DM__LengthShift);
      base_   += x;                                            // update interval
      length_  = M.distribution[data+1] * length_ - x;
    }
             
    if (init_base > base_) propagate_carry();                 // overflow = carry

    if (length_ < AC__MinLength) renorm_enc_interval();        // renormalization
  }

  unsigned arithmetic_codec::decode(static_data_model & M) {
#ifndef NDEBUG
    if (mode_ != 2) SL_FAIL("decoder not initialized");
#endif

    unsigned n, s, x, y = length_;

    if (M.decoder_table) {              // use table look-up for faster decoding

      unsigned dv = value_ / (length_ >>= DM__LengthShift);
      unsigned t = dv >> M.table_shift;

      s = M.decoder_table[t];         // initial decision based on table look-up
      n = M.decoder_table[t+1] + 1;

      while (n > s + 1) {                        // finish with bisection search
        unsigned m = (s + n) >> 1;
        if (M.distribution[m] > dv) n = m; else s = m;
      }
      // compute products
      x = M.distribution[s] * length_;
      if (s != M.last_symbol) y = M.distribution[s+1] * length_;
    }

    else {                                  // decode using only multiplications

      x = s = 0;
      length_ >>= DM__LengthShift;
      unsigned m = (n = M.data_symbols) >> 1;
      // decode via bisection search
      do {
        unsigned z = length_ * M.distribution[m];
        if (z > value_) {
          n = m;
          y = z;                                             // value is smaller
        }
        else {
          s = m;
          x = z;                                     // value is larger or equal
        }
      } while ((m = (s + n) >> 1) != s);
    }

    value_ -= x;                                               // update interval
    length_ = y - x;

    if (length_ < AC__MinLength) renorm_dec_interval();        // renormalization

    return s;
  }

  void arithmetic_codec::encode(unsigned data,
                                adaptive_data_model & M) {
#ifndef NDEBUG
    if (mode_ != 1) SL_FAIL("encoder not initialized");
    if (data >= M.data_symbols) SL_FAIL("invalid data symbol");
#endif

    unsigned x, init_base = base_;
    // compute products
    if (data == M.last_symbol) {
      x = M.distribution[data] * (length_ >> DM__LengthShift);
      base_   += x;                                            // update interval
      length_ -= x;                                          // no product needed
    }
    else {
      x = M.distribution[data] * (length_ >>= DM__LengthShift);
      base_   += x;                                            // update interval
      length_  = M.distribution[data+1] * length_ - x;
    }

    if (init_base > base_) propagate_carry();                 // overflow = carry

    if (length_ < AC__MinLength) renorm_enc_interval();        // renormalization

    ++M.symbol_count[data];
    if (--M.symbols_until_update == 0) M.update(true);  // periodic model update
  }

  unsigned arithmetic_codec::decode(adaptive_data_model & M) {
#ifndef NDEBUG
    if (mode_ != 2) SL_FAIL("decoder not initialized");
#endif

    unsigned n, s, x, y = length_;

    if (M.decoder_table) {              // use table look-up for faster decoding

      unsigned dv = value_ / (length_ >>= DM__LengthShift);
      unsigned t = dv >> M.table_shift;

      s = M.decoder_table[t];         // initial decision based on table look-up
      n = M.decoder_table[t+1] + 1;

      while (n > s + 1) {                        // finish with bisection search
        unsigned m = (s + n) >> 1;
        if (M.distribution[m] > dv) n = m; else s = m;
      }
      // compute products
      x = M.distribution[s] * length_;
      if (s != M.last_symbol) y = M.distribution[s+1] * length_;
    }

    else {                                  // decode using only multiplications

      x = s = 0;
      length_ >>= DM__LengthShift;
      unsigned m = (n = M.data_symbols) >> 1;
      // decode via bisection search
      do {
        unsigned z = length_ * M.distribution[m];
        if (z > value_) {
          n = m;
          y = z;                                             // value is smaller
        }
        else {
          s = m;
          x = z;                                     // value is larger or equal
        }
      } while ((m = (s + n) >> 1) != s);
    }

    value_ -= x;                                               // update interval
    length_ = y - x;

    if (length_ < AC__MinLength) renorm_dec_interval();        // renormalization

    ++M.symbol_count[s];
    if (--M.symbols_until_update == 0) M.update(false);  // periodic model update

    return s;
  }


  arithmetic_codec::arithmetic_codec() {
    mode_ = buffer_size_ = 0;
    new_buffer_ = code_buffer_ = 0;
  }

  arithmetic_codec::arithmetic_codec(unsigned max_code_bytes,
                                     unsigned char * user_buffer) {
    mode_ = buffer_size_ = 0;
    new_buffer_ = code_buffer_ = 0;
    set_buffer(max_code_bytes, user_buffer);
  }

  arithmetic_codec::~arithmetic_codec() {
    if (new_buffer_) delete [] new_buffer_;
    new_buffer_ = 0;
  }

  void arithmetic_codec::set_buffer(unsigned max_code_bytes,
                                    unsigned char * user_buffer) {
    // test for reasonable sizes
    if (max_code_bytes > 0x1000000U)
      SL_FAIL("invalid codec buffer size");
    if (mode_ != 0) SL_FAIL("cannot set buffer while encoding or decoding");

    if (user_buffer != 0) {                       // user provides memory buffer
      buffer_size_ = max_code_bytes;
      code_buffer_ = user_buffer;               // set buffer for compressed data
      if (new_buffer_) delete [] new_buffer_;                 // free anything previously assigned
      new_buffer_ = 0;
    } else if (max_code_bytes <= buffer_size_) {
      // enough available
    } else {
      buffer_size_ = max_code_bytes;                           // assign new memory
      if (new_buffer_) delete [] new_buffer_;                   // free anything previously assigned
      if ((new_buffer_ = new unsigned char[buffer_size_+16]) == 0) // 16 extra bytes
        SL_FAIL("cannot assign memory for compressed data buffer");
      code_buffer_ = new_buffer_;                  // set buffer for compressed data
    }
  }

  void arithmetic_codec::start_encoder() {
    if (mode_ != 0) SL_FAIL("cannot start encoder");
    if (buffer_size_ == 0) SL_FAIL("no code buffer set");
    if (buffer_size_ <16) SL_FAIL("buffer too small");
    
    mode_   = 1;
    base_   = 0;            // initialize encoder variables: interval and pointer
    length_ = AC__MaxLength;
    ac_pointer_ = code_buffer_;                       // pointer to next data byte
  }

  void arithmetic_codec::start_decoder() {
    if (mode_ != 0) SL_FAIL("cannot start decoder");
    if (buffer_size_ == 0) SL_FAIL("no code buffer set");

    // initialize decoder: interval, pointer, initial code value
    mode_   = 2;
    length_ = AC__MaxLength;
      
    ac_pointer_ = code_buffer_ + 3;
    value_ = (unsigned(code_buffer_[0]) << 24)|(unsigned(code_buffer_[1]) << 16) |
      (unsigned(code_buffer_[2]) <<  8)| unsigned(code_buffer_[3]);
  }

  unsigned arithmetic_codec::stop_encoder() {
    if (mode_ != 1) SL_FAIL("invalid to stop encoder");
    mode_ = 0;

    unsigned init_base = base_;            // done encoding: set final data bytes

    if (length_ > 2 * AC__MinLength) {
      base_  += AC__MinLength;                                     // base offset
      length_ = AC__MinLength >> 1;             // set new length for 1 more byte
    }
    else {
      base_  += AC__MinLength >> 1;                                // base offset
      length_ = AC__MinLength >> 9;            // set new length for 2 more bytes
    }

    if (init_base > base_) propagate_carry();                 // overflow = carry

    renorm_enc_interval();                // renormalization = output last bytes

    unsigned code_bytes = unsigned(ac_pointer_ - code_buffer_);
    while (code_bytes<4) { // FIXME: to ensure start_decoding works
      ++ac_pointer_;
      ++code_bytes;
    }
    
    if (code_bytes > buffer_size_) SL_FAIL("code buffer overflow");

    return code_bytes;                                   // number of bytes used
  }

  void arithmetic_codec::stop_decoder() {
    if (mode_ != 2) SL_FAIL("invalid to stop decoder");
    mode_ = 0;
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - Static bit model implementation - - - - - - - - - - - - - - - - - - - - -

  static_bit_model::static_bit_model() {
    bit_0_prob = 1U << (BM__LengthShift - 1);                        // p0 = 0.5
  }

  void static_bit_model::set_probability_0(double p0) {
    if ((p0 < 0.0001)||(p0 > 0.9999)) SL_FAIL("invalid bit probability");
    bit_0_prob = unsigned(p0 * (1 << BM__LengthShift));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - Adaptive bit model implementation - - - - - - - - - - - - - - - - - - - -

  adaptive_bit_model::adaptive_bit_model() {
    reset();
  }

  void adaptive_bit_model::reset() {
    // initialization to equiprobable model
    bit_0_count = 1;
    bit_count   = 2;
    bit_0_prob  = 1U << (BM__LengthShift - 1);
    update_cycle = bits_until_update = 4;         // start with frequent updates
  }

  void adaptive_bit_model::update() {
    // halve counts when a threshold is reached

    if ((bit_count += update_cycle) > BM__MaxCount) {
      bit_count = (bit_count + 1) >> 1;
      bit_0_count = (bit_0_count + 1) >> 1;
      if (bit_0_count == bit_count) ++bit_count;
    }
    // compute scaled bit 0 probability
    unsigned scale = 0x80000000U / bit_count;
    bit_0_prob = (bit_0_count * scale) >> (31 - BM__LengthShift);

    // set frequency of model updates
    update_cycle = (5 * update_cycle) >> 2;
    if (update_cycle > 64) update_cycle = 64;
    bits_until_update = update_cycle;
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - Static data model implementation  - - - - - - - - - - - - - - - - - - -

  static_data_model::static_data_model() {
    data_symbols = 0;
    distribution = 0;
  }

  static_data_model::~static_data_model() {
    delete [] distribution;
  }

  void static_data_model::set_distribution(unsigned number_of_symbols,
                                           const double probability[]) {
    if ((number_of_symbols < 2) || (number_of_symbols > (1 << 11)))
      SL_FAIL("invalid number of data symbols");

    if (data_symbols != number_of_symbols) {     // assign memory for data model
      data_symbols = number_of_symbols;
      last_symbol = data_symbols - 1;
      delete [] distribution;
      // define size of table for fast decoding
      if (data_symbols > 16) {
        unsigned table_bits = 3;
        while (data_symbols > (1U << (table_bits + 2))) ++table_bits;
        table_size  = 1 << table_bits;
        table_shift = DM__LengthShift - table_bits;
        distribution = new unsigned[data_symbols+table_size+2];
        decoder_table = distribution + data_symbols;
      }
      else {                                  // small alphabet: no table needed
        decoder_table = 0;
        table_size = table_shift = 0;
        distribution = new unsigned[data_symbols];
      }
      if (distribution == 0) SL_FAIL("cannot assign model memory");
    }
    // compute cumulative distribution, decoder table
    unsigned s = 0;
    double sum = 0.0, p = 1.0 / double(data_symbols);

    for (unsigned k = 0; k < data_symbols; k++) {
      if (probability) p = probability[k];
      if ((p < 0.0001) || (p > 0.9999)) SL_FAIL("invalid symbol probability");
      distribution[k] = unsigned(sum * (1 << DM__LengthShift));
      sum += p;
      if (table_size == 0) continue;
      unsigned w = distribution[k] >> table_shift;
      while (s < w) decoder_table[++s] = k - 1;
    }

    if (table_size != 0) {
      decoder_table[0] = 0;
      while (s <= table_size) decoder_table[++s] = data_symbols - 1;
    }

    if ((sum < 0.9999) || (sum > 1.0001)) SL_FAIL("invalid probabilities");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - Adaptive data model implementation  - - - - - - - - - - - - - - - - - -

  adaptive_data_model::adaptive_data_model() {
    data_symbols = 0;
    distribution = 0;
  }

  adaptive_data_model::adaptive_data_model(unsigned number_of_symbols) {
    data_symbols = 0;
    distribution = 0;
    set_alphabet(number_of_symbols);
  }

  adaptive_data_model::~adaptive_data_model() {
    delete [] distribution;
  }

  void adaptive_data_model::set_alphabet(unsigned number_of_symbols) {
    if ((number_of_symbols < 2) || (number_of_symbols > (1 << 11)))
      SL_FAIL("invalid number of data symbols");

    if (data_symbols != number_of_symbols) {     // assign memory for data model
      data_symbols = number_of_symbols;
      last_symbol = data_symbols - 1;
      delete [] distribution;
      // define size of table for fast decoding
      if (data_symbols > 16) {
        unsigned table_bits = 3;
        while (data_symbols > (1U << (table_bits + 2))) ++table_bits;
        table_size  = 1 << table_bits;
        table_shift = DM__LengthShift - table_bits;
        distribution = new unsigned[2*data_symbols+table_size+2];
        decoder_table = distribution + 2 * data_symbols;
      }
      else {                                  // small alphabet: no table needed
        decoder_table = 0;
        table_size = table_shift = 0;
        distribution = new unsigned[2*data_symbols];
      }
      symbol_count = distribution + data_symbols;
      if (distribution == 0) SL_FAIL("cannot assign model memory");
    }

    reset();                                                 // initialize model
  }

  void adaptive_data_model::update(bool from_encoder) {
    // halve counts when a threshold is reached

    if ((total_count += update_cycle) > DM__MaxCount) {
      total_count = 0;
      for (unsigned n = 0; n < data_symbols; n++)
        total_count += (symbol_count[n] = (symbol_count[n] + 1) >> 1);
    }
    // compute cumulative distribution, decoder table
    unsigned k, sum = 0, s = 0;
    unsigned scale = 0x80000000U / total_count;

    if (from_encoder || (table_size == 0))
      for (k = 0; k < data_symbols; k++) {
        distribution[k] = (scale * sum) >> (31 - DM__LengthShift);
        sum += symbol_count[k];
      }
    else {
      for (k = 0; k < data_symbols; k++) {
        distribution[k] = (scale * sum) >> (31 - DM__LengthShift);
        sum += symbol_count[k];
        unsigned w = distribution[k] >> table_shift;
        while (s < w) decoder_table[++s] = k - 1;
      }
      decoder_table[0] = 0;
      while (s <= table_size) decoder_table[++s] = data_symbols - 1;
    }
    // set frequency of model updates
    update_cycle = (5 * update_cycle) >> 2;
    unsigned max_cycle = (data_symbols + 6) << 3;
    if (update_cycle > max_cycle) update_cycle = max_cycle;
    symbols_until_update = update_cycle;
  }

  void adaptive_data_model::reset() {
    if (data_symbols == 0) return;

    // restore probability estimates to uniform distribution
    total_count = 0;
    update_cycle = data_symbols;
    for (unsigned k = 0; k < data_symbols; k++) symbol_count[k] = 1;
    update(false);
    symbols_until_update = update_cycle = (data_symbols + 6) >> 1;
  }

} // namespace sl
