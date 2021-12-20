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

#ifndef SL_ARITHMETIC_CODEC
#define SL_ARITHMETIC_CODEC

#include <cassert>

namespace sl {

  /// Static model for binary data
  class static_bit_model {
  public:

    static_bit_model();
    
    void set_probability_0(double);             // set probability of symbol '0'
    
  private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    unsigned bit_0_prob;
    friend class arithmetic_codec;
  };

  /// Static model for general data
  class static_data_model {
  public:

    static_data_model();
    ~static_data_model();
    
    unsigned model_symbols() { return data_symbols; }
    
    void set_distribution(unsigned number_of_symbols,
                          const double probability[] = 0);    // 0 means uniform
    
  private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    unsigned * distribution, * decoder_table;
    unsigned data_symbols, last_symbol, table_size, table_shift;
    friend class arithmetic_codec;
  };

  /// Adaptive model for binary data
  class adaptive_bit_model {
  public:
    
    adaptive_bit_model();         
    
    void reset();                             // reset to equiprobable model
    
  private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    void     update();
    unsigned update_cycle, bits_until_update;
    unsigned bit_0_prob, bit_0_count, bit_count;
    friend class arithmetic_codec;
  };

  // Adaptive model for general data
  class adaptive_data_model {
  public:
    
    adaptive_data_model();
    adaptive_data_model(unsigned number_of_symbols);
    ~adaptive_data_model();
    
    unsigned model_symbols() { return data_symbols; }
    
    void reset();                             // reset to equiprobable model
    void set_alphabet(unsigned number_of_symbols);
    
  private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    void     update(bool);
    unsigned * distribution, * symbol_count, * decoder_table;
    unsigned total_count, update_cycle, symbols_until_update;
    unsigned data_symbols, last_symbol, table_size, table_shift;
    friend class arithmetic_codec;
  };


  /**
   * Class with both the arithmetic encoder and decoder.  All compressed data is
   * saved to a memory buffer
   */
  class arithmetic_codec {
  public:

    arithmetic_codec();
    ~arithmetic_codec();
    arithmetic_codec(unsigned max_code_bytes,
                     unsigned char * user_buffer = 0);         // 0 = assign new
    
    unsigned char       *buffer() { return code_buffer_; }
    const unsigned char *buffer() const { return code_buffer_; }
    std::size_t          buffer_size() const { return buffer_size_; }
    
    void set_buffer(unsigned max_code_bytes,
                    unsigned char * user_buffer = 0);          // 0 = assign new
    
    void     start_encoder();
    void     start_decoder();
    
    unsigned stop_encoder();                 // returns number of bytes used
    void     stop_decoder();
    
    void     put_bit(unsigned bit);
    unsigned get_bit();
    
    void     put_bits(unsigned data, unsigned number_of_bits);
    unsigned get_bits(unsigned number_of_bits);
    
    void     encode(unsigned bit,
                    static_bit_model &);
    unsigned decode(static_bit_model &);
    
    void     encode(unsigned data,
                    static_data_model &);
    unsigned decode(static_data_model &);
    
    void     encode(unsigned bit,
                    adaptive_bit_model &);
    unsigned decode(adaptive_bit_model &);
    
    void     encode(unsigned data,
                    adaptive_data_model &);
    unsigned decode(adaptive_data_model &);

    inline unsigned current_encode_byte_count() const {
      return unsigned(ac_pointer_ - code_buffer_);
    }
    
  private:  //  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
    void propagate_carry();    
    void renorm_enc_interval();    
    void renorm_dec_interval();
    
    unsigned char * code_buffer_, * new_buffer_, * ac_pointer_;
    unsigned base_, value_, length_;                     // arithmetic coding state
    unsigned buffer_size_, mode_;     // mode: 0 = undef, 1 = encoder, 2 = decoder
  };

} // namespace sl

#endif
