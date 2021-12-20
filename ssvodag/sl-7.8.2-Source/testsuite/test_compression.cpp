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
/////// ALWAYS TEST IN DEBUG MODE
#if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////
#include <sl/index.hpp> //mbr: required under clang to have operator<< defined
#include <sl/tester.hpp>
#include <sl/dense_array.hpp>
#include "sl/arithmetic_codec.hpp"
#include "sl/wavelet_transform.hpp"
#include "sl/embedded_zerotree_codec.hpp"
#include "sl/micro_jpgls_array_codec.hpp"
#include "sl/hidwt_array_codec.hpp"
#include "sl/ezw_array_codec.hpp"
#include "sl/geometric_bandelet_array_codec.hpp"
#include "sl/quantized_array_codec.hpp"
#include <sl/clock.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
static std::size_t failed_test_count = 0;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void dump_matrix(const sl::dense_array<T,2,void>& m) {
  for (std::size_t i=0; i<m.extent()[0]; ++i) {
    for (std::size_t j=0; j<m.extent()[1]; ++j) {
      std::cerr << std::setw(10) << m(i,j) << ' ';
    }
    std::cerr << std::endl;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <class T>
sl::dense_array<T,2,void> make_zero_array(std::size_t height, std::size_t width) {
  typedef T value_t;
  typedef sl::dense_array<T,2,void> array_t;
  array_t result(height, width);
  for (std::size_t i=0; i<height; ++i) {
    for (std::size_t j=0; j<width; ++j) {
      result(i,j) = value_t(0);
    }
  }
  return result;
}
  
template <class T>
sl::dense_array<T,2,void> make_array(std::size_t height, std::size_t width) {
  typedef T value_t;
  typedef sl::dense_array<T,2,void> array_t;

  array_t result(height, width);
  for (std::size_t i=0; i<height; ++i) {
    for (std::size_t j=0; j<width; ++j) {
      result(i,j) = value_t(i*j+i*17%3);
    }
  }
  return result;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

std::size_t encode_data_buffer(const std::vector<sl::uint16_t>& data_buffer,
                               sl::rc_adaptive_int_model& model,
                               sl::fixed_rc_int_codec& codec) {
  codec.start_encoder();
  for (unsigned k = 0; k < data_buffer.size(); ++k) {
    codec.encode_int(model, sl::int32_t(data_buffer[k]));
  }
  codec.stop_encoder();
  // Damage buffer after end of stream
  std::size_t code_bytes = codec.current_byte_count();
  for (std::size_t i=0; i<16; ++i) {
    ((char*)(codec.buffer()))[code_bytes+i] = 13+i*17;
  }
  return code_bytes;
}


void decode_data_buffer(std::vector<sl::uint16_t>& data_buffer,
                        sl::rc_adaptive_int_model& model,
                        sl::fixed_rc_int_codec& codec) {
  codec.start_decoder();
  for (unsigned k = 0; k < data_buffer.size(); ++k)
    data_buffer[k] = sl::uint16_t(codec.decode_int(model));
  codec.stop_decoder();
}

void test_range_codec(std::size_t symbol_count) {
  sl::tester tester("range_codec with " + sl::to_string(symbol_count) + " symbols");

  for (std::size_t sz=1; sz<128; ++sz) {
    std::vector<sl::uint16_t> encoded_data(sz);
    std::vector<sl::uint16_t> decoded_data(encoded_data.size());

    std::vector<sl::uint8_t> buf(65536);
    sl::fixed_rc_int_codec codec;
    codec.set_buffer(&(buf[0]), 65536);
    sl::rc_adaptive_int_model adaptive_model;

    for (std::size_t i=0; i<encoded_data.size(); ++i) {
      encoded_data[i] = i % symbol_count;
    }
  
    adaptive_model.reset();
    std::size_t code_bytes = encode_data_buffer(encoded_data, adaptive_model, codec);
    for (std::size_t i=0; i<16; ++i) {
      ((char*)(codec.buffer()))[code_bytes+i] = 13+i*17;
    }
    
    adaptive_model.reset();
    decode_data_buffer(decoded_data, adaptive_model, codec);
    
    tester.test(std::string("Size-")+sl::to_string(sz)+": "+"decoded(encoded(equiprobable source)) == source", decoded_data == encoded_data, true);
    //tester.test(std::string("Size-")+sl::to_string(sz)+": "+"bit rate (equiprobable source)", ((8*code_bytes) / double(encoded_data.size()) < std::log(double(symbol_count))/std::log(2.0)+0.1), true);
    
    // -- 
    for (std::size_t i=0; i<encoded_data.size(); ++i) {
      encoded_data[i] = 0;
    }
    
    adaptive_model.reset();
    code_bytes = encode_data_buffer(encoded_data, adaptive_model, codec);
    
    adaptive_model.reset();
    decode_data_buffer(decoded_data, adaptive_model, codec);
    
    tester.test(std::string("Size-")+sl::to_string(sz)+": "+"decoded(encoded(zero source)) == source", decoded_data == encoded_data, true);
    //tester.test(std::string("Size-")+sl::to_string(sz)+": "+"bit rate (zero source)", ((8*code_bytes) / double(encoded_data.size()) < 1.0), true);
    
    failed_test_count += tester.failed_test_count();
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

std::size_t encode_data_buffer(const std::vector<sl::uint16_t>& data_buffer,
                               sl::adaptive_data_model& model,
                               sl::arithmetic_codec& codec) {
  codec.start_encoder();
  for (unsigned k = 0; k < data_buffer.size(); ++k) {
    codec.encode(unsigned(data_buffer[k]), model);
  }
  codec.stop_encoder();
  // Damage buffer after end of stream
  std::size_t code_bytes = codec.current_encode_byte_count();
  for (std::size_t i=0; i<16; ++i) {
    ((char*)(codec.buffer()))[code_bytes+i] = 13+i*17;
  }
  return code_bytes;
}


void decode_data_buffer(std::vector<sl::uint16_t>& data_buffer,
                        sl::adaptive_data_model& model,
                        sl::arithmetic_codec& codec) {
  codec.start_decoder();
  for (unsigned k = 0; k < data_buffer.size(); ++k)
    data_buffer[k] = sl::uint16_t(codec.decode(model));
  codec.stop_decoder();
}

void test_arithmetic_codec(std::size_t symbol_count) {
  sl::tester tester("arithmetic_codec with " + sl::to_string(symbol_count) + " symbols");

  for (std::size_t sz=1; sz<128; ++sz) {
    std::vector<sl::uint16_t> encoded_data(sz);
    std::vector<sl::uint16_t> decoded_data(encoded_data.size());

    sl::arithmetic_codec codec(1024+encoded_data.size() * 2);
    sl::adaptive_data_model adaptive_model(symbol_count);
    adaptive_model.set_alphabet(symbol_count);

    for (std::size_t i=0; i<encoded_data.size(); ++i) {
      encoded_data[i] = i % symbol_count;
    }
  
    adaptive_model.reset();
    std::size_t code_bytes = encode_data_buffer(encoded_data, adaptive_model, codec);
    for (std::size_t i=0; i<16; ++i) {
      ((char*)(codec.buffer()))[code_bytes+i] = 13+i*17;
    }
    
    adaptive_model.reset();
    decode_data_buffer(decoded_data, adaptive_model, codec);
    
    tester.test(std::string("Size-")+sl::to_string(sz)+": "+"decoded(encoded(equiprobable source)) == source", decoded_data == encoded_data, true);
    //tester.test(std::string("Size-")+sl::to_string(sz)+": "+"bit rate (equiprobable source)", ((8*code_bytes) / double(encoded_data.size()) < std::log(double(symbol_count))/std::log(2.0)+0.1), true);
    
    // -- 
    for (std::size_t i=0; i<encoded_data.size(); ++i) {
      encoded_data[i] = 0;
    }
    
    adaptive_model.reset();
    code_bytes = encode_data_buffer(encoded_data, adaptive_model, codec);
    
    adaptive_model.reset();
    decode_data_buffer(decoded_data, adaptive_model, codec);
    
    tester.test(std::string("Size-")+sl::to_string(sz)+": "+"decoded(encoded(zero source)) == source", decoded_data == encoded_data, true);
    //tester.test(std::string("Size-")+sl::to_string(sz)+": "+"bit rate (zero source)", ((8*code_bytes) / double(encoded_data.size()) < 1.0), true);
    
    failed_test_count += tester.failed_test_count();
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void test_wavelet1d(std::string wavelet_name,
                    sl::wavelet_transform<float>* xform) {
  sl::tester tester(wavelet_name + " 1D tests");

  typedef sl::dense_array<float, 1, void> float_array1_t;

  const std::size_t SZ=64;
  float_array1_t in_data(SZ);
  for (std::size_t i=0; i<SZ; ++i) {
    in_data(i) = std::sin(float(i));
  }
  float_array1_t xformed_data(SZ);
  float_array1_t in_data2(SZ);

  xform->forward1d(in_data, xformed_data, in_data.count());
  xform->backward1d(xformed_data, in_data2, in_data.count());

  tester.test("power of 2: backward1d(forward1d(signal)) == signal",
              sl::rms(in_data, in_data2), sl::intervalf("0.0000"));
  if (xform->is_orthonormal()) {
    tester.test("energy(spatial domain) == energy(wavelet domain)",
                float(sl::abs(sl::energy(in_data) - sl::energy(xformed_data))),sl::intervalf( "0.0000"));
  }

  in_data.resize(in_data.count()-1);
  xformed_data.resize(xformed_data.count()-1);
  in_data2.resize(in_data2.count()-1);
  xform->forward1d(in_data, xformed_data, in_data.count());
  xform->backward1d(xformed_data, in_data2, in_data.count());
  tester.test("non power of 2: backward1d(forward1d(signal)) ~ signal",
              sl::rms(in_data, in_data2), sl::intervalf("0.000"));
  
  failed_test_count += tester.failed_test_count();
}

void test_wavelet2d(std::string wavelet_name,
                    sl::wavelet_transform<float>* xform) {
  sl::tester tester(wavelet_name + " 2D tests");

  const std::size_t NR=64;
  const std::size_t NC=64;
  typedef sl::dense_array<float, 2, void> float_array2_t;

  float_array2_t in_data(NR,NC);
  for (std::size_t i=0; i<NR; ++i) {
    for (std::size_t j=0; j<NC; ++j) {
      in_data(i,j) = std::sin(float(i))*std::cos(float(j));
    }
  }
  float_array2_t xformed_data(NR,NC);
  float_array2_t in_data2(NR,NC);

  xform->forward2d(in_data, xformed_data, in_data.extent()[0], in_data.extent()[1]);
  xform->backward2d(xformed_data, in_data2, in_data.extent()[0], in_data.extent()[1]);

  tester.test("power of 2: backward2d(forward2d(signal)) == signal",
              sl::rms(in_data, in_data2), sl::intervalf("0.0000"));
  if (xform->is_orthonormal()) {
    tester.test("energy(spatial domain) == energy(wavelet domain)",
                float(sl::abs(sl::energy(in_data) - sl::energy(xformed_data))),sl::intervalf( "0.00"));
  }
  
  failed_test_count += tester.failed_test_count();
}

void test_embedded_zerotree_codec() {
  sl::tester tester("embedded_zerotree_codec<int16_t>");
  
  typedef sl::embedded_zerotree_codec<sl::int16_t> embedded_zerotree_codec_t;
  typedef embedded_zerotree_codec_t::int_t int_t;
  typedef embedded_zerotree_codec_t::int_matrix_t int_matrix_t;

  embedded_zerotree_codec_t embedded_zerotree_codec;
  int_matrix_t shapiro_test_case(8,8);
  shapiro_test_case =
     63,-34, 49, 10,  7, 13,-12,  7,
    -31, 23, 14,-13,  3,  4,  6, -1,
     15, 14,  3,-12,  5, -7,  3,  9,
    -9,  -7,-14,  8,  4, -2,  3,  2,
    -5,   9, -1, 47,  4,  6, -2,  2,
     3,   0, -3,  2,  3, -2,  0,  4,
     2,  -3,  6, -4,  3,  6,  3,  6,
     5,  11,  5,  6,  0,  3, -4,  4;
  
  std::vector<sl::uint8_t> buf(1024+shapiro_test_case.count()*sizeof(int_t));
  std::size_t buf_size = buf.size();
  std::size_t target_size = shapiro_test_case.count()*sizeof(int_t);
  std::size_t actual_size;
  double actual_rms;
  embedded_zerotree_codec.compress_to_target_size(shapiro_test_case,
                                                  target_size,
                                                  &buf[0],
                                                  buf_size,
                                                  &actual_size,
                                                  &actual_rms);
  int_matrix_t decoded_matrix;
  embedded_zerotree_codec.decompress(decoded_matrix, &buf[0], actual_size);

  tester.test("compressed size", actual_size <= target_size);
  tester.test("Lossless compression", shapiro_test_case == decoded_matrix, true);
  
  failed_test_count += tester.failed_test_count();

  if (tester.failed_test_count()) {
    std::cerr << "IN: " << std::endl;
    dump_matrix(shapiro_test_case);
    std::cerr << "OUT: " << std::endl;
    dump_matrix(decoded_matrix);
  }
}

template <class T>
void test_ezw_array_codec(const std::string& tpname,
                          const sl::dense_array<T,2,void>& array) {
  sl::tester tester("ezw_array_codec<" + tpname + "> for size " + sl::to_string(array.extent()[0]) + "x" + sl::to_string(array.extent()[1]));
  sl::ezw_array_codec codec;

  const std::size_t float_array_size = array.count()*sizeof(float);
  const std::size_t max_buf_size = 1024+float_array_size;
  char* buf = new char[max_buf_size];
  
  std::vector< std::pair<sl::ezw_array_codec::transform_kind_t , std::string> > xform_kinds;
  xform_kinds.push_back( std::make_pair( sl::ezw_array_codec::TRANSFORM_KIND_HAAR,                "Haar") );
  xform_kinds.push_back( std::make_pair( sl::ezw_array_codec::TRANSFORM_KIND_DAUBECHIES4,         "Daub4") );
  xform_kinds.push_back( std::make_pair( sl::ezw_array_codec::TRANSFORM_KIND_SYMMLET8,            "Symmlet8") );
  xform_kinds.push_back( std::make_pair( sl::ezw_array_codec::TRANSFORM_KIND_AUTO, "Auto") );
  
  for (std::size_t target_size=2; target_size<=float_array_size; target_size*=2) {
    for (std::size_t tkind = 0; tkind<xform_kinds.size(); ++tkind) {
      codec.set_transform_kind(xform_kinds[tkind].first);
      
      std::size_t buf_size = max_buf_size;
      std::size_t actual_size;
      double actual_rms;
      codec.compress_to_target_size(array, target_size, buf, buf_size, &actual_size, &actual_rms);

      sl::dense_array<T,2,void> array2;
      codec.decompress(array2, buf, buf_size);

      tester.test(xform_kinds[tkind].second + " - Compress to size " + sl::to_string(target_size) + ": Good extent", array.extent(), array2.extent());
      tester.test(xform_kinds[tkind].second + " - Compress to size " + sl::to_string(target_size) + ": Good size", actual_size <= target_size+6 || target_size < 24, true);
#if 1
      std::cerr << xform_kinds[tkind].second + " - Desired size: " << target_size << " - Actual size: " << actual_size << " RMS = " << sl::rms(array,array2) << " PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
    }
  }

  for (std::size_t tkind = 0; tkind<xform_kinds.size(); ++tkind) {
    codec.set_transform_kind(xform_kinds[tkind].first);
    std::size_t buf_size = max_buf_size;
    float target_rms = 1.0f;
    std::size_t actual_size;
    double actual_rms;
    codec.compress_to_target_rms_error(array, target_rms, buf, buf_size, &actual_size, &actual_rms);
    sl::dense_array<T,2,void> array2;
    codec.decompress(array2, buf, buf_size);
    tester.test(xform_kinds[tkind].second + " - Compress to rms " + sl::to_string(target_rms) + ": Good extent", array.extent(), array2.extent());
    tester.test(xform_kinds[tkind].second + " - Compress to rms " + sl::to_string(target_rms) + ": Good rms", actual_rms <= target_rms, true);
    if (sl::bitops<int>::is_power2(array.extent()[0]) && sl::bitops<int>::is_power2(array.extent()[1])) tester.test(xform_kinds[tkind].second + " - Compress to rms " + sl::to_string(target_rms) + ": Good rms approx", float(sl::abs(actual_rms - sl::rms(array,array2))), sl::intervalf(0.0f, 0.1f));
#if 1
    std::cerr << xform_kinds[tkind].second + " - Desired rms: " << target_rms << " - Actual size: " << actual_size << " RMS = " << sl::rms(array,array2) << " (real) vs. " << actual_rms << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
  }

  delete [] buf;

  failed_test_count += tester.failed_test_count();
}

void test_ezw_array_codec() {
  test_ezw_array_codec("uint8_t", make_array<sl::uint8_t>(64,64));
  test_ezw_array_codec("int8_t", make_array<sl::int8_t>(64,64));
  test_ezw_array_codec("uint16_t", make_array<sl::uint16_t>(64,64));
  test_ezw_array_codec("int16_t", make_array<sl::int16_t>(64,64));
  test_ezw_array_codec("float", make_array<float>(64,64));
  test_ezw_array_codec("float", make_array<float>(7,7));
  test_ezw_array_codec("float", make_array<float>(9,9));
}

template <class T>
void test_geometric_bandelet_array_codec(const std::string& tpname,
                          const sl::dense_array<T,2,void>& array) {
  sl::tester tester("geometric_bandelet_array_codec<" + tpname + "> for size " + sl::to_string(array.extent()[0]) + "x" + sl::to_string(array.extent()[1]));
  sl::geometric_bandelet_array_codec codec;

  const std::size_t float_array_size = array.count()*sizeof(float);
  const std::size_t max_buf_size = 1024+float_array_size;
  char* buf = new char[max_buf_size];
  
  std::vector< std::pair<sl::geometric_bandelet_array_codec::transform_kind_t , std::string> > xform_kinds;
  xform_kinds.push_back( std::make_pair( sl::geometric_bandelet_array_codec::TRANSFORM_KIND_HAAR,                "Haar") );
  xform_kinds.push_back( std::make_pair( sl::geometric_bandelet_array_codec::TRANSFORM_KIND_DAUBECHIES4,         "Daub4") );
  xform_kinds.push_back( std::make_pair( sl::geometric_bandelet_array_codec::TRANSFORM_KIND_SYMMLET8,            "Symmlet8") );
  xform_kinds.push_back( std::make_pair( sl::geometric_bandelet_array_codec::TRANSFORM_KIND_AUTO, "Auto") );
  
  for (std::size_t target_size=2; target_size<=float_array_size; target_size*=2) {
    for (std::size_t tkind = 0; tkind<xform_kinds.size(); ++tkind) {
      codec.set_transform_kind(xform_kinds[tkind].first);
      
      std::size_t buf_size = max_buf_size;
      std::size_t actual_size;
      double actual_rms;
      codec.compress_to_target_size(array, target_size, buf, buf_size, &actual_size, &actual_rms);

      sl::dense_array<T,2,void> array2;
      codec.decompress(array2, buf, buf_size);
      
      tester.test(xform_kinds[tkind].second + " - Compress to size " + sl::to_string(target_size) + ": Good extent", array.extent(), array2.extent());
      tester.test(xform_kinds[tkind].second + " - Compress to size " + sl::to_string(target_size) + ": Good size", actual_size <= target_size+6 || target_size < 24, true);
#if 1
      std::cerr << xform_kinds[tkind].second + " - Desired size: " << target_size << " - Actual size: " << actual_size << " RMS = " << sl::rms(array,array2) << " PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
    }
  }

  for (std::size_t tkind = 0; tkind<xform_kinds.size(); ++tkind) {
    codec.set_transform_kind(xform_kinds[tkind].first);
    std::size_t buf_size = max_buf_size;
    float target_rms = 1.0f;
    std::size_t actual_size;
    double actual_rms;
    codec.compress_to_target_rms_error(array, target_rms, buf, buf_size, &actual_size, &actual_rms);
    sl::dense_array<T,2,void> array2;
    codec.decompress(array2, buf, buf_size);
    tester.test(xform_kinds[tkind].second + " - Compress to rms " + sl::to_string(target_rms) + ": Good extent", array.extent(), array2.extent());
    tester.test(xform_kinds[tkind].second + " - Compress to rms " + sl::to_string(target_rms) + ": Good rms", actual_rms <= target_rms, true);
    if (sl::bitops<int>::is_power2(array.extent()[0]) && sl::bitops<int>::is_power2(array.extent()[1])) tester.test(xform_kinds[tkind].second + " - Compress to rms " + sl::to_string(target_rms) + ": Good rms approx", float(sl::abs(actual_rms - sl::rms(array,array2))), sl::intervalf(0.0f, 0.1f));
#if 1
    std::cerr << xform_kinds[tkind].second + " - Desired rms: " << target_rms << " - Actual size: " << actual_size << " RMS = " << sl::rms(array,array2) << " (real) vs. " << actual_rms << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
  }

  delete [] buf;

  failed_test_count += tester.failed_test_count();
}

void test_geometric_bandelet_array_codec() {
  test_geometric_bandelet_array_codec("uint8_t", make_array<sl::uint8_t>(64,64));
  test_geometric_bandelet_array_codec("int8_t", make_array<sl::int8_t>(64,64));
  test_geometric_bandelet_array_codec("uint16_t", make_array<sl::uint16_t>(64,64));
  test_geometric_bandelet_array_codec("int16_t", make_array<sl::int16_t>(64,64));
  test_geometric_bandelet_array_codec("float", make_array<float>(64,64));
  test_geometric_bandelet_array_codec("float", make_array<float>(7,7));
  test_geometric_bandelet_array_codec("float", make_array<float>(9,9));
}

template <class T>
void test_micro_jpgls_array_codec(const std::string& tpname,
                          const sl::dense_array<T,2,void>& array) {
  sl::tester tester("micro_jpgls_array_codec<" + tpname + "> for size " + sl::to_string(array.extent()[0]) + "x" + sl::to_string(array.extent()[1]));
  sl::micro_jpgls_array_codec codec;

  const std::size_t float_array_size = array.count()*sizeof(float);
  const std::size_t max_buf_size = 1024+2*float_array_size;
  char* buf = new char[max_buf_size];
        
  std::size_t buf_size = max_buf_size;
  std::size_t actual_size;
  double actual_rms;
  codec.compress(array, buf, buf_size, &actual_size, &actual_rms);
  
  sl::dense_array<T,2,void> array2;
  codec.decompress(array2, buf, buf_size);

  tester.test("Good extent", array.extent(), array2.extent());
  tester.test("Good rms", float(actual_rms), sl::intervalf(0.0f, 0.1f));
  tester.test("Good rms approx", float(sl::abs(actual_rms - sl::rms(array,array2))), sl::intervalf(0.0f, 0.005f));
#if 1
  std::cerr << "Compressed size: " << actual_size << " bpp "  << actual_size*8/float(array.count()) << " RMS = " << sl::rms(array,array2) << " (real) vs. " << actual_rms << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
      
  delete [] buf;

  failed_test_count += tester.failed_test_count();
}

void test_micro_jpgls_array_codec() {
  test_micro_jpgls_array_codec("uint8_t", make_array<sl::uint8_t>(64,64));
  test_micro_jpgls_array_codec("int8_t", make_array<sl::int8_t>(64,64));
  test_micro_jpgls_array_codec("uint16_t", make_array<sl::uint16_t>(64,64));
  test_micro_jpgls_array_codec("int16_t", make_array<sl::int16_t>(64,64));
  test_micro_jpgls_array_codec("float", make_array<float>(64,64));
  test_micro_jpgls_array_codec("float", make_array<float>(7,7));
  test_micro_jpgls_array_codec("float", make_array<float>(9,9));

  test_micro_jpgls_array_codec("zero uint8_t", make_zero_array<sl::uint8_t>(64,64));
  test_micro_jpgls_array_codec("zero int8_t", make_zero_array<sl::int8_t>(64,64));
  test_micro_jpgls_array_codec("zero uint16_t", make_zero_array<sl::uint16_t>(64,64));
  test_micro_jpgls_array_codec("zero int16_t", make_zero_array<sl::int16_t>(64,64));
  test_micro_jpgls_array_codec("zero float", make_zero_array<float>(64,64));
  test_micro_jpgls_array_codec("zero float", make_zero_array<float>(7,7));
  test_micro_jpgls_array_codec("zero float", make_zero_array<float>(9,9));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <class T>
void test_hidwt_array_codec(const std::string& tpname,
                          const sl::dense_array<T,2,void>& array) {
  sl::tester tester("hidwt_array_codec<" + tpname + "> for size " + sl::to_string(array.extent()[0]) + "x" + sl::to_string(array.extent()[1]));
  sl::hidwt_array_codec codec;

  const std::size_t float_array_size = array.count()*sizeof(float);
  const std::size_t max_buf_size = 1024+2*float_array_size;
  char* buf = new char[max_buf_size];

  for (float target_rms = 0.0f; target_rms<sl::amax(array); target_rms = std::max(2.0f*target_rms, 0.1f)) {
    std::size_t buf_size = max_buf_size;
    std::size_t actual_size;
    double actual_rms;
    codec.compress_to_target_rms_error(array,
                                       target_rms,
                                       buf,
                                       buf_size, &actual_size, &actual_rms);
    
    sl::dense_array<T,2,void> array2;
    codec.decompress(array2, buf, buf_size);
    
    tester.test("Good extent", array.extent(), array2.extent());
    if (target_rms==0.0f) {
      tester.test("Good rms", float(actual_rms), sl::intervalf(0.0f, 0.05f));
    } 

    if (sl::bitops<int>::is_power2(array.extent()[0]) && sl::bitops<int>::is_power2(array.extent()[1])) {
      // FIXME 
      std::size_t old_failed_count = tester.failed_test_count();

      tester.test("Good rms approx", float(sl::abs(actual_rms - sl::rms(array,array2))), sl::intervalf(0.0f, 0.3f*std::max(0.2f,sl::rms(array,array2))));

      if ( tester.failed_test_count() != old_failed_count ) {
	std::cerr << "ORIGINAL" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DECOMPRESSED" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array2(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DELTA" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << (array(i,j)-array2(i,j));
	  }
	  std::cerr << std::endl;
	}
      }
    }
#if 1
      std::cerr << "Target RMS: " << target_rms << " Compressed size: " << actual_size << " bpp "  << actual_size*8/float(array.count()) << " RMS = " << sl::rms(array,array2) << " (real) vs. " << actual_rms << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
  }
  
  for (float target_amax = 0.0f; target_amax<sl::amax(array); target_amax = std::max(2.0f*target_amax, 0.1f)) {
    std::size_t buf_size = max_buf_size;
    std::size_t actual_size;
    double actual_amax;
    codec.compress_to_target_amax_error(array,
                                       target_amax,
                                       buf,
                                       buf_size, &actual_size, &actual_amax);
    
    sl::dense_array<T,2,void> array2;
    codec.decompress(array2, buf, buf_size);
    
    tester.test("Good extent", array.extent(), array2.extent());
    if (target_amax==0.0f) {
      tester.test("Good amax", float(actual_amax), sl::intervalf(0.0f, 0.05f));
    } 

    if (sl::bitops<int>::is_power2(array.extent()[0]) && sl::bitops<int>::is_power2(array.extent()[1])) {
      // FIXME 
      std::size_t old_failed_count = tester.failed_test_count();

      tester.test("Good amax approx", float(sl::abs(actual_amax - sl::amax_diff(array,array2))), sl::intervalf(0.0f, 0.3f*std::max(0.2f,sl::amax_diff(array,array2))));

      if ( tester.failed_test_count() != old_failed_count ) {
	std::cerr << "ORIGINAL" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DECOMPRESSED" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array2(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DELTA" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << (array(i,j)-array2(i,j));
	  }
	  std::cerr << std::endl;
	}
      }
    }
#if 1
      std::cerr << "Target AMAX: " << target_amax << " Compressed size: " << actual_size << " bpp "  << actual_size*8/float(array.count()) << " AMAX = " << sl::amax_diff(array,array2) << " (real) vs. " << actual_amax << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
  }

  delete [] buf;

  failed_test_count += tester.failed_test_count();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void test_hidwt_array_codec() {
  test_hidwt_array_codec("uint8_t", make_array<sl::uint8_t>(64,64));
  test_hidwt_array_codec("int8_t", make_array<sl::int8_t>(64,64));
  test_hidwt_array_codec("uint16_t", make_array<sl::uint16_t>(64,64));
  test_hidwt_array_codec("int16_t", make_array<sl::int16_t>(64,64));
  test_hidwt_array_codec("float", make_array<float>(64,64));
  test_hidwt_array_codec("float", make_array<float>(7,7));
  test_hidwt_array_codec("float", make_array<float>(9,9));

  test_hidwt_array_codec("zero uint8_t", make_zero_array<sl::uint8_t>(64,64));
  test_hidwt_array_codec("zero int8_t", make_zero_array<sl::int8_t>(64,64));
  test_hidwt_array_codec("zero uint16_t", make_zero_array<sl::uint16_t>(64,64));
  test_hidwt_array_codec("zero int16_t", make_zero_array<sl::int16_t>(64,64));
  test_hidwt_array_codec("zero float", make_zero_array<float>(64,64));
  test_hidwt_array_codec("zero float", make_zero_array<float>(7,7));
  test_hidwt_array_codec("zero float", make_zero_array<float>(9,9));
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <class T>
void test_quantized_array_codec(const std::string& tpname,
				const sl::dense_array<T,2,void>& array) {
  sl::tester tester("quantized_array_codec<" + tpname + "> for size " + sl::to_string(array.extent()[0]) + "x" + sl::to_string(array.extent()[1]));
  sl::quantized_array_codec codec;
  codec.set_is_compressing_header(true);

  const std::size_t float_array_size = array.count()*sizeof(float);
  const std::size_t max_buf_size = 1024+2*float_array_size;
  char* buf = new char[max_buf_size];

  for (float target_rms = 0.0f; target_rms<sl::amax(array); target_rms = std::max(2.0f*target_rms, 0.1f)) {
    std::size_t buf_size = max_buf_size;
    std::size_t actual_size;
    double actual_rms;
    codec.compress_to_target_rms_error(array,
                                       target_rms,
                                       buf,
                                       buf_size, &actual_size, &actual_rms);
    
    sl::dense_array<T,2,void> array2;
    codec.decompress(array2, buf, buf_size);
    
    tester.test("Good extent", array.extent(), array2.extent());
    if (target_rms==0.0f) {
      tester.test("Good rms", float(actual_rms), sl::intervalf(0.0f, 0.05f));
    } 

    if (sl::bitops<int>::is_power2(array.extent()[0]) && sl::bitops<int>::is_power2(array.extent()[1])) {
      // FIXME 
      std::size_t old_failed_count = tester.failed_test_count();

      tester.test("Good rms approx", float(sl::abs(actual_rms - sl::rms(array,array2))), sl::intervalf(0.0f, 0.3f*std::max(0.2f,sl::rms(array,array2))));

      if ( tester.failed_test_count() != old_failed_count ) {
	std::cerr << "ORIGINAL" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DECOMPRESSED" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array2(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DELTA" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << (array(i,j)-array2(i,j));
	  }
	  std::cerr << std::endl;
	}
      }
    }
#if 1
      std::cerr << "Target RMS: " << target_rms << " Compressed size: " << actual_size << " bpp "  << actual_size*8/float(array.count()) << " RMS = " << sl::rms(array,array2) << " (real) vs. " << actual_rms << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
  }
  
  for (float target_amax = 0.0f; target_amax<sl::amax(array); target_amax = std::max(2.0f*target_amax, 0.1f)) {
    std::size_t buf_size = max_buf_size;
    std::size_t actual_size;
    double actual_amax;
    codec.compress_to_target_amax_error(array,
                                       target_amax,
                                       buf,
                                       buf_size, &actual_size, &actual_amax);
    
    sl::dense_array<T,2,void> array2;
    codec.decompress(array2, buf, buf_size);
    
    tester.test("Good extent", array.extent(), array2.extent());
    if (target_amax==0.0f) {
      tester.test("Good amax", float(actual_amax), sl::intervalf(0.0f, 0.05f));
    } 

    if (sl::bitops<int>::is_power2(array.extent()[0]) && sl::bitops<int>::is_power2(array.extent()[1])) {
      // FIXME 
      std::size_t old_failed_count = tester.failed_test_count();

      tester.test("Good amax approx", float(sl::abs(actual_amax - sl::amax_diff(array,array2))), sl::intervalf(0.0f, 0.3f*std::max(0.2f,sl::amax_diff(array,array2))));

      if ( tester.failed_test_count() != old_failed_count ) {
	std::cerr << "ORIGINAL" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DECOMPRESSED" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << array2(i,j);
	  }
	  std::cerr << std::endl;
	}
	std::cerr << "DELTA" << std::endl;
	for (std::size_t i=0; i< array.extent()[0]; ++i) {
	  for (std::size_t j=0; j<array.extent()[1]; ++j) {
	    std::cerr << std::setw(8) << (array(i,j)-array2(i,j));
	  }
	  std::cerr << std::endl;
	}
      }
    }
#if 1
      std::cerr << "Target AMAX: " << target_amax << " Compressed size: " << actual_size << " bpp "  << actual_size*8/float(array.count()) << " AMAX = " << sl::amax_diff(array,array2) << " (real) vs. " << actual_amax << " (approx) PSNR = " << sl::psnr(array, array2) << " dB" << std::endl;
#endif
  }

  delete [] buf;

  failed_test_count += tester.failed_test_count();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void test_quantized_array_codec() {
  test_quantized_array_codec("uint8_t", make_array<sl::uint8_t>(64,64));
  test_quantized_array_codec("int8_t", make_array<sl::int8_t>(64,64));
  test_quantized_array_codec("uint16_t", make_array<sl::uint16_t>(64,64));
  test_quantized_array_codec("int16_t", make_array<sl::int16_t>(64,64));
  test_quantized_array_codec("float", make_array<float>(64,64));
  test_quantized_array_codec("float", make_array<float>(7,7));
  test_quantized_array_codec("float", make_array<float>(9,9));

  test_quantized_array_codec("zero uint8_t", make_zero_array<sl::uint8_t>(64,64));
  test_quantized_array_codec("zero int8_t", make_zero_array<sl::int8_t>(64,64));
  test_quantized_array_codec("zero uint16_t", make_zero_array<sl::uint16_t>(64,64));
  test_quantized_array_codec("zero int16_t", make_zero_array<sl::int16_t>(64,64));
  test_quantized_array_codec("zero float", make_zero_array<float>(64,64));
  test_quantized_array_codec("zero float", make_zero_array<float>(7,7));
  test_quantized_array_codec("zero float", make_zero_array<float>(9,9));
}
// -----------------------------------------------------------------------------



int main() {
#if 0
  test_hidwt_array_codec("float", make_array<float>(64,64));
#else
  test_arithmetic_codec(2);
  test_arithmetic_codec(8);
  test_arithmetic_codec(16);
  test_range_codec(2);
  test_range_codec(8);
  test_range_codec(16);

  test_wavelet1d("haar_wavelet_transform", new sl::haar_wavelet_transform<float>());
  test_wavelet1d("normalized_haar_wavelet_transform", new sl::normalized_haar_wavelet_transform<float>());
  test_wavelet1d("ts_wavelet_transform", new sl::ts_wavelet_transform<float>());
  test_wavelet1d("normalized_daubechies4_wavelet_transform", new sl::normalized_daubechies4_wavelet_transform<float>());
  test_wavelet1d("normalized_symmlet8_wavelet_transform", new sl::normalized_symmlet8_wavelet_transform<float>());
  test_wavelet2d("haar_wavelet_transform", new sl::haar_wavelet_transform<float>());
  test_wavelet2d("normalized_haar_wavelet_transform", new sl::normalized_haar_wavelet_transform<float>());
  test_wavelet2d("ts_wavelet_transform", new sl::ts_wavelet_transform<float>());
  test_wavelet2d("normalized_daubechies4_wavelet_transform", new sl::normalized_daubechies4_wavelet_transform<float>());
  test_wavelet2d("normalized_symmlet8_wavelet_transform", new sl::normalized_symmlet8_wavelet_transform<float>());
  test_embedded_zerotree_codec();

  test_micro_jpgls_array_codec();

  test_hidwt_array_codec();

  test_ezw_array_codec();

  test_quantized_array_codec();

  //test_geometric_bandelet_array_codec();
#endif  
  return (int)failed_test_count;
}
