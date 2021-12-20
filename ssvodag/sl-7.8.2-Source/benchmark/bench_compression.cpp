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

#include "sl/quantized_array_codec.hpp"
#include "sl/geometric_bandelet_array_codec.hpp"
#include "sl/micro_jpgls_array_codec.hpp"
#include "sl/hidwt_array_codec.hpp"
#include "sl/ezw_array_codec.hpp"
#include <sl/clock.hpp>
#include <stdio.h>

// Perlin noise stuff

class perlin {
public:

  /// A random value between -1 and +1
  static float noise(sl::int32_t x) {
    x = (x<<13) ^ x;
    return ( 1.0f - ((float)( (x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff)) / 1073741824.0f);    
  }

  /// A random value between -1 and +1
  static float noise(sl::int32_t x, sl::int32_t y) {
    return noise(x+y+57);
  }

  static float smoothed_noise(sl::int32_t x) {
    return 0.5f*noise(x) + 0.25f*noise(x-1) + 0.25f*noise(x+1);
  }

  static float smoothed_noise(sl::int32_t x, sl::int32_t y) {
    float corners = 0.0625f * (noise(x-1,y-1)+noise(x+1,y-1)+noise(x-1,y+1)+noise(x+1,y+1));
    float sides   = 0.1250f * (noise(x-1,y  )+noise(x+1,y  )+noise(x,y-1  )+noise(x,y+1  ));
    float center  = 0.2500f * noise(x,y);
    return corners+sides+center;
  }
  
  static float smoothed_noise(float x) {
    sl::int32_t ix = sl::int32_t(x);
    float       fx = x-float(ix);
    float v1 = smoothed_noise(ix);
    float v2 = smoothed_noise(ix+1);
    float f = (1.0f - std::cos(fx * 3.14159265f)) * 0.5f;
    return  v1*(1-f) + v2*f;
  }

  static float smoothed_noise(float x, float y) {
    sl::int32_t ix = sl::int32_t(x);
    float       fx = x-float(ix);
    sl::int32_t iy = sl::int32_t(y);
    float       fy = y-float(iy);
    
    float v1 = smoothed_noise(ix,   iy);
    float v2 = smoothed_noise(ix+1, iy);
    float v3 = smoothed_noise(ix,   iy+1);
    float v4 = smoothed_noise(ix+1, iy+1);
    
    float f1 = (1.0f - std::cos(fx * 3.14159265f)) * 0.5f;
    float f2 = (1.0f - std::cos(fy * 3.14159265f)) * 0.5f;

    float i1 = v1*(1-f1)+v2*f1;
    float i2 = v3*(1-f1)+v4*f1;

    return i1*(1-f2)+i2*f2;
  }
  
  static float noise(float x,
                     float persistence,
                     std::size_t octave_count) {
    float result = 0;
    float p = persistence;
    float frequency = 1;
    float amplitude = 1;
    for (std::size_t i=0; i<octave_count; ++i) {
      result += amplitude * smoothed_noise(x*frequency);
      frequency *= 2.0f;
      amplitude *= p;
    }
    return result;
  }

  static float noise(float x, float y,
                     float persistence,
                     std::size_t octave_count) {
    float result = 0;
    float p = persistence;
    float frequency = 1;
    float amplitude = 1;
    for (std::size_t i=0; i<octave_count; ++i) {
      result += amplitude * smoothed_noise(x*frequency, y*frequency);
      frequency *= 2.0f;
      amplitude *= p;
    }
    return result;
  }
};

template <class T>
sl::dense_array<T,2,void> make_array(std::size_t height, std::size_t width) {
  typedef T value_t;
  typedef sl::dense_array<T,2,void> array_t;
  
  array_t result(height, width);
  for (std::size_t i=0; i<height; ++i) {
    for (std::size_t j=0; j<width; ++j) {
#if 1
      // Emulate a c-bdam delta distribution
      
      float x00    = perlin::noise(float(i+0)/float(height),float(j+0)/float(width), 0.25f, std::size_t(4));
      float x01    = perlin::noise(float(i+0)/float(height),float(j+1)/float(width), 0.25f, std::size_t(4));
      float x10    = perlin::noise(float(i+1)/float(height),float(j+0)/float(width), 0.25f, std::size_t(4));
      float x11    = perlin::noise(float(i+1)/float(height),float(j+1)/float(width), 0.25f, std::size_t(4));
      float xpred  = 0.25f * (x00+x01+x10+x11);
      float xhires = perlin::noise(float(i+0.5f)/float(height),float(j+0.5f)/float(width), 0.25f, std::size_t(4));
      float dx     = xhires-xpred;
      dx *= 32767.0f;
      //std::cerr << "noise(" << i << " " << j << ") = " << dx << std::endl;
      result(i,j) = value_t(dx);
#else
      result(i,j) = (i*13+j*17)%16384;
#endif
    }
  }
  return result;
}

// ======================= RMS Benchmarks

template <class T>
void output_array_codec_performance_table_header_rms(const std::string& tpname,
                                                 const sl::dense_array<T,2,void>& array,
                                                 float target_rms) {
  fprintf(stdout, "========================================================================\n");
  fprintf(stdout, "Codec performance stats for %dx%d array of %s, Desired RMSE=%.2f\n",
          (int)array.extent()[0], (int)array.extent()[1], tpname.c_str(), target_rms);
  fprintf(stdout, "========================================================================\n");
  fprintf(stdout, "%18s %8s %8s %8s %8s %8s\n",
          "CODEC/Base", "RMS", "SZ", "BPP", "enc/s", "dec/s");
  fprintf(stdout, "------------------------------------------------------------------------\n");
}

void output_array_codec_performance_table_footer_rms() {
  fprintf(stdout, "------------------------------------------------------------------------\n\n");
}

template <class T>
void output_array_codec_performance_table_line_rms(sl::array_codec* codec,
						   const sl::dense_array<T,2,void>& array,
						   float target_rms) {
  const std::size_t array_size = array.count()*sizeof(T);
  const std::size_t max_buf_size = 1024+2*array_size;
  char* buf = new char[max_buf_size];

  std::size_t buf_size = max_buf_size;
  std::size_t actual_size = 0;
  double      actual_rms = 0;
  double      compress_rate = 0;
  double      decompress_rate = 0;
  {
    std::size_t N=20;
    sl::cpu_time_clock ck;
    ck.restart();
    for (std::size_t i=0; i<N; ++i) {
      codec->compress_to_target_rms_error(array,
                                          target_rms,
                                          buf,
                                          buf_size,
                                          &actual_size,
                                          &actual_rms);
    }
    double s = ck.elapsed().as_microseconds()/1.0E6;
    compress_rate = double(N)/s;
  }
  {
    std::size_t N=200;
    sl::dense_array<T,2,void> array2;
    sl::cpu_time_clock ck;
    ck.restart();
    for (std::size_t i=0; i<N; ++i) {
      codec->decompress(array2, buf, actual_size);
    }
    double s = ck.elapsed().as_microseconds()/1.0E6;
    decompress_rate = double(N)/s;
  }

  fprintf(stdout, "%18s %8.3f %8d %8.3f %8.2f %8.2f\n",
          codec->description().c_str(),
          (float)actual_rms,
          (int)actual_size,
          (float)(actual_size*8.0f/array.count()),
          (float)compress_rate,
          (float)decompress_rate);

  delete[] buf;
}

template <class T>
void output_array_codec_performance_table_rms(const std::string& tpname,
					      const sl::dense_array<T,2,void>& array,
					      float target_rms) {
  output_array_codec_performance_table_header_rms(tpname, array, target_rms);
  
  {
    sl::quantized_array_codec codec;
    output_array_codec_performance_table_line_rms(&codec, array, target_rms);
  }
  
  {
    sl::micro_jpgls_array_codec codec;
    output_array_codec_performance_table_line_rms(&codec, array, target_rms);
  }

  {
    sl::hidwt_array_codec codec;
    output_array_codec_performance_table_line_rms(&codec, array, target_rms);
  }

  {
    sl::ezw_array_codec codec;
    for (sl::int32_t tpkind = sl::int32_t(codec.first_transform_kind())+1;
         tpkind <= sl::int32_t(codec.last_transform_kind());
         ++tpkind) {
      codec.set_transform_kind(sl::wavelet_array_codec::transform_kind_t(tpkind));
      output_array_codec_performance_table_line_rms(&codec, array, target_rms);
    }
  }

#if 0
  {
    sl::geometric_bandelet_array_codec codec;
    for (sl::int32_t tpkind = sl::int32_t(codec.first_transform_kind())+1;
         tpkind <= sl::int32_t(codec.last_transform_kind());
         ++tpkind) {
      codec.set_transform_kind(sl::wavelet_array_codec::transform_kind_t(tpkind));
      output_array_codec_performance_table_line_rms(&codec, array, target_rms);
    }
  }
#endif
  
  output_array_codec_performance_table_footer_rms();
}

// ======================= AMAX Benchmarks

template <class T>
void output_array_codec_performance_table_header_amax(const std::string& tpname,
                                                 const sl::dense_array<T,2,void>& array,
                                                 float target_amax) {
  fprintf(stdout, "========================================================================\n");
  fprintf(stdout, "Codec performance stats for %dx%d array of %s, Desired AMAXE=%.2f\n",
          (int)array.extent()[0], (int)array.extent()[1], tpname.c_str(), target_amax);
  fprintf(stdout, "========================================================================\n");
  fprintf(stdout, "%18s %8s %8s %8s %8s %8s\n",
          "CODEC/Base", "AMAX", "SZ", "BPP", "enc/s", "dec/s");
  fprintf(stdout, "------------------------------------------------------------------------\n");
}

void output_array_codec_performance_table_footer_amax() {
  fprintf(stdout, "------------------------------------------------------------------------\n\n");
}

template <class T>
void output_array_codec_performance_table_line_amax(sl::array_codec* codec,
						   const sl::dense_array<T,2,void>& array,
						   float target_amax) {
  const std::size_t array_size = array.count()*sizeof(T);
  const std::size_t max_buf_size = 1024+2*array_size;
  char* buf= new char[max_buf_size];

  std::size_t buf_size = max_buf_size;
  std::size_t actual_size = 0;
  double      actual_amax = 0;
  double      compress_rate = 0;
  double      decompress_rate = 0;
  {
    std::size_t N=20;
    sl::cpu_time_clock ck;
    ck.restart();
    for (std::size_t i=0; i<N; ++i) {
      codec->compress_to_target_amax_error(array,
                                          target_amax,
                                          buf,
                                          buf_size,
                                          &actual_size,
                                          &actual_amax);
    }
    double s = ck.elapsed().as_microseconds()/1.0E6;
    compress_rate = double(N)/s;
  }
  {
    std::size_t N=200;
    sl::dense_array<T,2,void> array2;
    sl::cpu_time_clock ck;
    ck.restart();
    for (std::size_t i=0; i<N; ++i) {
      codec->decompress(array2, buf, actual_size);
    }
    double s = ck.elapsed().as_microseconds()/1.0E6;
    decompress_rate = double(N)/s;
  }

  fprintf(stdout, "%18s %8.3f %8d %8.3f %8.2f %8.2f\n",
          codec->description().c_str(),
          (float)actual_amax,
          (int)actual_size,
          (float)(actual_size*8.0f/array.count()),
          (float)compress_rate,
          (float)decompress_rate);

  delete [] buf;
}

template <class T>
void output_array_codec_performance_table_amax(const std::string& tpname,
					      const sl::dense_array<T,2,void>& array,
					      float target_amax) {
  output_array_codec_performance_table_header_amax(tpname, array, target_amax);
  
  {
    sl::quantized_array_codec codec;
    output_array_codec_performance_table_line_amax(&codec, array, target_amax);
  }
  
  {
    sl::micro_jpgls_array_codec codec;
    output_array_codec_performance_table_line_amax(&codec, array, target_amax);
  }

  {
    sl::hidwt_array_codec codec;
    output_array_codec_performance_table_line_amax(&codec, array, target_amax);
  }

  {
    sl::ezw_array_codec codec;
    for (sl::int32_t tpkind = sl::int32_t(codec.first_transform_kind())+1;
         tpkind <= sl::int32_t(codec.last_transform_kind());
         ++tpkind) {
      codec.set_transform_kind(sl::wavelet_array_codec::transform_kind_t(tpkind));
      output_array_codec_performance_table_line_amax(&codec, array, target_amax);
    }
  }

#if 0
  {
    sl::geometric_bandelet_array_codec codec;
    for (int32_t sl::tpkind = sl::int32_t(codec.first_transform_kind())+1;
         tpkind <= sl::int32_t(codec.last_transform_kind());
         ++tpkind) {
      codec.set_transform_kind(sl::wavelet_array_codec::transform_kind_t(tpkind));
      output_array_codec_performance_table_line_amax(&codec, array, target_amax);
    }
  }
#endif
  
  output_array_codec_performance_table_footer_amax();
}

int main() {
  output_array_codec_performance_table_rms("float", make_array<float>(64,64), 0.00);
  output_array_codec_performance_table_rms("float", make_array<float>(64,64), 1.00);
  output_array_codec_performance_table_rms("float", make_array<float>(64,64),10.00);

  output_array_codec_performance_table_rms("int32_t", make_array<sl::int32_t>(64,64), 0.00);
  output_array_codec_performance_table_rms("int32_t", make_array<sl::int32_t>(64,64), 1.00);
  output_array_codec_performance_table_rms("int32_t", make_array<sl::int32_t>(64,64),10.00);
  
  output_array_codec_performance_table_rms("int16_t", make_array<sl::int16_t>(64,64), 0.00);
  output_array_codec_performance_table_rms("int16_t", make_array<sl::int16_t>(64,64), 1.00);
  output_array_codec_performance_table_rms("int16_t", make_array<sl::int16_t>(64,64),10.00);

  output_array_codec_performance_table_amax("float", make_array<float>(64,64), 0.00);
  output_array_codec_performance_table_amax("float", make_array<float>(64,64), 1.00);
  output_array_codec_performance_table_amax("float", make_array<float>(64,64),10.00);

  output_array_codec_performance_table_amax("int32_t", make_array<sl::int32_t>(64,64), 0.00);
  output_array_codec_performance_table_amax("int32_t", make_array<sl::int32_t>(64,64), 1.00);
  output_array_codec_performance_table_amax("int32_t", make_array<sl::int32_t>(64,64),10.00);
  
  output_array_codec_performance_table_amax("int16_t", make_array<sl::int16_t>(64,64), 0.00);
  output_array_codec_performance_table_amax("int16_t", make_array<sl::int16_t>(64,64), 1.00);
  output_array_codec_performance_table_amax("int16_t", make_array<sl::int16_t>(64,64),10.00);

  return 0;
}
