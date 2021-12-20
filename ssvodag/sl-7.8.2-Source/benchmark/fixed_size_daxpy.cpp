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
// sl::fixed_size_vector<sl::column_orientation,N,T> DAXPY benchmark

#include <sl/fixed_size_vector.hpp>
#include <sl/clock.hpp>
#include <iomanip>

static unsigned int SEED = 93186752;

static float random (float l, float h)  {
  static unsigned int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901;

  SEED = a*(SEED % q) - r*(SEED / q);
  float result = l + ((float)SEED / (float)m) * (h-l);
  
  return result;
}

template<class T>
void optimizationSink(T&);

template<size_t N_rank, class value_t>
void fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,N_rank,value_t>, 
		       size_t iters, value_t a) {
    sl::cpu_time_clock timer;

    sl::fixed_size_vector<sl::column_orientation,N_rank,value_t> ta, tb, tc, td, te, tf, tg, th, ti, tj;
    for (size_t i=0; i < N_rank; ++i) {
      ta[i] = value_t(random(1.0f,2.0f));
      tb[i] = value_t(random(1.0f,2.0f));
      tc[i] = value_t(random(1.0f,2.0f));
      td[i] = value_t(random(1.0f,2.0f));
      te[i] = value_t(random(1.0f,2.0f));
      tf[i] = value_t(random(1.0f,2.0f));
      tg[i] = value_t(random(1.0f,2.0f));
      th[i] = value_t(random(1.0f,2.0f));
      ti[i] = value_t(random(1.0f,2.0f));
      tj[i] = value_t(random(1.0f,2.0f));
    }

    value_t b = -a;

    double numFlops = 0;
    double elapsedSeconds = 0;

    if (N_rank < 20) {
      timer.restart();
      for (size_t i=0; i < iters; ++i) {
        ta += a * tb;
        tc += a * td;
        te += a * tf;
        tg += a * th;
        ti += a * tj;
        tb += b * ta;
        td += b * tc;
        tf += b * te;
        th += b * tg;
        tj += b * ti;
        ta += a * tb;
        tc += a * td;
        te += a * tf;
        tg += a * th;
        ti += a * tj;
        tb += b * ta;
        td += b * tc;
        tf += b * te;
        th += b * tg;
        tj += b * ti;
      }
      elapsedSeconds = timer.elapsed().as_microseconds() / 1E6;
      numFlops = 40.0 * N_rank * double(iters);
    } else {
      timer.restart();
      for (size_t i=0; i < iters; ++i) {
        ta += a * tb;
        tb += b * ta;
      }
      elapsedSeconds = timer.elapsed().as_microseconds() / 1E6;
      numFlops = 4.0 * N_rank * double(iters);
    }
    
    optimizationSink(ta);
    optimizationSink(tb);
    optimizationSink(tc);
    optimizationSink(td);
    optimizationSink(te);
    optimizationSink(tf);
    optimizationSink(tg);
    optimizationSink(th);
    optimizationSink(ti);
    optimizationSink(tj);
    
    if (elapsedSeconds > 0) {
      double Mflops = numFlops / (1.0e+6) / elapsedSeconds;
      std::cout << std::setw(5) << N_rank << '\t' << Mflops << std::endl;
    }
}    
 
double a = 0.3429843;
 
template<class T>
void optimizationSink(T&)
{
}
 
int main() {
    std::cout << "sl::fixed_size_vector<sl::column_orientation,N,float> DAXPY benchmark" << std::endl
         << std::setw(5) << "N" << '\t' << "Mflops/s" << std::endl;
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,1,float>(), 800000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,2,float>(), 800000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,3,float>(), 800000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,4,float>(), 700000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,5,float>(), 600000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,6,float>(), 500000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,7,float>(), 500000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,8,float>(), 500000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,9,float>(), 500000, float(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,10,float>(), 500000, float(a));

    std::cout << "sl::fixed_size_vector<sl::column_orientation,N,double> DAXPY benchmark" << std::endl
         << std::setw(5) << "N" << '\t' << "Mflops/s" << std::endl;
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,1,double>(), 800000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,2,double>(), 800000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,3,double>(), 800000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,4,double>(), 700000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,5,double>(), 600000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,6,double>(), 500000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,7,double>(), 500000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,8,double>(), 500000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,9,double>(), 500000, double(a));
    fsvDAXPYBenchmark(sl::fixed_size_vector<sl::column_orientation,10,double>(), 500000, double(a));
    
    return 0;
}   




