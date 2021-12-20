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
// sl::dense_vector<T,sl::xpr::vec_column_orientation_tag> DAXPY benchmark

#include <sl/dense_vector.hpp>
#include <sl/clock.hpp>
#include <iomanip>

static unsigned int SEED = 93186752;

static double random (double l, double h)  {
  static unsigned int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901;

  SEED = a*(SEED % q) - r*(SEED / q);
  double result = l + ((double)SEED / (double)m) * (h-l);
  
  return result;
}

template<class T>
void optimizationSink(T&);

template<class value_t>
void denseDAXPYBenchmark(sl::dense_vector<value_t,sl::xpr::vec_column_orientation_tag>, 
			 int N_rank,
			 size_t iters, value_t a) {
    sl::cpu_time_clock timer;

    sl::dense_vector<value_t, sl::xpr::vec_column_orientation_tag> 
      ta(N_rank), tb(N_rank), tc(N_rank), td(N_rank), te(N_rank), tf(N_rank), tg(N_rank), th(N_rank), ti(N_rank), tj(N_rank);
    for (size_t i=0; i < N_rank; ++i) {
      ta[i] = value_t(random(1.,2.));
      tb[i] = value_t(random(1.,2));
      tc[i] = value_t(random(1.,2.));
      td[i] = value_t(random(1.,2.));
      te[i] = value_t(random(1.,2.));
      tf[i] = value_t(random(1.,2.));
      tg[i] = value_t(random(1.,2.));
      th[i] = value_t(random(1.,2.));
      ti[i] = value_t(random(1.,2.));
      tj[i] = value_t(random(1.,2.));
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
    std::cout << "sl::dense_vector<sl::xpr::vec_column_orientation_tag,float>(N) DAXPY benchmark" << std::endl
         << std::setw(5) << "N" << '\t' << "Mflops/s" << std::endl;
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (1), 800000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (2), 800000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (3), 800000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (4), 700000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (5), 600000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (6), 500000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (7), 500000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (8), 500000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (9), 500000, float(a));
    denseDAXPYBenchmark(sl::dense_vector<float, sl::xpr::vec_column_orientation_tag>(), (10), 500000, float(a));

    std::cout << "sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(N) DAXPY benchmark" << std::endl
         << std::setw(5) << "N" << '\t' << "Mflops/s" << std::endl;
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (1), 800000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (2), 800000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (3), 800000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (4), 700000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (5), 600000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (6), 500000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (7), 500000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (8), 500000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (9), 500000, double(a));
    denseDAXPYBenchmark(sl::dense_vector<double, sl::xpr::vec_column_orientation_tag>(), (10), 500000, double(a));
    
    return 0;
}   




