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
#include <sl/fixed_size_point.hpp>
#include <sl/clock.hpp>
#include <sl/kdtree.hpp>
#include <sl/utility.hpp>
#include <vector>
#include <stdlib.h>

// -------------------------------------------

static inline void random_seed(unsigned int seed) {
  srand(seed);
}

static inline sl::point3f random_uniform_point() {
  float x = (rand() / (float)(RAND_MAX));
  float y = (rand() / (float)(RAND_MAX));
  float z = (rand() / (float)(RAND_MAX));
  return sl::point3f(x,y,z);
}


// -------------------------------------------

std::vector<sl::point3f> uniform_points(std::size_t sz) {
  std::vector<sl::point3f> result;

  random_seed(1234567);
  for (std::size_t i=0; i<sz; ++i) {
    sl::point3f p = random_uniform_point(); 
    result.push_back(p);
  }
  return result;
}

// -------------------------------------------

std::vector<sl::point3f> ordered_points(std::size_t sz) {
  std::vector<sl::point3f> result;

  sl::point3f o(0.,0.,0.);
  sl::point3f e(1.,1.,1.);

  random_seed(1234567);
  for (std::size_t i=0; i<sz; ++i) {
    float t = float(i) / float((sz>1) ? (sz-1) : 1);
    sl::point3f p;
    p = o + t * (e-o);
    result.push_back(p);
  }
  return result;
}

// -------------------------------------------

std::vector<sl::point3f> z_ordered_points(std::size_t sz) {
  std::vector<sl::point3f> result;

  random_seed(1234567);
  for (std::size_t i=0; i<sz; ++i) {
    float t = float(i) / float((sz>1) ? (sz-1) : 1);
    sl::point3f p = random_uniform_point();
    p[2] = t;
    result.push_back(p);
  }
  return result;
}

// -------------------------------------------

static void benchmark_kdt(const char* id,
			  sl::kdtree<3,float,sl::point3f>& kdt,
			  const std::vector<sl::point3f>& points) {
  typedef sl::kdtree<3,float,sl::point3f> kdt_t;

  sl::point3f pmin = points[0];
  sl::point3f pmax = points[0];
  for (std::size_t i=1; i<points.size(); ++i) {
    pmin[0] = std::min(points[i][0], pmin[0]);
    pmin[1] = std::min(points[i][1], pmin[1]);
    pmin[2] = std::min(points[i][2], pmin[2]);
    pmax[0] = std::max(points[i][0], pmax[0]);
    pmax[1] = std::max(points[i][1], pmax[1]);
    pmax[2] = std::max(points[i][2], pmax[2]);
  }
  
  kdt.clear();

  int kdt_create_ms = 0;
  {
    sl::cpu_time_clock clk;
    clk.restart();
    for (std::size_t i=0; i<points.size(); ++i) {
      kdt.insert(points[i]);
    }
    kdt_create_ms = (int)(clk.elapsed().as_microseconds()/1000);
  }

  std::cout << std::endl;
  std::cout << "== " << id << std::endl;
  {
    std::size_t kddeepest;
    std::size_t kdshortest;
    kdt.tree_depth_in(&kdshortest, &kddeepest);
    std::size_t kdcount = kdt.count();
    
    std::size_t kdinternal = kdt.internal_path_length();
    std::size_t kdexternal = kdt.external_path_length();
    
    double kdcost         = (kdinternal + kdexternal) / (kdcount == 0 ? 1.0 : (double)kdcount);
    double kdoptimum_cost = kdcount == 0 ? 0.0 : std::log((double)kdcount)/std::log(2.0);

    std::cout << "   count    : " << kdcount << std::endl;
    std::cout << "   depth    : " << kdshortest << " (min) / " << kddeepest << " (max)" << std::endl;
    std::cout << "   cost     : " << kdcost << " (est.) / " << kdoptimum_cost << " (opt)" << std::endl;
  }

  int cps = int(float(points.size()) / float(kdt_create_ms/1000.0f));
  std::cout << "   create   : " << kdt_create_ms << " ms / " << sl::human_readable_quantity(int(cps)) << " insert/second" << std::endl;
  {
    std::cout << "---------- Queries for points in kdtree " << std::endl;
    std::vector<sl::point3f> queries(1000);
    for (std::size_t k = 0;
	 k < queries.size();
	 ++k) {
      int j = rand() % (points.size());
      queries[k] = points[j];      
    }

    for (std::size_t q = 0; q<2; ++q) {
      std::vector<sl::point3f> queries(1000);
      if (q==0) {
	std::cout << "---------- Queries for points in kdtree " << std::endl;
	for (std::size_t k = 0;
	     k < queries.size();
	     ++k) {
	  int j = rand() % (points.size());
	  queries[k] = points[j];      
	}
      } else {
	std::cout << "---------- Queries for points far from kdtree bbox " << std::endl;
	for (std::size_t k = 0;
	     k < queries.size();
	     ++k) {
	  queries[k] = pmin + 2.0f*(pmax-pmin);
	}
      }

      for (std::size_t k = 1;
	   k <= 32;
	   k *= 2) {
	for (std::size_t e=0; e<5; e+=2) {
	  
	  for (std::size_t mode = 0; mode<3; ++mode) {
	    bool is_prioritized = mode > 0;
	    std::size_t max_visits = mode < 2 ? std::size_t(-1) : std::size_t(100);
	    float eps = 0.1*float(e);
	    float bot = 1e30f;
	    std::vector<kdt_t::value_t>    knn_val(k);
	    
	    kdt.set_is_priority_search_enabled(is_prioritized);
	    kdt.set_priority_search_max_visited_count(max_visits); 
	    kdt.stat_reset();
	    
	    std::size_t cnt = 0;
	    std::size_t maxv = 0;
	    float maxd = 0.0f;
	    sl::cpu_time_clock clk;
	    clk.restart();
	    while (clk.elapsed() < 0.5*sl::time_duration::one_second()) {
	      
	      for (std::size_t i=0; i<queries.size(); i++) {
		++cnt;
		std::size_t v0 =  kdt.stat_knn_visited_count();
		kdt.k_approximately_nearest_neighbors_in(knn_val,
							 queries[i],
							 eps,
							 bot,
							 k);
		maxv = std::max(maxv, kdt.stat_knn_visited_count()-v0);
		if (!knn_val.empty()) maxd = std::max(maxd,
						      knn_val.back().distance_to(queries[i]));
	      }
	    }
	    
	    float qps = cnt / (clk.elapsed().as_microseconds()/1E6f);
	    float vpq = kdt.stat_knn_visited_count()/float(kdt.stat_knn_query_count());
	    
	    std::cout << "   knn query: " << (is_prioritized ? "Prioritized" : "Depth-first") << (max_visits ==std::size_t(-1) ? std::string(" No Limit ") : std::string(" Limit=")+sl::to_string(max_visits)) << " EPS = " << eps << " K = " << k << " /  " << sl::human_readable_quantity(int(qps)) << " queries/second / " << sl::human_readable_quantity(int(vpq)) << " visits/query (max=" << maxv << ")/ maxd=" << maxd << std::endl;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }
  
  int kdt_clear_ms = 0;
  {
    sl::cpu_time_clock clk;
    clk.restart();
    kdt.clear();
    kdt_clear_ms = (int)(clk.elapsed().as_microseconds()/1000);
  }
  std::cout << "   clear    : " << kdt_clear_ms << " ms" << std::endl;
}


static void benchmark_distribution(const char* id,
				   const std::vector<sl::point3f>& points) {
  sl::kdtree<3,float,sl::point3f>            kdt;

  std::cout << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "DISTRIBUTION: " << id << " (" << points.size() << " points)" << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  benchmark_kdt("RANDOMIZED KD TREE", kdt, points);
}

static void test_features(const char* id,
		 const std::vector<sl::point3f>& points) {
  typedef sl::kdtree<3,float,sl::point3f> kdt_t;
  
  kdt_t kdt;

  std::cout << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "TEST: " << id << " (" << points.size() << " points)" << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  
  kdt.insert(points.begin(), points.end());

  std::size_t count = 0;
  for (kdt_t::iterator it = kdt.begin(); it!= kdt.end(); ++it) {
    std::cout << count << ": " << *it << std::endl;
    ++count;
  }
  std::cout << "Count = " << count << " vs. " << kdt.count() << std::endl;
}

int main() {
  test_features("Test", uniform_points(10));
#if 1
  benchmark_distribution("Uniform", uniform_points(1000000));
  benchmark_distribution("Z-Ordered", z_ordered_points(1000000));
  benchmark_distribution("Ordered", ordered_points(1000000));
#endif
}
