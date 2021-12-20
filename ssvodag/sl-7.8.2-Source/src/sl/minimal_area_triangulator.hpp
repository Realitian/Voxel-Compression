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

#ifndef SL_MINIMAL_AREA_TRIANGULATOR_HPP
#define SL_MINIMAL_AREA_TRIANGULATOR_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/connectivity.hpp>
#include <sl/math.hpp>
#include <vector>
#include <cassert>

namespace sl {

  /**
   *  A minimal area triangulator of (small) 
   *  polygons
   */
  template <std::size_t DIMENSION, class T>
  class minimal_area_triangulator {
  public:
    enum { dimension = DIMENSION };
    typedef T value_t;
    typedef sl::fixed_size_point<dimension,value_t> point_t;
    typedef sl::fixed_size_vector<sl::column_orientation,dimension,value_t> vector_t;
    typedef triangle_connectivity triidx_t;
  public:

    minimal_area_triangulator() {
    }

    ~minimal_area_triangulator() {
    }

    void triangulation_in(std::vector<triidx_t>& triangles,
			  const std::vector<point_t>& vertices) const {
      const std::size_t N = vertices.size();
      if (N<3) {
	// Degenerate
	triangles.clear();
      } else if (N == 3) {
	// Triangle
	triangles.resize(1);
	triangles[0] = triidx_t(0,1,2);
      } else if (N == 4) {
	// Square
	triidx_t trial_triangles[2][2];
	value_t  trial_area[2];
	trial_triangles[0][0] = triidx_t(0,1,2);
	trial_triangles[0][1] = triidx_t(2,3,0);
	trial_triangles[1][0] = triidx_t(0,1,3);
	trial_triangles[1][1] = triidx_t(3,1,2);
	for(int i=0; i<2; ++i) {
	  trial_area[i] = 0.0;
	  for(int j=0; j<2; ++j) {
	    trial_area[i] += sl::triangle_area(vertices[trial_triangles[i][j][0]],
					       vertices[trial_triangles[i][j][1]],
					       vertices[trial_triangles[i][j][2]]);
	  }
	}	
	triangles.resize(2);
	if (trial_area[0]<trial_area[1]) {
	  triangles[0] = trial_triangles[0][0];
	  triangles[1] = trial_triangles[0][1];
	} else {
	  triangles[0] = trial_triangles[1][0];
	  triangles[1] = trial_triangles[1][1];
	}
      } else {
	// General case
	std::vector<value_t> best_area(N*N, value_t(-1.0));
	std::vector<int>     best_idx(N*N, int(-1)); 

	update_area_idx_in(best_area,
			   best_idx,
			   0,1,
			   vertices);
	triangles.clear();
	append_triangulation_in(triangles,
				best_idx,
				N, 0, 1);
      }
    }
    
  protected:
 
    void update_area_idx_in(std::vector<value_t>& best_area,
			    std::vector<int>& best_idx,
			    std::size_t i,
			    std::size_t j,
			    const std::vector<point_t>& vertices) const {
      const std::size_t N=vertices.size();
      const std::size_t idx=i*N+j;

      std::size_t ii=(i<j)?i+N:i;

      if(j+1>=ii) {
	best_area[idx]=0;
	//best_idx[idx]=-1; // FIXME
      } else if (best_idx[idx]!=-1) {
	// Already computed
      } else {
	// Compute
	value_t best_a=sl::scalar_math<value_t>::finite_upper_bound();
	int     best_i=-1;
	for(std::size_t r=j+1;r<ii;++r) {
	  std::size_t rr=r%N;
	  std::size_t idx1=i*N+rr;
	  std::size_t idx2=rr*N+j;
	  
	  value_t this_a = sl::triangle_area(vertices[rr], 
					     vertices[i], 
					     vertices[j]);
	  if (best_area[idx1]>=0) {
	    this_a+=best_area[idx1];

	    if (this_a<best_a) {
	      if (best_area[idx2]<0) {
		update_area_idx_in(best_area,
				   best_idx,
				   rr, j,
				   vertices);
	      }
	      this_a+=best_area[idx2];
	    }
	  } else {
	     if(best_area[idx2]<0) {
	       update_area_idx_in(best_area,
				 best_idx,
				 rr, j,
				 vertices);
	    }
	    this_a+=best_area[idx2];

	    if (this_a<best_a) {
	      update_area_idx_in(best_area,
				 best_idx,
				 i, rr,
				 vertices);
	      this_a+=best_area[idx1];
	    }
	  }

	  if(this_a<best_a) {
	    best_a=this_a;
	    best_i=int(rr);
	  }
	}

	best_area[idx]=best_a;
	best_idx[idx]=best_i;
      }  
    }

    void append_triangulation_in(std::vector<triidx_t>& triangles,
				 const std::vector<int>& best_idx,
				 std::size_t N,
				 std::size_t i,
				 std::size_t j) const {
      const std::size_t ii=(i<j)?(i+N):i;
      if (j+1<ii) {
	const std::size_t idx=i*N+j;
	const int k=best_idx[idx];
	if (k>=0) {
	  triangles.push_back(triidx_t(i,j,k));
	  append_triangulation_in(triangles, best_idx, N, i, k);
	  append_triangulation_in(triangles, best_idx, N, k, j);
	}
      }
    }
  }; // class minimal_area_triangulator

} // namespace sl

#endif
