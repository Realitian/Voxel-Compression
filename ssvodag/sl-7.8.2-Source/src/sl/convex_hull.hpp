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
#ifndef SL_CONVEX_HULL_HPP
#define SL_CONVEX_HULL_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/utility.hpp>
#include <vector>
#include <list>
#include <set>
#include <algorithm>

namespace sl {

  /// Connectivity of an n-dimensional simplex
  template <size_t DIMENSION> 
  class simplex_connectivity: public fixed_size_array<DIMENSION, size_t> {
  public:
    enum { dimension = DIMENSION };

    typedef simplex_connectivity<dimension>    self_t;
    typedef fixed_size_array<dimension,size_t> super_t;
    typedef size_t                             value_t;
    
  public: // Type constraints
    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  public: // Creation & destruction
    
    /// Default init (zero)
    inline simplex_connectivity() {
    }

    /// Fast init (garbage), handle with care!
    inline simplex_connectivity(const tags::not_initialized tag): super_t(tag) {
      // Garbage memory contents!
    }

    /// Explicit init from components (1D)
    inline simplex_connectivity(size_t i0): super_t(tags::not_initialized()) {
      SL_REQUIRE("Good dimension", dimension == 1);
      (this->storage_)[0] = i0;
    }

    /// Explicit init from components (2D)
    inline simplex_connectivity(size_t i0, 
				size_t i1): super_t(tags::not_initialized()) {
      SL_REQUIRE("Good dimension", dimension == 2);
      (this->storage_)[0] = i0;
      (this->storage_)[1] = i1;
    }

    /// Explicit init from components (3D)
    inline simplex_connectivity(size_t i0, 
				size_t i1,
				size_t i2): super_t(tags::not_initialized()) {
      SL_REQUIRE("Good dimension", dimension == 3);
      (this->storage_)[0] = i0;
      (this->storage_)[1] = i1;
      (this->storage_)[2] = i2;
    }

    /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1; (fill) and 
     *  v = 1, 2, 3; (component init)
     */
    inline manifest_array1d_initializer<self_t&,value_t> operator= (value_t v) {
      return manifest_array1d_initializer<self_t&,value_t>(*this, dimension, v);
    }
    
  }; // simplex_connectivity

} // namespace sl

namespace sl {

  /**
   *  Wrapper for converting a simplex-point representation of
   *  a 2D polygon/polyline to an array representation.
   */
  template <class T_VALUE, class T_POINT_ARRAY>
  class simplex_to_poly_wrapper {
  public:
    enum { dimension = 2 };

    typedef T_VALUE                             value_t;
    typedef fixed_size_point<dimension,value_t> point_t; 

    typedef simplex_connectivity<dimension>      edge_t;

  protected:
    const std::vector<edge_t>&          edges_;
    const T_POINT_ARRAY&                points_;
  public:
    
    /**
     *  Create a wrapper for a polygon defined by edges and points.
     *  Internally, only *references* to those structures are kept;
     *  make sure that the lifesize of edges and points exceeds the
     *  lifesize of this wrapper!
     */
    inline simplex_to_poly_wrapper(const std::vector<edge_t>& edges,
				    const T_POINT_ARRAY& points):
      edges_(edges), points_(points) {
    }

    /// The number of points
    inline size_t size() const {
      return edges_.size();
    }

    /// The i-th point
    inline point_t operator[](size_t i) const {
      SL_REQUIRE("Good index", i<size());
      SL_REQUIRE("Consistent structure", edges_[i][0] < points_.size());

      return points_[edges_[i][0]];
    }
      
  };

} // namespace sl

namespace sl {

  /**
   *  Base class for N-dimensional convex hull builders.
   *  The point array class should support size(), returning
   *  the number of points, and operator[], returning
   *  a point_t providing the i-th vertex
   *  position.
   */
  template <size_t DIMENSION, class T, class T_POINT_ARRAY>
  class convex_hull_builder_base {
  public:
    enum { dimension = DIMENSION };

    typedef T                                      value_t; 
    typedef convex_hull_builder_base<DIMENSION, T, T_POINT_ARRAY> self_t;

    typedef fixed_size_point<DIMENSION,T>          point_t;
    typedef simplex_connectivity<DIMENSION>        simplex_t;
    typedef std::vector<simplex_t>                 hull_t;
    typedef T_POINT_ARRAY                          point_array_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  protected: // State

    hull_t last_hull_;

  public: // Creation & Destruction

    convex_hull_builder_base() { 
    }

    virtual inline ~convex_hull_builder_base() {
    }

  public: // Construction
      
    /**
     *  Build the hull from the given polygon. 
     */
    virtual void build(const point_array_t& points) = 0;
    
    /// The hull constructed by last build()
    virtual const hull_t& last_hull() { 
      return last_hull_; 
    }

  }; // class convex_hull_builder_base


  /// N-dimensional convex hull builder - currently not implemented...
  template <size_t DIMENSION, class T, class T_POINT_ARRAY> 
  class convex_hull_builder: public convex_hull_builder_base<DIMENSION,T, T_POINT_ARRAY> {
  public:

    typedef T                                      value_t; 
    typedef convex_hull_builder_base<DIMENSION, T, T_POINT_ARRAY> super_t;
    typedef convex_hull_builder<DIMENSION, T, T_POINT_ARRAY> self_t;

    typedef fixed_size_point<DIMENSION,T>          point_t;
    typedef simplex_connectivity<DIMENSION>        simplex_t;
    typedef std::vector<simplex_t>                 hull_t;
    typedef T_POINT_ARRAY                          point_array_t;

  public:

    /// Default init
    inline convex_hull_builder() {}
 
    /// Destructor
    virtual inline ~convex_hull_builder() {
    }

    /**
     *  Build the hull from the given point array.
     */
    virtual void build(const point_array_t& points) {
      SL_FAIL("Convex hull not implemented for general dimension");
    }

  }; // convex_hull_builder

  
  /**
   *  Two-dimensional CCW-ordered convex hull builder. 
   *  Algorithms are generally O(N) in storage and 
   *  O(N*log(N)) in time complexity.
   */
  template <class T, class T_POINT_ARRAY> 
  class convex_hull_builder<2,T, T_POINT_ARRAY>: public convex_hull_builder_base<2,T,T_POINT_ARRAY> {
  public:

    enum { dimension = 2 };

    typedef T                                      value_t; 
    typedef convex_hull_builder_base<dimension, T, T_POINT_ARRAY> super_t;
    typedef convex_hull_builder<dimension, T, T_POINT_ARRAY>      self_t;

    typedef fixed_size_point<dimension,T>          point_t;
    typedef simplex_connectivity<dimension>        simplex_t;
    typedef std::vector<simplex_t>                 hull_t;
    typedef T_POINT_ARRAY                          point_array_t;

  public: // Constraints

    /// Default init
    inline convex_hull_builder() {
    }

    /// Destructor
    virtual inline ~convex_hull_builder() {
    }

    /**
     *  Build the hull from the given polygon. The
     *  point array class should support size(), returning
     *  the number of points, and operator[], returning
     *  a point_t providing the i-th vertex
     *  position.
     */
    virtual void build(const point_array_t& points) {
#if 1
      build_chain(points);
#else
      build_incremental(points);
#endif
      remove_collinear(points);
    }

  protected: // Helpers
    
    /**
     *  Ordering of points, first by x, then by y position
     */
    class xy_ordering: public std::binary_function<size_t, size_t, bool> {
    protected:
      const point_array_t& points_;
      const value_t eps_;
    public:
      xy_ordering(const point_array_t& pts,
		  const value_t& eps): points_(pts), eps_(eps) {
      }

      inline bool operator()(size_t i,
			     size_t j) const {
#if 1
	bool result = false;

	const value_t x_i = points_[i][0];
	const value_t x_j = points_[j][0];

	if (x_i + eps_ < x_j) {
	  result = true;
	} else if (x_i < x_j + eps_) {
	  const value_t y_i = points_[i][1];
	  const value_t y_j = points_[j][1];

	  result = (y_i + eps_ < y_j);
	}

	return result;
#else
	const value_t x_i = points_[i][0];
	const value_t x_j = points_[j][0];
	return 
	  (x_i + eps_ < x_j) ||
	  ((x_i < x_j + eps_) &&
	   (points_[i][1] + eps_ < points_[j][1]));
#endif
      }
    };

    /**
     *  Epsilon-equality of points
     */
    class xy_equal: public std::binary_function<size_t, size_t, bool> {
    protected:
      const point_array_t& points_;
      const value_t eps_;
    public:
      xy_equal(const point_array_t& pts,
	       const value_t& eps): points_(pts), eps_(eps) {
      }

      inline bool operator()(size_t i,
			     size_t j) const {

	const value_t x_i = points_[i][0];
	const value_t x_j = points_[j][0];

	bool result = (x_i + eps_ >= x_j) && (x_i <= x_j + eps_);
	if (result) {
	  const value_t y_i = points_[i][1];
	  const value_t y_j = points_[j][1];

	  result = (y_i + eps_ >= y_j) && (y_i <= y_j + eps_);
	}

	return result;
      }
    };

    // tests if point P2 is Left|On|Right of the line P0 to P1.
    // returns: >0 for left, 0 for on, and <0 for right of the line.
    inline value_t is_left(const point_t& P0, 
			   const point_t& P1, 
			   const point_t& P2) const {
      return (P1[0] - P0[0])*(P2[1] - P0[1]) - (P2[0] - P0[0])*(P1[1] - P0[1]);
    }

    
  public:

    /**
     *  Build the 2D hull using the O(N*log(N)) monotone chain algorithm
     *  A.M. Andrew, "Another Efficient Algorithm for Convex Hulls in Two Dimensions", 
     *  Info. Proc. Letters 9, 216-219 (1979)
     */
    virtual void build_chain(const point_array_t& points) {
      const value_t eps = scalar_math<value_t>::epsilon();

      (this->last_hull_).clear();

      size_t n = points.size();
      if (n > 0) {
	// Build index array and sort by x component
	// Note: most of the time is spent in the sorting routine!
	std::vector<size_t> point_directory(n);
	for (size_t i=0; i<n; ++i) {
	  point_directory[i] = i;
	}
	std::sort(point_directory.begin(), point_directory.end(), xy_ordering(points, eps));

	// Remove duplicate points
	point_directory.erase(std::unique(point_directory.begin(), 
					  point_directory.end(), 
					  xy_equal(points, eps)),
			      point_directory.end());

	// Get actual point number and build hull
	const size_t N = point_directory.size();
	SL_CHECK("At least one point",N > 0);
	std::vector<size_t> hull_vertex_stack;
	
	const value_t xmin = points[point_directory[0]][0];
	const value_t xmax = points[point_directory[N-1]][0];

	// Get the indices of points with min x-coord and min|max y-coord
	size_t minmin = 0;
	size_t minmax = N-1;
	for (size_t i=1; i<N; ++i) {
	  if (points[point_directory[i]][0] != xmin) {
	    minmax = (i-1);
	    break;
	  }
	}
	if (minmax == N-1) {
	  // degenerate case: all x-coords == xmin
	  hull_vertex_stack.push_back(point_directory[minmin]);
	  if (points[point_directory[minmax]][1] != points[point_directory[minmin]][1]) {
	    hull_vertex_stack.push_back(point_directory[minmax]);
	  }
	} else {
	  // Get the indices of points with max x-coord and min|max y-coord
	  size_t maxmax = N-1;
	  size_t maxmin = 0;
	  //for (size_t i=N-2; i>=0; --i) {
      SL_CHECK("Not exceeded max point count", N < ((uint32_t)-1)>>1);
      for (int32_t i=N-2; i>=0; --i) {
	    if (points[point_directory[i]][0] != xmax) {
	      maxmin = (size_t)(i+1);
	      break;
	    }
	  }
	  
	  // Compute the lower hull 
	  hull_vertex_stack.push_back(point_directory[minmin]);
	  for (size_t i=minmax+1; i <= maxmin; ++i) {
	    if ((i < maxmin) && 
		(sl::is_non_negative(is_left(points[point_directory[minmin]],
					     points[point_directory[maxmin]],
					     points[point_directory[i]])))) {
	      // Ignore: above or on lower line
	    } else {
	      size_t sz = hull_vertex_stack.size();
	      while ( sz >= 2 ) {
		if (sl::is_positive(is_left(points[hull_vertex_stack[sz-2]],
					    points[hull_vertex_stack[sz-1]],
					    points[point_directory[i]]))) {
		  break;
		} else {
		  --sz;
		  hull_vertex_stack.pop_back();
		}
	      }
	      hull_vertex_stack.push_back(point_directory[i]);
	    }
	  }
	  
	  // Next, compute the upper hull on the stack above the bottom hull
	  if (maxmax != maxmin) {  
	    hull_vertex_stack.push_back(point_directory[maxmax]);
	  }
	  size_t up_sz = hull_vertex_stack.size();
	  for (int i=(int)maxmin-1; i>=(int)minmax; --i) {
	    if ((i > (int)minmax) && 
		(sl::is_non_negative(is_left(points[point_directory[maxmax]],
					     points[point_directory[minmax]],
					     points[point_directory[i]])))) {
	      // Ignore: below or on the upper line
	    } else {
	      size_t sz = hull_vertex_stack.size();
	      while ( sz > up_sz) {
		if (sl::is_positive(is_left(points[hull_vertex_stack[sz-2]],
					    points[hull_vertex_stack[sz-1]],
					    points[point_directory[i]]))) {
		  break;
		} else {
		  --sz;
		  hull_vertex_stack.pop_back();
		}
	      }
	      hull_vertex_stack.push_back(point_directory[i]);
	    }
	  }
	  if (minmax != minmin) {      // if distinct xmax points
	    if (hull_vertex_stack.back() != point_directory[minmax]) {
	      hull_vertex_stack.push_back(point_directory[minmax]);
	    }
	  }

	  // Cleanup
	  while (hull_vertex_stack.back() == hull_vertex_stack.front()) {
	    hull_vertex_stack.pop_back();
	  }
	}
	
	// Convert hull to simplex form
	std::vector<size_t>::iterator       it     = hull_vertex_stack.begin();
	const std::vector<size_t>::iterator it_end = hull_vertex_stack.end();
	while (it != it_end) {
	  const size_t i0 = *it;
	  ++it; 
	  const size_t i1 = ((it == it_end) ? (hull_vertex_stack.front()) : (*it));
	  (this->last_hull_).push_back(simplex_t(i0, i1));
	}
      }
    }
    

    /**
     *  Build the 2D hull using the O(N*log(N)) incremental
     *  algorithm.
     *  M. Kallay, "The Complexity of Incremental Convex Hull 
     *  Algorithms in Rd", Info. Proc. Letters 19, 197 (1984)
     */
    virtual void build_incremental(const point_array_t& points) {
      const value_t eps = scalar_math<value_t>::epsilon();

      (this->last_hull_).clear();

      size_t n = points.size();
      if (n > 0) {
	// Build index array and sort by x component
	// Note: most of the time is spent in the sorting routine!
#if 1
	// Build an array and sort it
	typedef std::vector<size_t> point_directory_t;

	point_directory_t point_directory(n);
	for (size_t i=0; i<n; ++i) {
	  point_directory[i] = i;
	}
	std::sort(point_directory.begin(), point_directory.end(), xy_ordering(points, eps));

	// Remove duplicate points
	point_directory.erase(std::unique(point_directory.begin(), 
					  point_directory.end(), 
					  xy_equal(points, eps)),
			      point_directory.end());
#else
	// Construct a sorted set - this seems actually slower
	// than sorting the array
	typedef std::set<size_t,xy_ordering> point_directory_t;
	point_directory_t point_directory(xy_ordering(points,eps));
	for (size_t i=0; i<n; ++i) {
	  point_directory.insert(i);
	}
#endif

	// Add points to hull incrementally
	std::list<size_t> hull_vertex_list;
	{
	  point_directory_t::iterator       point_it   = point_directory.begin();
	  const point_directory_t::iterator point_it_end = point_directory.end();

	  SL_CHECK("At least one point", point_it != point_it_end);
	  hull_vertex_list.push_back(*point_it);
	  ++point_it;
	  if (point_it != point_it_end) {
	    hull_vertex_list.push_back(*point_it);
	    ++point_it;
	    while (hull_vertex_list.size() == 2 && point_it != point_it_end) {
	      merge_1d(*point_it, points, hull_vertex_list);
	      ++point_it;
	    }
	    while (point_it != point_it_end) {
	      merge_2d(*point_it, points, hull_vertex_list);
	      ++point_it;
	      if (SL_TRACE_LEVEL > 1) {
		SL_TRACE_OUT(1) << "current_hull: ";
		for (std::list<size_t>::iterator it = hull_vertex_list.begin();
		     it !=  hull_vertex_list.end();
		     ++it) {
		  SL_TRACE_OUT(1) << " " << *it;
		}
		SL_TRACE_OUT(1) << std::endl;
	      }
	    }
	  }
	}


	// Convert hull to simplex form
	std::list<size_t>::iterator       it     = hull_vertex_list.begin();
	const std::list<size_t>::iterator it_end = hull_vertex_list.end();
	while (it != it_end) {
	  const size_t i0 = *it;
	  ++it; 
	  const size_t i1 = ((it == it_end) ? (hull_vertex_list.front()) : (*it));
	  (this->last_hull_).push_back(simplex_t(i0, i1));
	}
      }
    }

  protected: 

    void remove_collinear(const point_array_t& points) {
      const size_t N = (this->last_hull_).size();
      const value_t eps = scalar_math<value_t>::epsilon();
      if (N > 2 ) {
	std::vector<size_t> new_vertices;
	new_vertices.reserve(N);
	for (size_t i0 = N-1, i1 = 0, i2 = 1; i1 < N; /**/) {
	  switch(order(points[(this->last_hull_)[i0][0]], 
		       points[(this->last_hull_)[i1][0]], 
		       points[(this->last_hull_)[i2][0]], 
		       eps)) {
	  case ORDER_CURVED_CCW:
	  case ORDER_CURVED_CW:
	    new_vertices.push_back((this->last_hull_)[i1][0]);
	    break;
	  default:
	    // COLLINEAR, REMOVE
	    break;
	  };
	  i0 = i1;
	  ++i1;
	  ++i2;
	  if (i2 == N) i2 = 0;
	}
	const size_t N2 = new_vertices.size();
	if (N2 != N) {
	  SL_TRACE_OUT(1) << 
	    " Removing " << N-N2 << " collinear vertices from hull of size " << N << std::endl;
	  (this->last_hull_).clear();
	  (this->last_hull_).reserve(N2);
	  for (size_t i=0; i<N2; ++i) {
	    (this->last_hull_).push_back(simplex_t(new_vertices[i], 
                                                 new_vertices[(i+1)%N2]));
	  }
	}
      }
    }

    void merge_1d(size_t  idx,
		  const point_array_t& points,
		  std::list<size_t>& hull_vertex_list) {
      SL_REQUIRE("Good hull size", hull_vertex_list.size() == 2);
      
      const value_t eps = scalar_math<value_t>::epsilon();

      std::list<size_t>::iterator it_0 = hull_vertex_list.begin();
      std::list<size_t>::iterator it_1 = ++hull_vertex_list.begin();
      
      switch (order(points[idx],
		    points[*it_0],
		    points[*it_1],
		    eps)) {
      case ORDER_CURVED_CCW:
	hull_vertex_list.push_back(*it_1);
	*it_1 = *it_0;
	*it_0 = idx;
	break;
      case ORDER_CURVED_CW:
	hull_vertex_list.push_back(*it_0);
	*it_0 = idx;
	break;
      case ORDER_COLLINEAR_BEFORE:
	*it_0 = idx;
	break;
      case ORDER_COLLINEAR_AFTER:
	*it_1 = idx;
	break;
      case ORDER_COLLINEAR_INBETWEEN:
	break;
      default:
	SL_FAIL("Bad ordering tag");
	break;
      }
      SL_ENSURE("Good hull size", hull_vertex_list.size() >= 2);
    }

    void merge_2d(size_t  idx,
		  const point_array_t& points,
		  std::list<size_t>& hull_vertex_list) {
      SL_REQUIRE("Good hull size", hull_vertex_list.size() >= 3);
      
      const value_t eps = scalar_math<value_t>::epsilon();

      const std::list<size_t>::iterator it_first = hull_vertex_list.begin();
      std::list<size_t>::iterator it_last  = hull_vertex_list.end(); --it_last;
      
      // Search counterclockwise for last visible vertex
      std::list<size_t>::iterator       it_hi  = cyclic_pred(it_first, it_first, it_last);
      {
	bool found = false;
	bool done  = false;

	while (!(found || done)) {
	  it_hi = cyclic_succ(it_hi, it_first, it_last);
	  point_ordering_t o = order(points[idx],
				     points[*it_hi],
				     points[*cyclic_succ(it_hi, it_first, it_last)],
				     eps);
	  
	  found = (o == ORDER_CURVED_CCW || o == ORDER_COLLINEAR_BEFORE);
	  done  = (!found) && (o != ORDER_CURVED_CW);
	}
	if (done) return;
      }
	  
      // Search clockwise for first visible vertex
      std::list<size_t>::iterator       it_lo  = cyclic_succ(it_first, it_first, it_last);
      {
	bool found = false;
	bool done  = false;
	while (!(found || done)) {
	  it_lo = cyclic_pred(it_lo, it_first, it_last);
	  point_ordering_t o = order(points[idx],
				     points[*cyclic_pred(it_lo, it_first, it_last)],
				     points[*it_lo],
				     eps);
	  
	  found = (o == ORDER_CURVED_CCW || o == ORDER_COLLINEAR_AFTER);
	  done  = (!found) && (o != ORDER_CURVED_CW);
	}
	if (done) return;
      }

      SL_TRACE_OUT(1) << "LO : " << *it_lo << std::endl;
      SL_TRACE_OUT(1) << "HI : " << *it_hi << std::endl;
	
      // Construct CCW hull by removing invisible vertices and adding new one
      if (it_lo == it_first) {
	++it_lo;
	if (it_hi != it_first) {
	  const size_t old_first = *it_first;
	  hull_vertex_list.erase(it_first, it_hi);
	  hull_vertex_list.push_back(old_first);
	}
      } else {
	if (it_hi != it_first) {
	  hull_vertex_list.erase(it_first, it_hi);
	}
	++it_lo;
	if (it_lo != hull_vertex_list.end()) {
	  hull_vertex_list.erase(it_lo, hull_vertex_list.end());
	}
      }
      hull_vertex_list.push_front(idx);
 
      SL_ENSURE("Good hull size", hull_vertex_list.size() >= 3);
    }

  }; // convex_hull_builder
  

}; // namespace sl


#endif
