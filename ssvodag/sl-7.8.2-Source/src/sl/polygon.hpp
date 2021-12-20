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
#ifndef SL_POLYGON_HPP
#define SL_POLYGON_HPP

#include <sl/fixed_size_point.hpp>

namespace sl {

  /**
   *  Is poly a CCW-ordered convex polygon?
   */
  template <class T_POLYGON>
  bool is_ccw_convex_polygon(const T_POLYGON& poly) {
    bool result = false;
    const size_t n = poly.size();
    if (n>=2) {
      result = true;
      for (size_t i =0; i<n && result; ++i) {
	point_ordering_t o = 
	  order(poly[(i+n+0)%n], 
		poly[(i+n-1)%n],
		poly[(i+n+1)%n],
		epsilon(poly[0][0]));
	  
	result = (o == ORDER_CURVED_CW) || (o == ORDER_COLLINEAR_INBETWEEN);
      }
    }
    return result;
  }
 
  /**
   *  Is poly a CW-ordered convex polygon?
   */
  template <class T_POLYGON>
  bool is_cw_convex_polygon(const T_POLYGON& poly) {
    bool result = false;
    const size_t n = poly.size();
    if (n>=2) {
      result = true;
      for (size_t i =0; i<n && result; ++i) {
	point_ordering_t o = 
	  order(poly[(i+n+0)%n], 
		poly[(i+n-1)%n],
		poly[(i+n+1)%n],
		epsilon(poly[0][0]));
	result = (o == ORDER_CURVED_CCW || o == ORDER_COLLINEAR_INBETWEEN);
      }
    }
    return result;
  }

  /**
   *  Is poly a CW-ordered convex polygon?
   */
  template <class T_POLYGON>
  bool is_convex_polygon(const T_POLYGON& poly) {
    return 
      is_ccw_convex_polygon(poly) || 
      is_cw_convex_polygon(poly);
  }

  /**
   *  Is point p inside the CCW-ordered convex polygon p?
   */
  template <class T_POLYGON,  class T_VALUE>
  bool point_in_ccw_convex_polygon(const T_POLYGON& poly,
				   const sl::fixed_size_point<2,T_VALUE>& p) {
    typedef T_VALUE value_t;

    bool result = false;
    const size_t n = poly.size();
    if (n>=3) {
      result = true;
      for (size_t i1 = 0, i0 = n-1; i1 < n && result; i0 = i1, ++i1) {
	const value_t n_x = poly[i1][1] - poly[i0][1];
	const value_t n_y = poly[i0][0] - poly[i1][0];
	const value_t d_x = p[0] - poly[i0][0];
	const value_t d_y = p[1] - poly[i0][1];
	const bool out = sl::is_positive(n_x*d_x + n_y*d_y);
	result = !out;
      }
    }
    return result;
  }
    
} // namespace sl

#endif

