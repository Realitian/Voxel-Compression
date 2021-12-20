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
#ifndef SL_INTERPOLATION_HPP
#define SL_INTERPOLATION_HPP

#include <sl/numeric_traits.hpp>
#include <sl/math.hpp>
#include <functional>
#include <map>
#include <vector>

namespace sl {

  /// linear interpolation of (p1, p2) evaluated at t
  template <class T_CONTROL_POINT, class T_PARAMETER>
  inline T_CONTROL_POINT lerp(const T_CONTROL_POINT& p1, 
			      const T_CONTROL_POINT& p2, 
			      const T_PARAMETER& t) {
    // Default: we assume that the control point class 
    // implements it!
    return p1.lerp(p2, t); 
  }

#define SL_DECLARE_LERP(T1,T2) \
  inline T1 lerp(const T1& p1, const T1& p2, const T2& t) { \
    typedef SL_FLOATTYPE(T1) float_t; \
    return T1(float_t(p1) + (float_t(p2) - float_t(p1)) * float_t(t)); \
  }

  SL_DECLARE_LERP(float,float)
  SL_DECLARE_LERP(float,double)
  SL_DECLARE_LERP(double,float)
  SL_DECLARE_LERP(double,double)
  SL_DECLARE_LERP(int8_t,float)
  SL_DECLARE_LERP(int8_t,double)
  SL_DECLARE_LERP(int16_t,float)
  SL_DECLARE_LERP(int16_t,double)
  SL_DECLARE_LERP(int32_t,float)
  SL_DECLARE_LERP(int32_t,double)
#if HAVE_INTMAX_BITS >= 64
  SL_DECLARE_LERP(int64_t,float)
  SL_DECLARE_LERP(int64_t,double)
#endif

  // ... and so on for other basic data types ...

  /// Function object interface to linear interpolation
  template <class T_CONTROL_POINT, class T_PARAMETER>
  class lerp_op : 
    public std::binary_function<T_CONTROL_POINT,T_CONTROL_POINT,T_CONTROL_POINT> {
  protected:
    const T_PARAMETER t;
  public:
    inline lerp_op(const T_PARAMETER& t_arg = T_PARAMETER(0)): t(t_arg) {};
    inline T_CONTROL_POINT operator()(const T_CONTROL_POINT& x,
                                      const T_CONTROL_POINT& y) const { 
      return lerp(x,y,t); 
    }
  }; // lerp_op

  /// Cubic Bezier interpolation
  template <class T_CONTROL_POINT, class T_PARAMETER>
  T_CONTROL_POINT bezier_interpolation(const T_CONTROL_POINT &p00, 
				       const T_CONTROL_POINT &p01, 
				       const T_CONTROL_POINT &p02, 
				       const T_CONTROL_POINT &p03 , 
				       const T_PARAMETER& t) {
    T_CONTROL_POINT p10 = lerp(p00, p01, t);
    T_CONTROL_POINT p11 = lerp(p01, p02, t);
    T_CONTROL_POINT p12 = lerp(p02, p03, t);
    
    T_CONTROL_POINT p20 = lerp(p10, p11, t);
    T_CONTROL_POINT p21 = lerp(p11, p12, t);
    
    return lerp(p20, p21, t);
  }

  /// Conversion of a Catmull-Rom control point to a bezier one
  template <class T_CONTROL_POINT, class T_PARAMETER>
  T_CONTROL_POINT catmull_rom_to_bezier_control_point(const T_CONTROL_POINT &pa, const T_PARAMETER& ta,
						      const T_CONTROL_POINT &p0, const T_PARAMETER& t0,
						      const T_PARAMETER& t1) {
    T_PARAMETER ta0 = t0 - ta;
    if (sl::is_zero(ta0)) {
      return p0;
    } else {
      T_PARAMETER three_ta0 = ta0+ta0+ta0;
      return lerp(pa, p0, (three_ta0 + t1 - t0)/three_ta0);
    }
  }

  /// Linear interpolation 
  template <class T_CONTROL_POINT, class T_PARAMETER>
  T_CONTROL_POINT  linear_interpolation(const T_CONTROL_POINT &p0, const T_PARAMETER& t0,
                                        const T_CONTROL_POINT &p1, const T_PARAMETER& t1,
                                        const T_PARAMETER& t) {
    T_PARAMETER t01 = t1 - t0;
    if (sl::is_zero(t01)) {
      return p0;
    } else {
      return lerp(p0,p1,(t - t0)/t01);
    }
  }

  /// Catmull-Rom interpolation 
  template <class T_CONTROL_POINT, class T_PARAMETER>
  T_CONTROL_POINT  catmull_rom_interpolation(const T_CONTROL_POINT &pa, const T_PARAMETER& ta,
                                             const T_CONTROL_POINT &p0, const T_PARAMETER& t0,
                                             const T_CONTROL_POINT &p1, const T_PARAMETER& t1,
                                             const T_CONTROL_POINT &pb, const T_PARAMETER& tb,
                                             const T_PARAMETER& t) {
    T_PARAMETER t01 = t1 - t0;
    if (sl::is_zero(t01)) {
      return p0;
    } else {
#if 1
      return bezier_interpolation
        (p0,
         catmull_rom_to_bezier_control_point(pa, ta, p0, t0, t1), 
         catmull_rom_to_bezier_control_point(pb, tb, p1, t1, t0),
         p1,
         (t - t0)/t01);
#else
      T_PARAMETER t_prime = 
	bezier_interpolation
        (t0,
         lerp(t0, lerp(lerp(ta, t0, T_PARAMETER(2.0)), t1, T_PARAMETER(0.5)), T_PARAMETER(1.0/3.0)),
         lerp(t1, lerp(t0, lerp(tb, t1, T_PARAMETER(2.0)), T_PARAMETER(0.5)), T_PARAMETER(1.0/3.0)),
	 t1,
         (t - t0)/t01);
      return bezier_interpolation
        (p0,
         lerp(p0, lerp(lerp(pa, p0, T_PARAMETER(2.0)), p1, T_PARAMETER(0.5)), T_PARAMETER(1.0/3.0)),
         lerp(p1, lerp(p0, lerp(pb, p1, T_PARAMETER(2.0)), T_PARAMETER(0.5)), T_PARAMETER(1.0/3.0)),
	 p1,
         (t_prime - t0)/t01);
#endif
    }
  }

}; // namespace sl

namespace sl {

  /**
   *  Quick and dirty parameter track with interpolation.
   *  The track is represented by a map time to control point value.
   *  Standard std::map features are used to traverse and query
   *  the control point list. Feature value_at is used for
   *  interpolation.
   */
  template <class T_KNOT, class T_VALUE>
  class interpolation_track: public std::map<T_KNOT, T_VALUE> {
  public:
    typedef T_KNOT  knot_t;
    typedef T_VALUE value_t;
    typedef interpolation_track<knot_t, value_t> this_t;
    typedef std::map<knot_t, value_t> super_t;
    typedef typename super_t::key_type key_type;
    typedef typename super_t::value_type value_type;
    typedef typename super_t::key_compare key_compare;
    typedef typename super_t::value_compare value_compare;
    // not supported on win32:  typedef typename super_t::pointer pointer;
    typedef typename super_t::size_type size_type;
    typedef typename super_t::difference_type difference_type;
    typedef typename super_t::reference reference;
    typedef typename super_t::const_reference const_reference;
    typedef typename super_t::iterator iterator;
    typedef typename super_t::const_iterator const_iterator;
    typedef typename super_t::reverse_iterator reverse_iterator;
    typedef typename super_t::const_reverse_iterator const_reverse_iterator;
    
  protected:
    std::size_t continuity_;         // C_n
    bool        is_interpolating_;   // True: interpolating spling; False: approximating spline
  public: // Creation
    
    inline interpolation_track() : continuity_(2), is_interpolating_(true) {
    }

    inline ~interpolation_track() {
    }

    /// Set the continuity of the curve (0: position; 1: velocity; 2: acceleration...)
    inline void set_continuity(std::size_t Cn) {
      continuity_ = Cn;
    }

    /// Set whether you want to have more than positional continuity
    inline void set_smooth(bool b) {
      continuity_ = b ? 2 : 0;
    }
    
  public: // Queries

    std::size_t continuity() const {
      return continuity_;
    }
    
    inline bool is_smooth() const {
      return continuity_ > 0;
    }

    inline bool is_interpolating() const {
      return is_interpolating_ > 0;
    }
    
    inline knot_t start_time() const {
      SL_REQUIRE("Not empty", !this->empty());
      return this->begin()->first;
    }
    
    inline knot_t end_time() const {
      SL_REQUIRE("Not empty", !this->empty());
      return this->rbegin()->first;
    }
    
    value_t value_at(const knot_t& t) const {
      SL_REQUIRE("Not empty", !this->empty());

      if (t>= end_time()) {
	// Return last value
        return this->rbegin()->second;
      } else if (t<= start_time()) {
	// Return first value
        return this->begin()->second;
      }
      
      // Default case
      // Quick'n'dirty adaptation of Shoemake's DialaASpline. 
      // FIXME **** Very slow... *** Rewrite without copying
      // entire map to temporary storage!!
      std::vector<const_iterator> a; // Temporary, to emulate random access 
      a.reserve(this->size());
      for (const_iterator it=this->begin(); it!=this->end(); ++it) {
        a.push_back(it);
      }
      int n = a.size()-1;
      if (n==0) {
        return a[0]->second;
      }
      int Cn = std::min(n-1, (int)continuity());
      
      // Find enclosing knot interval
      int k=0;
      while (t>a[k]->first) { ++k; }
      int h=k;
      while (t==a[k]->first) { ++k; }     
      if (k>n) {k = n; if (h>k) h = k;}
      h = 1+Cn - (k-h); k--;
      int lo = k-Cn;
      int hi = k+1+Cn;

      std::vector<value_t> work;
      if (is_interpolating()) {
        // Lagrange interpolation steps
        int drop=0;
        if (lo<0) {
          lo = 0; drop += Cn-k;
          if (hi-lo<Cn) {
            drop += Cn-hi;
            hi = Cn;
          }
        }
        if (hi>n) {
          hi = n; drop += k+1+Cn-n;
          if (hi-lo<Cn) {
            drop += lo-(n-Cn);
            lo = n-Cn;
          }
        }
	if (hi-lo+1>0) work.reserve(hi-lo+1);
        for (int i=lo; i<=hi; i++) {
          work.push_back(a[i]->second);
        }
        for (int j=1; j<=Cn; j++) {
          for (int i=lo; i<=hi-j; i++) {
            work[i-lo] = linear_interpolation(work[i-lo], a[i]->first,
                                              work[i+1-lo], a[i+j]->first,
                                              t);
          }
        }
        h = 1+Cn-drop;
      } else {
        // Prepare for B-spline steps
        if (lo<0) {h += lo; lo = 0;}
        if (h+1>0) work.reserve(h+1);
	for (int i=lo; i<=lo+h; i++) {
          work.push_back(a[i]->second);
        }
        if (h<0) h = 0;
      }

      // B-spline steps
      for (int j=0; j<h; j++) {
        int tmp = 1+Cn-j;
        for (int i=h-1; i>=j; i--) {
          if (lo+i+tmp <=n) {
            work[i+1] = linear_interpolation(work[i], a[lo+i]->first,
                                             work[i+1], a[lo+i+tmp]->first,
                                             t);
          } else {
            work[i+1] = work[i];
          }
        }
      }

      return work[h];
    }
    
  };
}

#endif

