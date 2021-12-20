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
#ifndef SL_CIE_HPP
#define SL_CIE_HPP

#include <sl/config.hpp>

namespace sl {

  /// Color conversion based on the CIE model
  template <class T_value>
  class cie {
  public:
    typedef T_value value_t;
  protected:
    value_t x_r_, y_r_, 
          x_g_, y_g_, 
          x_b_, y_b_,
          x_w_, y_w_;
    value_t rf_, gf_, bf_;

    value_t xyz2rgbmat_[3][3];
    value_t rgb2xyzmat_[3][3];

    value_t l_;

  public:

    cie(); 

    /// Set reference primaries for XYZ->RGB conversion
    void set_reference_primaries(value_t xr, value_t yr,
				 value_t xg, value_t yg,
				 value_t xb, value_t yb,
				 value_t xw, value_t yw);

    /// Set the reference primaries for a standard NTSC monitor
    void set_ntsc_standard_monitor();

    /// Set the reference primaries for a standard PAL monitor
    void set_pal_standard_monitor();

    /// Set the value used for tristimulus white efficacy [lm/W]
    inline void set_luminous_efficacy(value_t x) {
      l_ = x;
    }

    /** 
     *  The Luminous efficacy [lm/W] for conversion between radiometric 
     *  and photometric units (e.g., for conversion between radiance and
     *  luminance).
     */
    inline value_t luminous_efficacy() const {
      return l_;
    }
    
  public: // Parameters

    inline value_t x_r() const { return x_r_; }
    inline value_t x_g() const { return x_g_; }
    inline value_t x_b() const { return x_b_; }
    inline value_t x_w() const { return x_w_; }

    inline value_t y_r() const { return y_r_; }
    inline value_t y_g() const { return y_g_; }
    inline value_t y_b() const { return y_b_; }
    inline value_t y_w() const { return y_w_; }

    inline value_t D() const { return 
			       x_r()*(y_g() - y_b()) + 
			       x_g()*(y_b() - y_r()) + 
			       x_b()*(y_r() - y_g()); }

    inline value_t C_rD() const { return 
				  (((value_t)1.0)/y_w()) * 
				  ( x_w()*(y_g() - y_b()) -
				    y_w()*(x_g() - x_b()) + 
				    x_g()*y_b() - x_b()*y_g()); }
    inline value_t C_gD() const { return 
				  (((value_t)1.0)/y_w()) *
				  ( x_w()*(y_b() - y_r()) -
				    y_w()*(x_b() - x_r()) -
				    x_r()*y_b() + x_b()*y_r() ); }
    inline value_t C_bD() const { return
				  (((value_t)1.0)/y_w()) * 
				  ( x_w()*(y_r() - y_g()) - 
				    y_w()*(x_r() - x_g()) + 
				    x_r()*y_g() - x_g()*y_r() ); }

    inline value_t rf() const { return rf_; }
    inline value_t gf() const { return gf_; }
    inline value_t bf() const { return bf_; }

  public: // Conversion

    /// Convert rgb to xyz using the current primary definition
    inline void rgb_to_xyz(value_t  r, value_t  g, value_t  b,
			   value_t& x, value_t& y, value_t& z) const {
      x = rgb2xyzmat_[0][0] * r + rgb2xyzmat_[0][1] * g + rgb2xyzmat_[0][2] * b;
      y = rgb2xyzmat_[1][0] * r + rgb2xyzmat_[1][1] * g + rgb2xyzmat_[1][2] * b;
      z = rgb2xyzmat_[2][0] * r + rgb2xyzmat_[2][1] * g + rgb2xyzmat_[2][2] * b;
    }

    /// Convert cie to rgb using the current primary definition
    inline void xyz_to_rgb(value_t  x, value_t  y, value_t   z,
			   value_t& r, value_t& g, value_t& b) const {
      r = xyz2rgbmat_[0][0] * x + xyz2rgbmat_[0][1] * y + xyz2rgbmat_[0][2] * z;
      g = xyz2rgbmat_[1][0] * x + xyz2rgbmat_[1][1] * y + xyz2rgbmat_[1][2] * z;
      b = xyz2rgbmat_[2][0] * x + xyz2rgbmat_[2][1] * y + xyz2rgbmat_[2][2] * z;    
    }

    /// An achromatic value representing the spectral quantity [W]
    inline value_t achromatic_from_rgb(value_t r, value_t g, value_t b) const {
      return rf() * r + gf() * g + bf() * b;
    }

    /// The luminance [lm], a photometric quantity corresponding to the radiance of the given spectrum.
    inline value_t luminance_from_rgb(value_t r, value_t g, value_t b) const {
      return luminous_efficacy() * achromatic_from_rgb(r,g,b);
    }

  }; 

  // -- Inline implementation

  template<class T_value>
  cie<T_value>::cie() {
    //    const value_t MAX_EFFICACY   = 683.002;  // photopic maximum for 555 nm
    const value_t WHITE_EFFICACY = 183.07;   // CIE 1988 uniform equal energy 1W white light

    l_ = WHITE_EFFICACY;

    set_pal_standard_monitor();
  }

  template<class T_value>
  void cie<T_value>::set_ntsc_standard_monitor() {
    // NTSC standard NTSC primaries
    static const value_t  CIE_NTSC_x_r =                0.670; 
    static const value_t  CIE_NTSC_y_r =                0.330;
    static const value_t  CIE_NTSC_x_g =                0.210;
    static const value_t  CIE_NTSC_y_g =                0.710;
    static const value_t  CIE_NTSC_x_b =                0.140;
    static const value_t  CIE_NTSC_y_b =                0.080;
 
    /// NTSC standard monitor white point
    static const value_t  CIE_NTSC_x_w =                0.3333333333;
    static const value_t  CIE_NTSC_y_w =                0.3333333333;

    set_reference_primaries(CIE_NTSC_x_r, CIE_NTSC_y_r,
			  CIE_NTSC_x_g, CIE_NTSC_y_g,
			  CIE_NTSC_x_b, CIE_NTSC_y_b,
			  CIE_NTSC_x_w, CIE_NTSC_y_w);		  
  }

  template<class T_value>
  void cie<T_value>::set_pal_standard_monitor() {
    // PAL nominal CRT primaries
    static const value_t  CIE_PAL_x_r =                0.640;
    static const value_t  CIE_PAL_y_r =                0.330;
    static const value_t  CIE_PAL_x_g =                0.290;
    static const value_t  CIE_PAL_y_g =                0.600;
    static const value_t  CIE_PAL_x_b =                0.150;
    static const value_t  CIE_PAL_y_b =                0.060;
     
    /// PAL standard monitor white point
    static const value_t  CIE_PAL_x_w =                0.3333333333; 
    static const value_t  CIE_PAL_y_w =                0.3333333333;

    set_reference_primaries(CIE_PAL_x_r, CIE_PAL_y_r,
			  CIE_PAL_x_g, CIE_PAL_y_g,
			  CIE_PAL_x_b, CIE_PAL_y_b,
			  CIE_PAL_x_w, CIE_PAL_y_w);		  
  }

  template<class T_value>
  void cie<T_value>::set_reference_primaries(value_t xr, value_t yr,
					     value_t xg, value_t yg,
					     value_t xb, value_t yb,
					     value_t xw, value_t yw) {
    x_r_ = xr;  y_r_ = yr;
    x_g_ = xg;  y_g_ = yg;
    x_b_ = xb;  y_b_ = yb;
    x_w_ = xw;  y_w_ = yw;

    rf_ = (y_r()*C_rD()/D());
    gf_ = (y_g()*C_gD()/D());
    bf_ = (y_b()*C_bD()/D());
      
    xyz2rgbmat_[0][0] = (y_g() - y_b() - x_b()*y_g() + y_b()*x_g())/C_rD();
    xyz2rgbmat_[0][1] = (x_b() - x_g() - x_b()*y_g() + x_g()*y_b())/C_rD();
    xyz2rgbmat_[0][2] = (x_g()*y_b() - x_b()*y_g())/C_rD();

    xyz2rgbmat_[1][0] = (y_b() - y_r() - y_b()*x_r() + y_r()*x_b())/C_gD();
    xyz2rgbmat_[1][1] = (x_r() - x_b() - x_r()*y_b() + x_b()*y_r())/C_gD();
    xyz2rgbmat_[1][2] = (x_b()*y_r() - x_r()*y_b())/C_gD();

    xyz2rgbmat_[2][0] = (y_r() - y_g() - y_r()*x_g() + y_g()*x_r())/C_bD();
    xyz2rgbmat_[2][1] = (x_g() - x_r() - x_g()*y_r() + x_r()*y_g())/C_bD();
    xyz2rgbmat_[2][2] = (x_r()*y_g() - x_g()*y_r())/C_bD();

    rgb2xyzmat_[0][0] = x_r()*C_rD()/D();
    rgb2xyzmat_[0][1] = x_g()*C_gD()/D();
    rgb2xyzmat_[0][2] = x_b()*C_bD()/D();

    rgb2xyzmat_[1][0] = y_r()*C_rD()/D();
    rgb2xyzmat_[1][1] = y_g()*C_gD()/D();
    rgb2xyzmat_[1][2] = y_b()*C_bD()/D();

    rgb2xyzmat_[2][0] = (((value_t)1.0)-x_r()-y_r())*C_rD()/D();
    rgb2xyzmat_[2][1] = (((value_t)1.0)-x_g()-y_g())*C_gD()/D();
    rgb2xyzmat_[2][2] = (((value_t)1.0)-x_b()-y_b())*C_bD()/D();
  }

}

#endif
