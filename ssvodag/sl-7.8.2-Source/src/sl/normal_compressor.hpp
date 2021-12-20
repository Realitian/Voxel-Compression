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
#ifndef SL_NORMAL_COMPRESSOR_HPP
#define SL_NORMAL_COMPRESSOR_HPP

#include <sl/integer.hpp>
#include <sl/math.hpp>
#include <sl/utility.hpp>
#include <cassert>

namespace sl {

  /**
   *  A compressor of normal vectors to the specified
   *  number of bits. 3 bits are used to identify the octant
   *  containing the vector, the other bits are used for
   *  the mantissa.
   *
   *  The code uses a radial projection from sphere to plane
   *  to obtain a nice sampling (see: R. Batista, Higher Accuracy
   *  quantized normals http://www.gamedev.net/reference/articles/article1252.asp)
   *  Quantization errors:
   *   16 bits: 1.91897 degrees
   *   32 bits: 0.00737 degrees
   */
  template< int G_Bits > 
  class normal_compressor {
  public:
    static const uint32_t Bits = G_Bits;
    typedef typename uint_t<Bits>::least compressed_t;  //< The type used for storing compressed normals
    static const uint32_t N = (8*sizeof(compressed_t)); //< The actual number of bits
    static const uint32_t N_octant     = 3;
    static const uint32_t N_mantissa   = (8*sizeof(compressed_t)-3);
    static const uint32_t N_mantissa_x = (8*sizeof(compressed_t)-3)/2;
    static const uint32_t N_mantissa_y = (8*sizeof(compressed_t)-3)-((8*sizeof(compressed_t)-3)/2);

  protected:
    
    const compressed_t Mask_octant_x;
    const compressed_t Mask_octant_y;
    const compressed_t Mask_octant_z;
    const compressed_t Mask_octant;
    const compressed_t Mask_mantissa;
    const compressed_t Mask_mantissa_x;
    const compressed_t Mask_mantissa_y;
    
    const double Scale_mantissa_x;
    const double Scale_mantissa_y;
    const double InvScale_mantissa_x;
    const double InvScale_mantissa_y;

  public:
    
    inline normal_compressor() :
        Mask_octant_x(compressed_t(1)<< N_mantissa),
        Mask_octant_y(compressed_t(1)<< (N_mantissa+1)),
        Mask_octant_z(compressed_t(1)<< (N_mantissa+2)),
        Mask_octant  ((compressed_t(1+2+4)) << N_mantissa),
        Mask_mantissa(compressed_t((~(compressed_t(1+2+4)) << N_mantissa))),
        Mask_mantissa_x((compressed_t(1)<< N_mantissa_x) -compressed_t(1)),
        Mask_mantissa_y((compressed_t(1)<< N_mantissa_y) -compressed_t(1)),
        Scale_mantissa_x(double((compressed_t(1)<< N_mantissa_x) -compressed_t(1))),
        Scale_mantissa_y(double((compressed_t(1)<< N_mantissa_y) -compressed_t(1))),
        InvScale_mantissa_x(1.0/double((compressed_t(1)<< N_mantissa_x) -compressed_t(1))),
	InvScale_mantissa_y(1.0/double((compressed_t(1)<< N_mantissa_y) -compressed_t(1))) {
#if 0
      std::cerr << "N_mantissa_x: " << N_mantissa_x << std::endl;
      std::cerr << "N_mantissa_y: " << N_mantissa_y << std::endl;
      std::cerr << "Mask_mantissa_x: " << (int64_t)Mask_mantissa_x << std::endl;
      std::cerr << "Mask_mantissa_y: " << (int64_t)Mask_mantissa_y << std::endl;
      std::cerr << "Scale_mantissa_x: " << Scale_mantissa_x << std::endl;
      std::cerr << "Scale_mantissa_y: " << Scale_mantissa_y << std::endl;
      std::cerr << "InvScale_mantissa_x: " << InvScale_mantissa_x << std::endl;
      std::cerr << "InvScale_mantissa_y: " << InvScale_mantissa_y << std::endl;
#endif
    }

    inline ~normal_compressor() {

    }

    /// Compress unit vector nx,ny,nz to cn  
    void compress_normal_in(compressed_t& cn,
                            double nx, double ny, double nz) const {
  
      cn = 0;
      
      // determine octant
      if (nx < 0.0) { cn |= Mask_octant_x; nx = -nx; }
      if (ny < 0.0) { cn |= Mask_octant_y; ny = -ny; }
      if (nz < 0.0) { cn |= Mask_octant_z; nz = -nz; }

      // project the normal onto the plane that goes through
      // X0=(1,0,0),Y0=(0,1,0),Z0=(0,0,1).
      // on that plane we choose an (projective!) coordinate system
      // such that X0->(0,0), Y0->(126,0), Z0->(0,126),(0,0,0)->Infinity
      assert( N_mantissa_x <= N_mantissa_y);
      double w = (Scale_mantissa_y-1.0) / ( nx + ny + nz );
      compressed_t cx = (compressed_t)( nx * w );
      compressed_t cy = (compressed_t)( ny * w );
      assert( cx <  Mask_mantissa_y );
      assert( cy <  Mask_mantissa_y );
      if (cx > Mask_mantissa_x) { 
        cx = Mask_mantissa_y - cx;
        cy = Mask_mantissa_y - cy;
      }
      assert( cx <= Mask_mantissa_x );
      assert( cy <= Mask_mantissa_y );

      // determine mantissa
      compressed_t cxy = cx | (cy << N_mantissa_x);
      // combine octant and mantissa
      cn |= cxy;
    }
    
    /// Decompress unit vector nx,ny,nz from cn  
    void uncompress_normal_in(double& nx, double& ny, double& nz,
                              compressed_t cn) const {
      // extract mantissa
      compressed_t cx = (cn & Mask_mantissa_x);
      compressed_t cy = (cn >> N_mantissa_x) & Mask_mantissa_y;
      assert( N_mantissa_x <= N_mantissa_y);
      if (cx+cy >= Mask_mantissa_y) {
        cx = Mask_mantissa_y - cx;
        cy = Mask_mantissa_y - cy;
      }
      assert( cx <  Mask_mantissa_y );
      assert( cy <  Mask_mantissa_y );
      
      // build approximate normal
      nx = double(cx);
      ny = double(cy);
      nz = (Scale_mantissa_y-1.0) - nx - ny;
      double len = std::sqrt(nx*nx+ny*ny+nz*nz);
      nx /= len;
      ny /= len;
      nz /= len;
      
      // determine octant
      if (cn & Mask_octant_x) { nx = -nx; }
      if (cn & Mask_octant_y) { ny = -ny; }
      if (cn & Mask_octant_z) { nz = -nz; }
    }
    
    /// Decompress unit vector nx,ny,nz from cn  
   void uncompress_normal_in(float& nx, float& ny, float& nz,
			     compressed_t cn) const {
     // extract mantissa
     compressed_t cx = (cn & Mask_mantissa_x);
     compressed_t cy = (cn >> N_mantissa_x) & Mask_mantissa_y;
     assert( N_mantissa_x <= N_mantissa_y);
     if (cx+cy >= Mask_mantissa_y) {
       cx = Mask_mantissa_y - cx;
       cy = Mask_mantissa_y - cy;
     }
     assert( cx <  Mask_mantissa_y );
     assert( cy <  Mask_mantissa_y );
      
     // build approximate normal
     nx = float(cx);
     ny = float(cy);
     nz = (Scale_mantissa_y-1.0f) - nx - ny;
     float len = std::sqrt(nx*nx+ny*ny+nz*nz);
     nx /= len;
     ny /= len;
     nz /= len;
     
     // determine octant
     if (cn & Mask_octant_x) { nx = -nx; }
     if (cn & Mask_octant_y) { ny = -ny; }
     if (cn & Mask_octant_z) { nz = -nz; }
    }
    
    /// Decompress unit vector nx,ny,nz from cn  
    void uncompress_normal_in(int32_t& nx, int32_t& ny, int32_t& nz,
                              compressed_t cn) const {
      double dnx, dny, dnz;
      uncompress_normal_in(dnx,dny,dnz,cn);
      nx = int32_t(dnx*2147483647.0);
      ny = int32_t(dny*2147483647.0);
      nz = int32_t(dnz*2147483647.0);
    }
    
    /// Decompress unit vector nx,ny,nz from cn  
    void uncompress_normal_in(int16_t& nx, int16_t& ny, int16_t& nz,
                              compressed_t cn) const {
      double dnx, dny, dnz;
      uncompress_normal_in(dnx,dny,dnz,cn);
      nx = int16_t(dnx*32767.0);
      ny = int16_t(dny*32767.0);
      nz = int16_t(dnz*32767.0);
    }

    /// Decompress unit vector nx,ny,nz from cn  
    void uncompress_normal_in(int8_t& nx, int8_t& ny, int8_t& nz,
			     compressed_t cn) const {
      double dnx, dny, dnz;
      uncompress_normal_in(dnx,dny,dnz,cn);
      nx = int8_t(dnx*127.0);
      ny = int8_t(dny*127.0);
      nz = int8_t(dnz*127.0);
    }
    
    /// An estimate of compressed representation accuracy (in radians)
    double accuracy() const {
      const int iS = 1024;
      double dDotMin = 1.0;

      for (int sX = -1; sX <= 1; sX+=2) {
        for (int sY = -1; sY <= 1; sY+=2) {
          for (int sZ = -1; sZ <= 1; sZ+=2) {
            for (int iY = 0; iY < iS; iY++) {
              double ny0 = sY * iY/(double)iS;
              for (int iX = 0; iX < iS; iX++) {
                double nx0 = sX * iX/(double)iS;
                double nz0 = 1.0 - nx0*nx0 - ny0*ny0;
                if ( nz0 >= 0.0 ) {
                  nz0 = sZ * std::sqrt(nz0);
                  compressed_t cn;
                  compress_normal_in(cn,nx0,ny0,nz0);             
                  double nx1, ny1, nz1;
                  uncompress_normal_in(nx1,ny1,nz1, cn);
                  
                  double dDot = nx0*nx1+ny0*ny1+nz0*nz1;
                  if ( dDot < dDotMin ) {
                    dDotMin = dDot;
                  }
                }
              }
            }
          }
        }
      }

      return std::acos(dDotMin);
    }
    
  };

}
    
#endif
