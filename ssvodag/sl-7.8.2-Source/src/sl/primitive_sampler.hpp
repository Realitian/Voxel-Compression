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
#ifndef SL_PRIMITIVE_SAMPLER_HPP
#define SL_PRIMITIVE_SAMPLER_HPP

#include <sl/fixed_size_point.hpp>

namespace sl {
    
  /**
   *  Uniform stratified sampling over the
   *  unit square
   */
  template <class T>
  class square_sampler {
  public:
    typedef T                     value_t; 
    typedef square_sampler<T>     self_t;
    typedef fixed_size_point<2,T> point_t;
  protected:
    enum { N = 1024U };
    static value_t samples_[N*2];
    size_t k_;
  public:

    /// Create a new generator initialized with the given seed
    inline square_sampler(size_t seed = 0)
      : k_(seed%N) {
    }
    
    /// The current sampling point (2D coordinates within the unit square)
    inline point_t current() const {
      return point_t(samples_[k_*2+0], samples_[k_*2+1]);
    }

    /// Advance to the next sampling point
    inline void forth() {
      k_ = (k_ + 1) % N;
    }
  };

  /**
   *  Uniform stratified sampling over the
   *  unit triangle
   */
  template <class T>
  class triangle_sampler {
  public:
    typedef T                     value_t; 
    typedef triangle_sampler<T>     self_t;
    typedef fixed_size_point<2,T> point_t;
  protected:
    enum { N = 1024U };
    static value_t samples_[N*2];
    size_t k_;
  public:

    /// Create a new generator initialized with the given seed
    inline triangle_sampler(size_t seed = 0)
      : k_(seed%N) {
    }
    
    /// The current sampling point (2D barycentric coordinates within the unit square
    inline point_t current() const {
      return point_t(samples_[k_*2+0], samples_[k_*2+1]);
    }

    /// Advance to the next sampling point
    inline void forth() {
      k_ = (k_ + 1) % N;
    }
  };

  /**
   *  Uniform stratified sampling over the
   *  unit circle
   */
  template <class T>
  class circle_sampler {
  public:
    typedef T                     value_t; 
    typedef circle_sampler<T>     self_t;
    typedef fixed_size_point<2,T> point_t;
  protected:
    enum { N = 1024U };
    static value_t samples_[N*2];
    size_t k_;
  public:

    /// Create a new generator initialized with the given seed
    inline circle_sampler(size_t seed = 0)
      : k_(seed%N) {
    }
    
    /// The current sampling point (2D cartesian coordinates within the circle)
    inline point_t current() const {
      return point_t(samples_[k_*2+0], samples_[k_*2+1]);
    }

    /// Advance to the next sampling point
    inline void forth() {
      k_ = (k_ + 1) % N;
    }
  };

  /**
   *  Uniform stratified sampling over the
   *  unit cube
   */
  template <class T>
  class cube_sampler {
  public:
    typedef T                     value_t; 
    typedef cube_sampler<T>     self_t;
    typedef fixed_size_point<3,T> point_t;
  protected:
    enum { N = 1024U };
    static value_t samples_[N*3];
    size_t k_;
  public:

    /// Create a new generator initialized with the given seed
    inline cube_sampler(size_t seed = 0)
      : k_(seed%N) {
    }
    
    /// The current sampling point (3D cartesian coordinates)
    inline point_t current() const {
      return point_t(samples_[k_*3+0], samples_[k_*3+1], samples_[k_*3+2]);
    }

    /// Advance to the next sampling point
    inline void forth() {
      k_ = (k_ + 1) % N;
    }
  };

  /**
   *  Uniform stratified sampling over the
   *  unit sphere
   */
  template <class T>
  class sphere_sampler {
  public:
    typedef T                     value_t; 
    typedef sphere_sampler<T>     self_t;
    typedef fixed_size_point<3,T> point_t;
  protected:
    enum { N = 1024U };
    static value_t samples_[N*3];
    size_t k_;
  public:

    /// Create a new generator initialized with the given seed
    inline sphere_sampler(size_t seed = 0)
      : k_(seed%N) {
    }
    
    /// The current sampling point (3D cartesian coordinates)
    inline point_t current() const {
      return point_t(samples_[k_*3+0], samples_[k_*3+1], samples_[k_*3+2]);
    }

    /// Advance to the next sampling point
    inline void forth() {
      k_ = (k_ + 1) % N;
    }
  };

}

#endif
