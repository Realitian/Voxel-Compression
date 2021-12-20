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
#ifndef SL_SHAFT_HPP
#define SL_SHAFT_HPP

#include <sl/axis_aligned_box.hpp>
#include <sl/oriented_box.hpp>
#include <vector>

namespace sl {

  /**
   *  Shaft, i.e. the volume between two axis-aligned bounding boxes.
   *  Implementation based on "Shaft Culling for Efficient Ray-Traced Radiosity,"
   *  Eric A. Haines and John R. Wallace, Photorealistic Rendering in
   *  Computer Graphics (Proceedings of the Second Eurographics Workshop
   *  on Rendering), Springer-Verlag, New York, 1994, p.122-138.
   *
   *  To test a given bounding volume against the shaft, a typical sequence is:
   *
   *      if (!shaft.culls_out(test_volume)) {
   *          if (shaft.contains(test_volume)) {
   *              ... test_volume is fully inside shaft ...
   *          } else {
   *              ... test_volume overlaps shaft, but not fully inside ...
   *          }
   *      } else {
   *          ... test_volume is fully outside shaft ...
   *      }
   */
  template <class T>
  class shaft {
  public:
    enum { dimension = 3 };

    typedef T                     value_t; 
    typedef shaft<T>              self_t;

    typedef axis_aligned_box<3,T>                        box_t;
    typedef oriented_box<3,T>                            oriented_box_t;

    typedef fixed_size_plane<3,T>                        plane_t;
    typedef typename              plane_t::point_t       point_t;
    typedef typename              plane_t::vector_t      vector_t;
    typedef typename              plane_t::dual_vector_t dual_vector_t;

    typedef typename std::vector<plane_t>::iterator       plane_iter_t;
    typedef typename std::vector<plane_t>::const_iterator const_plane_iter_t;
 
  protected:
    box_t  box_; 
    std::vector<plane_t> planes_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << box_ << planes_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> box_ >> planes_;
    }
  public:

    /// Build a shaft from default boxes 
    shaft();

    /// Build a shaft from two axis aligned boxes
    shaft(const box_t& b1,
	  const box_t& b2);

    /// Build a shaft from two oriented boxes
    shaft(const oriented_box_t& b1,
	  const oriented_box_t& b2);

    /// Add an extra plane to the shaft
    void add_plane(const plane_t& hp);

    /// Does the shaft contain p?
    bool contains(const point_t& p) const;

    /// Is the point p outside the shaft?
    bool culls_out(const point_t& p) const;

    /// Is the box fully outside the shaft?
    bool culls_out(const box_t& b) const;

    /// Is the box fully inside the shaft?
    bool contains(const box_t& b) const;

    /// Is the sphere fully outside the shaft?
    bool culls_out(const point_t& sph_center,
		   value_t sph_radius) const;

    /// Is the sphere fully inside the shaft?
    bool contains(const point_t& sph_center,
		  value_t sph_radius) const;

   /// Is the box fully outside the shaft?
    bool culls_out(const oriented_box_t& b) const;

    /// Is the box fully inside the shaft?
    bool contains(const oriented_box_t& b) const;
  };

}

//---- Shaft implementation

namespace sl {

  template <class T>
  inline shaft<T>::shaft() {
  }

  template <class T>
  inline void shaft<T>::add_plane(const plane_t& hp) {
    planes_.push_back(hp);
  }

  template <class T>
  shaft<T>::shaft(const box_t& box0,
		  const box_t& box1) {
    const int LO_X = 0;
    const int HI_Z = 5;

    int match[2][dimension];
    int face_tally[2][dimension];
    
    /* Store union of the two bounding boxes in the shaft's box structure.
     * Also set up "match", which tells whether two box coordinates are
     * the same value, and "face_tally", which says whether box box0 or box1 is
     * the one which defined the shaft's box (i.e. was the farther one out
     * in the given X/Y/Z -/+ coordinate direction).
     */
    for (int dim=0; dim<dimension; ++dim) {
      // Lower bounds
      match[0][dim] = (box0[0][dim] == box1[0][dim]);
      if (box0[0][dim] < box1[0][dim]) {
	match[0][dim] = false;
	box_[0][dim] = box0[0][dim];
	face_tally[0][dim] = 0;
      } else {
	match[0][dim] = (box0[0][dim] == box1[0][dim]);
	box_[0][dim] = box1[0][dim];
	face_tally[0][dim] = 1;
      }
      // Upper bounds
      if (box0[1][dim] > box1[1][dim]) {
	match[1][dim] = false;
	box_[1][dim] = box0[1][dim];
	face_tally[1][dim] = 0;
      } else {
	match[1][dim] = (box0[1][dim] == box1[1][dim]);
	box_[1][dim] = box1[1][dim];
	face_tally[1][dim] = 1;
      }
    }

    /* Search through all adjacent cube faces and see if they
     * should be joined by a plane. If faceTally differs, then
     * a plane (could be) added.
     */
    dual_vector_t pn = tags::not_initialized();
    for (int i1=LO_X; i1 <= HI_Z-1 /* HI_Z tested below */ ; ++i1) {
      int i1n = i1%3; int i1l = i1/3;
      if (!match[i1l][i1n]) {
	// loop through cube faces above current face */
	for (int i2 = i1+1; i2 <= HI_Z; ++i2) {
	  int i2n = i2 % 3;  int i2l = i2/3;
	  if ((!match[i2l][i2n]) && 
	      (i1+3 != i2 ) &&  // Do not bother with opposite faces
	      (face_tally[i1l][i1n] != face_tally[i2l][i2n])) {
		// A real split exists 

	    int i3n= ((i1n+1)%3 == i2n )? (i2n+1)%3 : (i1n+1)%3;

	    /* Plane's normal: made by rotating the plane
	     * joining the two faces 90 degrees clockwise.
	     */
	    pn[i1n] = box0[i2l][i2n] - box1[i2l][i2n];
	    
	    /* if the face is negative and the normal component points positive,
	     * or face/normal is positive/negative, then the normal must be flipped.
	     */
	    if ( ( i1l == 0 ) != (sl::is_negative(pn[i1n])) ) {
	      pn[i1n] = -pn[i1n] ;
	      pn[i2n] = box0[i1l][i1n] - box1[i1l][i1n];
	    } else {
	      pn[i2n] = box1[i1l][i1n] - box0[i1l][i1n];
	    }
	    pn[i3n] = sl::zero(value_t());
	    
	    // for sphere testing, the normal must be normalized!
	    value_t len = std::sqrt(sl::sqr(pn[i1n]) + 
				    sl::sqr(pn[i2n]));
	    pn[i1n] /= len ;
	    pn[i2n] /= len ;
	    
	    // compute d, offset of plane from origin and add plane to shaft
	    add_plane(plane_t(pn, -(box0[i1l][i1n] * pn[i1n] + box0[i2l][i2n] * pn[i2n])));
	  }
	}
      }
    }
  }

  template <class T>
  shaft<T>::shaft(const oriented_box_t& box0,
		  const oriented_box_t& box1) {
    // Sloppy implementation: the oriented boxes
    // are converted to axis aligned boxes before
    // constructing the shaft!
    self_t other(box0.bounding_axis_aligned_box(),
		 box1.bounding_axis_aligned_box());
    box_ = other.box_;
    planes_ = other.planes_;
  }

  template <class T>
  bool shaft<T>::contains(const point_t& p) const {
    bool result = box_.contains(p);

    for (const_plane_iter_t it = planes_.begin();
	 !(!result || it == planes_.end());
	 ++it) {
      result = sl::is_non_positive((*it).value(p));
    }

    return result;
  }

  template <class T>
  bool shaft<T>::culls_out(const point_t& p) const {
    return !contains(p);
  }
    
  template <class T>
  bool shaft<T>::culls_out(const box_t& b) const {
    bool result = 
      (b[0][0] > box_[1][0] ||
       b[0][1] > box_[1][1] ||
       b[0][2] > box_[1][2] ||
       b[1][0] < box_[0][0] ||
       b[1][1] < box_[0][1] ||
       b[1][2] < box_[0][2]);

    for (const_plane_iter_t it = planes_.begin();
	 !(result || it == planes_.end());
	 ++it) {
      
      result = b.is_fully_above((*it));
    }

    return result;
  }

  template <class T>
  bool shaft<T>::contains(const box_t& b) const {
    bool result = 
      !(b[0][0] < box_[0][0] ||
	b[0][1] < box_[0][1] ||
	b[0][2] < box_[0][2] ||
	b[1][0] > box_[1][0] ||
	b[1][1] > box_[1][1] ||
	b[1][2] > box_[1][2]);

    for (const_plane_iter_t it = planes_.begin();
	 !(!result || it == planes_.end());
	 ++it) {
      
      result = b.is_fully_below((*it));
    }

    return result;
  }

  template <class T>
  bool shaft<T>::culls_out(const point_t& sph_center,
			   value_t sph_radius) const {
    bool result = 
      ( sph_center[0] - sph_radius > box_[1][0] ||
	sph_center[1] - sph_radius > box_[1][1] ||
	sph_center[2] - sph_radius > box_[1][2] ||
	sph_center[0] + sph_radius < box_[0][0] ||
	sph_center[1] + sph_radius < box_[0][1] ||
	sph_center[2] + sph_radius < box_[0][2] );

    for (const_plane_iter_t it = planes_.begin();
	 !(result || it == planes_.end());
	 ++it) {
      
      result = (*it).value(sph_center) > sph_radius;
    }
    return result;
  }

  template <class T>
  bool shaft<T>::contains(const point_t& sph_center,
			  value_t sph_radius) const {
    bool result = 
      !( sph_center[0] - sph_radius > box_[0][0] ||
	 sph_center[1] - sph_radius > box_[0][1] ||
	 sph_center[2] - sph_radius > box_[0][2] ||
	 sph_center[0] + sph_radius < box_[1][0] ||
	 sph_center[1] + sph_radius < box_[1][1] ||
	 sph_center[2] + sph_radius < box_[1][2] );

    for (const_plane_iter_t it = planes_.begin();
	 !(!result || it == planes_.end());
	 ++it) {      
      result = (*it).value(sph_center) > -sph_radius;
    }
    return result;
  }


  template <class T>
  bool shaft<T>::culls_out(const oriented_box_t& b) const {
    bool result = false;
    for (size_t d=0; !(d==dimension || result); ++d) {
      interval<value_t> d_range = b.signed_distance_range(d);
      result = d_range[1] < box_[0][d] || d_range[0] > box_[1][d];
    }
    for (const_plane_iter_t it = planes_.begin();
	 !(result || it == planes_.end());
	 ++it) {
      
      result = b.is_fully_above((*it));
    }
    return result;
  }

  template <class T>
  bool shaft<T>::contains(const oriented_box_t& b) const {
    bool result = true;
    for (size_t d=0; !(d==dimension || !result); ++d) {
      interval<value_t> d_range = b.signed_distance_range(d);
      result = box_[0][d] <= d_range[0] && d_range[1] <= box_[1][d];
    }
    for (const_plane_iter_t it = planes_.begin();
	 !(!result || it == planes_.end());
	 ++it) {
      
      result = b.is_fully_below((*it));
    }
    return result;
  }

}

#endif
