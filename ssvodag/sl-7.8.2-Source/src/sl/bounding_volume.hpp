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
#ifndef SL_BOUNDING_VOLUME_HPP
#define SL_BOUNDING_VOLUME_HPP

#include <sl/utility.hpp> // std::pair, SL_DECLARE_GENERIC_SUPERCLASS_*, ...
#include <sl/fixed_size_ray.hpp>
#include <sl/fixed_size_line.hpp>
#include <sl/fixed_size_line_segment.hpp>

namespace sl {

  /// Generic base class for bounding volumes in N dimensions
  template <size_t DIMENSION, class T, class T_DERIVED, class T_DATA> 
  class bounding_volume_base {
  public:
    typedef bounding_volume_base<DIMENSION, T, T_DERIVED, T_DATA> this_t;
    typedef T_DERIVED                                             derived_t;

    enum { dimension = DIMENSION };
    typedef T_DATA                                                data_t;
    typedef T                                                     value_t; 
    typedef fixed_size_point<DIMENSION,T>                         point_t;
    typedef fixed_size_line<DIMENSION,T>                          line_t;
    typedef fixed_size_ray<DIMENSION,T>                           ray_t;
    typedef fixed_size_line_segment<DIMENSION,T>                  line_segment_t;

    /// Features to access this as a derived_t pointer
    SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  protected:  // Data

    T_DATA data_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << data_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> data_;
    }

  public: // Creation & Destruction

    /// default init
    inline bounding_volume_base() {}

    /// data init
    inline bounding_volume_base(const data_t& data): data_(data) {}

    /// default destructor
    inline ~bounding_volume_base() {}

  public: // Queries
				
    /// the volume of the object
    inline value_t volume() const {
      return derived_ref().volume();
    }

  public: // Containment

    /// does this contain point p?
    inline bool contains(const point_t& p) const {
      return derived_ref().contains(p);
    }

  public: // Ray-tracing

    /// True iff the ray intersects the volume
    inline bool intersection_exists(const ray_t& r) const {
      return derived_ref().intersection_exists(r);
    }

    /// True iff the line intersects the volume
    inline bool intersection_exists(const line_t& r) const {
      return derived_ref().intersection_exists(r);
    }

    /// True iff the line segment intersects the volume
    inline bool intersection_exists(const line_segment_t& r) const {
      return derived_ref().intersection_exists(r);
    }

    /**
     * ray - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest ray parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const ray_t& r) const {
      return derived_ref().intersection(r);
    }

    /**
     * line - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest line parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const line_t& r) const {
      return derived_ref().intersection(r);
    }

    /**
     * segment - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const line_segment_t& r) const {
      return derived_ref().intersection(r);
    }

    /**
     * ray - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest ray parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const ray_t& r) const {
      if (contains(r.origin())) {
	return std::pair<bool,value_t>(true,scalar_math<value_t>::zero());
      } else {
	return intersection(r);
      }
    }

    /**
     * line - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest line parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const line_t& r) const {
      return intersection(r);
    }

    /**
     * segment - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const line_segment_t& r) const {
      if (contains(r.origin())) {
	return std::pair<bool,value_t>(true,scalar_math<value_t>::zero());
      } else if (contains(r.extremity())) {
	return std::pair<bool,value_t>(true,scalar_math<value_t>::one());
      } else {
	return intersection(r);
      }
    }

  }; // bounding volume

}; // namespace sl

// --------------------------------------------------------------
// sl::bounding_volume_builder<DIMENSION,T>
// --------------------------------------------------------------

namespace sl {

  /// Abstract objects for constructing bounding volumes from point clouds
  template <class T_BOUNDING_VOLUME> 
  class bounding_volume_builder {
  public:

    typedef bounding_volume_builder<T_BOUNDING_VOLUME> this_t;
    typedef T_BOUNDING_VOLUME                          bounding_volume_t;

    enum { dimension = bounding_volume_t::dimension };
    typedef typename bounding_volume_t::value_t        value_t;
    typedef typename bounding_volume_t::point_t        point_t;

  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  protected: // State

    bool is_building_;

    bounding_volume_t last_bounding_volume_;

  public: // Creation & Destruction

    /// default constructor
    inline bounding_volume_builder() { is_building_ = false; }

    /// default destructor
    virtual ~bounding_volume_builder() {}

  public: // Construction

    /// is the object currently building a bounding volume?
    inline bool is_building() const {
      return is_building_;
    }

    /// start building a volume
    virtual inline void begin_model() {
      SL_REQUIRE("Not is building", !is_building());
      is_building_ = true;
      SL_ENSURE("Is building", is_building());
    }

    /// add point to current point cloud
    virtual void put_point(const point_t& p) = 0;

    /// end point cloud and build bounding volume
    virtual inline void end_model() {
      SL_REQUIRE("Is building", is_building());
      is_building_ = false;
      SL_ENSURE("Not is building", !is_building());
    }
      
    /// the bounding volume constructed by last end_model
    inline const bounding_volume_t& last_bounding_volume() const {
      return last_bounding_volume_;
    }

  }; // class bounding_volume_builder

}; // namespace sl

#endif
