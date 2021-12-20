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
#ifndef SL_UNIFORM_GRID_CONTAINER_HPP
#define SL_UNIFORM_GRID_CONTAINER_HPP

#include <sl/index.hpp>
#include <sl/axis_aligned_box.hpp>
#include <list>
#include <cassert>

namespace sl {

  /**
   * A spatial search structure for a accessing a container of
   * objects. It is based on a uniform grid overlayed over a protion of
   * space.  The grid partion the space into cells. Cells contains just
   * pointers to the object that are stored elsewhere.  The set of object
   * is meant to be static and pointer stable.  Useful for situation were
   * many space related query are issued over the same dataset (ray
   * tracing, measuring distances between meshes, re-detailing ecc.). Works
   * well for distributions that are reasonably uniform. 
   */
  template <size_t DIMENSION, class Scalar_T, class Object_T>
  class uniform_grid_container {
  public:
    enum {  dimension = DIMENSION };

    typedef uniform_grid_container<DIMENSION,Scalar_T, Object_T> this_t;
    typedef Object_T                                             object_t;
    typedef Scalar_T                                             value_t;
    
    typedef index<dimension>                                     subscript_t;
    
    typedef fixed_size_point<dimension,value_t>                         point_t;
    typedef fixed_size_vector<sl::column_orientation,dimension,value_t> vector_t;
    typedef axis_aligned_box<dimension,value_t>                         aabox_t;
    typedef fixed_size_ray<dimension,value_t>                           ray_t;

    typedef std::list<object_t>                                         object_list_t;
    typedef std::vector<object_list_t>                                  data_t;
    
  protected:

    subscript_t extent_;
    subscript_t stride_;
    data_t      data_;
    aabox_t     bbox_;
    vector_t    voxel_dim_;
    
  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << extent_ << stride_ << data_;
    }
    
    void retrieve_from(input_serializer& s) {
      s >> extent_ >> stride_ >> data_;
    }

  protected: // Helper

    inline void compute_strides() {
      if (dimension == 1) {
	stride_(0) = 1;
      } else {
	size_t s = 1;
	for (size_t n=0; n<dimension; ++n) {
	  stride_(n) = s;
	  s *= extent_(n);
	}
      }
    }

  public: // Construction / destruction

    /// Default init (no operation)
    explicit uniform_grid_container(subscript_t extent = subscript_t()) {
      SL_REQUIRE("Good extent", extent == subscript_t() || extent.element_count() > 0);
      extent_ = extent;
      data_.resize(extent.element_count());
      compute_strides();
    }

    ~uniform_grid_container() {
      // Smart pointer handles delete
    }

    void resize(const subscript_t& sz) {
      if (extent_ != sz) {
	size_t n = sz.element_count();
	if (n == 0) {
	  /// Wipe out
	  data_.clear();
	  extent_ = subscript_t();
	  stride_ = subscript_t();
	} else {
	  data_.resize(n);
	  extent_ = sz;
	  compute_strides();
	  clear();
	}
      }
    }
    
  public: // Indexed access

    /// Offset into array
    inline size_t offset(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      size_t result = 0;
      for (size_t r=0; r<dimension; ++r) {
	result += idx(r) * stride_(r);
      }
      return result;
    }

    /// Is i a valid rank index for this?
    inline bool good_rank_index(size_t i) const {
      return i<dimension;
    }
    /// Is idx a good subscript?
    inline bool good_subscript(const subscript_t& idx) const {
      return idx.all_lt(extent());
    }

    /// The total number of elements
    inline size_t count() const {
      return extent().element_count();
    }

    /// The number of indexed elements for each of the ranks
    inline const subscript_t& extent() const {
      return extent_;
    }

    /// The element referenced by subscript idx
    inline const object_list_t& item(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return data_[offset(idx)];
    }

    /// The element referenced by subscript idx
    inline object_list_t& item(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return data_[offset(idx)];
    }

    /// The element referenced by subscript idx
    inline const object_list_t& operator()(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return item(idx);
    }

    /// The element referenced by subscript idx (alias of operator())
    inline const object_list_t& operator[](const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return item(idx);
    }

    /// The element referenced by subscript idx
    inline object_list_t& operator()(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return item(idx);
    }

    /// The element referenced by subscript idx (alias of operator())
    inline object_list_t& operator[](const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return item(idx);
    }

  public:

    void clear() {
      data_.clear();
      data_.resize(extent().element_count(),
                   object_list_t());
    }

    void insert(const subscript_t& idx, 
                object_t& o) {
      if (good_subscript(idx)) {
        item(idx).push_front(o);
      }
    }
    
    void insert(const subscript_t& l,
                const subscript_t& h,
                object_t& o) {
      for (subscript_t idx = l; idx.all_le(h); idx.increment(l,h)) {
        insert(idx, o);
      }
    }
    
  public:

    const aabox_t& bounding_box() const {
      return bbox_;
    }

    const vector_t& voxel_dimensions() const {
      return voxel_dim_;
    }
    
    void set_bounding_box(const aabox_t& b) {
      bbox_ = b;
      for (std::size_t r=0; r<dimension; ++r) {
        voxel_dim_[r] = (bbox_[1][r]-bbox_[0][r]) / value_t(extent()[r]);
      }
    }

    subscript_t subscript(const point_t& p) const {
      subscript_t result;
      for (std::size_t r=0; r<dimension; ++r) {
        result[r] = (int)(value_t(extent()[r]) * (p[r]-bbox_[0][r])/(bbox_[1][r]-bbox_[0][r]));
      }
      return result;
    }
    
    std::pair<subscript_t,subscript_t> subscript_range(const aabox_t& b) const {
      return std::make_pair(subscript(b[0]), subscript(b[1]));
    }

    void insert(const aabox_t& b, object_t& o) {
      insert(subscript(b[0]), subscript(b[1]), o);
    }

    void insert(const point_t& p, object_t& o) {
      insert(subscript(p), o);
    }

    void resize(const aabox_t& b, const value_t& cell_sz) {
      const std::size_t MAX_GRID_CELLS = 256*256*256; // FIXME parameterize
      const std::size_t MIN_GRID_CELLS = std::size_t(std::pow(double(3),double(dimension))); // FIXME parameterize
      clear();

      subscript_t gsz;
      for (std::size_t i=0; i<dimension; ++i) {
        gsz[i] = std::max(int(b.diagonal()[i]/cell_sz), 1);
      }
      if (gsz.element_count() > MAX_GRID_CELLS) {
        resize(b, 2.0f*cell_sz);
      } else if (gsz.element_count() < MIN_GRID_CELLS) {
        // FIXME - use box size
        for (std::size_t i=0; i<dimension; ++i) {
          gsz[i] = 3;
        }
        resize(gsz);
      } else {
        resize(gsz);
      }
      
      // This also recomputes voxels
      set_bounding_box(b);
    }

    value_t squared_point_cell_distance(const point_t& p,
                                        const subscript_t& c) const {
      value_t result = scalar_math<value_t>::zero();
      for (size_t j=0; j<dimension; ++j) {
        value_t l = bbox_[0][j]+(c[j]+0)*voxel_dim_[j];
        value_t h = bbox_[0][j]+(c[j]+1)*voxel_dim_[j];
	if (p[j] < l) {
	  result += sl::sqr(p[j] - l);
	} else if (p[j] > h) {
	  result += sl::sqr(p[j] - h);
	}
      }
      return result;
    }
                
  }; //end class uniform_grid_container

} // end namespace

#endif
