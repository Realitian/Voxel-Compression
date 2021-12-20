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
#ifndef SL_OCTREE_INDEXING_HPP
#define SL_OCTREE_INDEXING_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/integer.hpp>
#include <sl/utility.hpp>
#include <cassert>

namespace sl {

  template< int Bits = 64 > 
  class octree_index {
  public:
    // Note: index is signed to support unlimited addressing of neighbors
    typedef typename sl::int_t<Bits>::least  value_t; 
    typedef sl::fixed_size_point<3,value_t>  location_t;
  protected:
    std::size_t level_;
    location_t  location_;
  public:
    inline octree_index(std::size_t level = 0,
			const location_t& location = location_t()) :
      level_(level), location_(location) {
    }

    inline std::size_t level() const {
      return level_;
    }

    inline const location_t& location() const {
      return location_;
    }

  };

  template< int Bits = 64 > 
  struct octree_location_hasher {
    typedef octree_index<Bits> index_t;
    typedef typename index_t::location_t location_t;

    inline std::size_t operator()(const location_t& x) const {
      return sl::hash_bytes<sizeof(x)>((const unsigned char*)&(x));
    }
  };

  template< int Bits = 64 > 
  struct octree_index_hasher { 
    typedef octree_index<Bits> index_t;
    typedef typename index_t::location_t location_t;

    inline std::size_t operator()(const index_t& x) const {
      return 
	sl::hash_bytes<sizeof(location_t)>((const unsigned char*)&(x.location())) +
	std::size_t(13)*sl::hash_bytes<sizeof(std::size_t)>((const unsigned char*)&(x.level()));
    }
  };

  /**
   * Indexing operations on octrees of a given max level nl.   
   * Cells are indexed by their min corner with integer values
   * going from 0 to 2^(nl-1).
   */
  template< int Bits = 64 > 
  class octree_indexing {
  public:
    typedef octree_index<Bits>           index_t; 
    typedef typename index_t::value_t    value_t;
    typedef typename index_t::location_t location_t;
  protected:
    std::size_t max_level_count_;
  public:

    inline octree_indexing(std::size_t nl) :
      max_level_count_(nl) {
      assert(nl>0);
    }

  public:

    /// The maximum number of ocree levels
    inline std::size_t max_level_count() const {
      return max_level_count_;
    }

    /// True iff l is a valid level index
    inline bool is_good_level(std::size_t l) const {
      return l>=0 && l<max_level_count_;
    }

    /// Number of octree cells in one direction at level l
    static inline value_t level_grid_width(std::size_t l) {
      return (value_t(1)<<value_t(l));
    }

    /// Number of octree cells at finest level
    inline value_t grid_width() const {
      return level_grid_width(max_level_count_-1);
    }

    /// Cell width at level l
    inline value_t level_cell_width(std::size_t l) const {
      return (value_t(1)<<value_t(max_level_count_-1-l));
    }

  public:

    /// True iff x is within bounds
    bool is_good_location(const location_t& x) const {
      const value_t lo = 0;
      const value_t hi = grid_width()-1;

      return 
	(x[0]>=lo && x[0]<=hi) &&
	(x[1]>=lo && x[1]<=hi) &&
	(x[2]>=lo && x[2]<=hi);
    }
	
    /// Clamp location to octree grid domain
    location_t clamped_location(const location_t& x) const {
      const value_t lo = 0;
      const value_t hi = grid_width()-1;

      return location_t(sl::median(x[0],lo,hi),
			sl::median(x[1],lo,hi),
			sl::median(x[2],lo,hi));
    }

  public:

    /// Snap coordinates to grid at level l
    inline location_t snap(const location_t& x, std::size_t l) const {
      assert(is_good_location(x));

      const value_t s = value_t(max_level_count_-l-1);
      return location_t((x[0]>>s)<<s,
			(x[1]>>s)<<s,
			(x[2]>>s)<<s);
    }

  public:

    /**
     * Child offset containing x for a cell at level l. 
     */
    inline std::size_t child_offset(const location_t& x, std::size_t l) {
      std::size_t child_bit = (std::size_t(1)<<std::size_t(max_level_count_-l-2));
      std::size_t result = 
	((x[2] & child_bit) != 0) << 2
	|
	((x[1] & child_bit) != 0) << 1
	|
	((x[0] & child_bit) != 0);
      return result;
    }

  public:
    
    /// Id of the parent of (l,x)
    inline index_t cell_parent_id(const index_t& lx) const {
      return index_t(lx.level()-1,snap(lx.location(),lx.level()-1));
    }

    /// Id of the child_offset's child of (l,x)
    inline index_t cell_child_id(const index_t& lx,
				 std::size_t child_offset) const {
      value_t dsz =  level_cell_width(lx.level()+1);
      location_t result = lx.location();
      result[0] += (child_offset&1) ? dsz : value_t(0);
      result[1] += (child_offset&2) ? dsz : value_t(0);
      result[2] += (child_offset&4) ? dsz : value_t(0);
      return index_t(lx.level()+1,result);
    }

    inline index_t cell_neighbor_id(const index_t& lx,
				    const location_t& neighbor_offset) const {
      value_t dsz = level_cell_width(lx.level());
      location_t result = lx.location();
      result[0] += (neighbor_offset[0]>0) ? dsz : -dsz;
      result[1] += (neighbor_offset[1]>0) ? dsz : -dsz;
      result[2] += (neighbor_offset[2]>0) ? dsz : -dsz;
      return index_t(lx.level(),result);
      
    }

  public:

    inline location_t cell_vertex_id(const index_t& lx,
				     std::size_t vertex_idx) const {
      value_t dsz = level_cell_width(lx.level());
      location_t result = lx.location();
      result[0] += ((vertex_idx&1)!=0) ? dsz : 0;
      result[1] += ((vertex_idx&2)!=0) ? dsz : 0;
      result[2] += ((vertex_idx&4)!=0) ? dsz : 0;
      return result;
    }

    inline location_t cell_vertex_id(const index_t& lx,
				     const location_t& vertex_offset) const {
      value_t dsz = level_cell_width(lx.level());
      location_t result = lx.location();
      result[0] += (vertex_offset[0]>0) ? dsz : 0;
      result[1] += (vertex_offset[1]>0) ? dsz : 0;
      result[2] += (vertex_offset[2]>0) ? dsz : 0;
      return result;
    }

    inline location_t cell_split_vertex_id(const index_t& lx) const {
      value_t dsz = level_cell_width(lx.level()+1);
      location_t result = lx.location();
      result[0] += dsz;
      result[1] += dsz;
      result[2] += dsz;
      return result;
    }

    inline location_t cell_face_split_vertex_id(const index_t& lx,
						const location_t& face_offset) const {
      value_t dsz = level_cell_width(lx.level()+1);
      location_t result = lx.location();
      result[0] += dsz + ((face_offset[0]>0) ? dsz : -dsz);
      result[1] += dsz + ((face_offset[1]>0) ? dsz : -dsz);
      result[2] += dsz + ((face_offset[2]>0) ? dsz : -dsz);
      return result;
    }

  public: // Bounds and box intersection tests

    inline std::pair<location_t,location_t> cell_bounds(const index_t& lx) const {
      value_t dsz = level_cell_width(lx.level());
      const location_t& v0 = lx.location();
      return std::make_pair(v0,
			    location_t(v0[0]+dsz, v0[1]+dsz, v0[2]+dsz));
    }

    /// True iff box xlo,xhi has no intersection with cell lx
    inline bool is_disjoint(const location_t& lx_lo,
			    const location_t& lx_hi,
			    const location_t& xlo,
			    const location_t& xhi) const {
      return
	(lx_hi[0]<=xlo[0] || xhi[0]<=lx_lo[0]) ||
	(lx_hi[1]<=xlo[1] || xhi[1]<=lx_lo[1]) ||
	(lx_hi[2]<=xlo[2] || xhi[2]<=lx_lo[2]);
    }

    /// True iff box xlo,xhi has some intersection with cell lx
    inline bool is_overlapping(const location_t& lx_lo,
			       const location_t& lx_hi,
			       const location_t& xlo,
			       const location_t& xhi) const {
      return !is_disjoint(lx_lo, lx_hi, xlo, xhi);
    }

    /// True iff box xlo,xhi is entirely within cell lx
    inline bool is_containing(const location_t& lx_lo,
			      const location_t& lx_hi,
			      const location_t& xlo,
			      const location_t& xhi) const {
      return
	(lx_lo[0]<=xlo[0] && xhi[0]<=lx_hi[0]) &&
	(lx_lo[1]<=xlo[1] && xhi[1]<=lx_hi[1]) &&
	(lx_lo[2]<=xlo[2] && xhi[2]<=lx_hi[2]);
    }

    /// True iff cell lx is entirely withing box xlo,xhi
    inline bool is_contained(const location_t& lx_lo,
			     const location_t& lx_hi,
			     const location_t& xlo,
			     const location_t& xhi) const {
      return
	(xlo[0]<=lx_lo[0] && lx_hi[0]<=xhi[0]) &&
	(xlo[1]<=lx_lo[1] && lx_hi[1]<=xhi[1]) &&
	(xlo[2]<=lx_lo[2] && lx_hi[2]<=xhi[2]);
    }
    

    /// True iff box xlo,xhi has no intersection with cell lx
    inline bool is_disjoint(const index_t& lx,
			    const location_t& xlo,
			    const location_t& xhi) const {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(lx);
      return is_disjoint(lx_bounds.first, lx_bounds.second, xlo, xhi);
    }

    /// True iff box xlo,xhi has some intersection with cell lx
    inline bool is_overlapping(const index_t& lx,
			       const location_t& xlo,
			       const location_t& xhi) const {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(lx);
      return is_overlapping(lx_bounds.first, lx_bounds.second, xlo, xhi);
    }

    /// True iff box xlo,xhi is entirely within cell lx
    inline bool is_containing(const index_t& lx,
			      const location_t& xlo,
			      const location_t& xhi) const {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(lx);
      return is_containing(lx_bounds.first, lx_bounds.second, xlo, xhi);
    }

    /// True iff cell lx is entirely withing box xlo,xhi
    inline bool is_contained(const index_t& lx,
			     const location_t& xlo,
			     const location_t& xhi) const {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(lx);
      return is_contained(lx_bounds.first, lx_bounds.second, xlo, xhi);
    }
 
  }; // class octree_indexing

} // namespace sl

namespace sl {

  typedef octree_index< 8> octree_index8_t; 
  typedef octree_index<16> octree_index16_t; 
  typedef octree_index<32> octree_index32_t; 
  typedef octree_index<64> octree_index64_t; 

  typedef octree_location_hasher< 8> octree_location_hasher8_t; 
  typedef octree_location_hasher<16> octree_location_hasher16_t; 
  typedef octree_location_hasher<32> octree_location_hasher32_t; 
  typedef octree_location_hasher<64> octree_location_hasher64_t; 

  typedef octree_index_hasher< 8> octree_index_hasher8_t; 
  typedef octree_index_hasher<16> octree_index_hasher16_t; 
  typedef octree_index_hasher<32> octree_index_hasher32_t; 
  typedef octree_index_hasher<64> octree_index_hasher64_t; 

  typedef octree_indexing< 8> octree_indexing8_t; 
  typedef octree_indexing<16> octree_indexing16_t; 
  typedef octree_indexing<32> octree_indexing32_t; 
  typedef octree_indexing<64> octree_indexing64_t; 
}

#endif
