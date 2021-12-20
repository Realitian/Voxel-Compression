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
#ifndef SL_CELL_VERTEX_OCTREE_HPP
#define SL_CELL_VERTEX_OCTREE_HPP

#include <sl/cell_octree.hpp>
#include <sl/utility.hpp>
#include <map>
#include <cassert>

namespace sl {  

  template <typename G_Vertex_Data, typename G_Derived>
  class octree_vertex_base {
  public:
    typedef octree_vertex_base<G_Vertex_Data, G_Derived> this_t;
    typedef G_Derived                                    derived_t;
    typedef G_Vertex_Data                                data_t;
  protected:
    mutable int32_t refcount_;
    data_t      data_;      // Application specific vertex data
  public:
    SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);
  public:
    
    octree_vertex_base() : refcount_(0) {
    }

  public:

    inline int32_t refcount() const {
      return refcount_;
    }

    inline void ref() const {
      ++refcount_;
    }

    inline void deref() const {
      --refcount_;
    }

  public:

    const data_t& data() const {
      return data_;
    }

    data_t& data() {
      return data_;
    }

    void set_data(const data_t& x) {
      data_ = x;
    }

  }; // octree_vertex_base
} // namespace sl


namespace sl {

  /// Vertex of an octree with cells and vertices
  template <typename G_Vertex_Data>
  class cell_vertex_octree_vertex: public octree_vertex_base<G_Vertex_Data, cell_vertex_octree_vertex<G_Vertex_Data> > {
  public:
    typedef cell_vertex_octree_vertex<G_Vertex_Data> this_t;
    typedef G_Vertex_Data                            data_t;
    typedef octree_vertex_base<data_t,this_t>        super_t;
  public:
  
    cell_vertex_octree_vertex() {
    }

  }; // cell_vertex_octree_vertex
} // namespace sl

namespace sl {

  /// Cell of an octree with cells and vertices
  template <typename G_Cell_Data, typename G_Vertex_Data>
  class cell_vertex_octree_cell: 
    public octree_cell_base<G_Cell_Data, cell_vertex_octree_cell<G_Cell_Data, G_Vertex_Data> > {
  public:
    typedef cell_vertex_octree_cell<G_Cell_Data, G_Vertex_Data>  this_t;
    typedef G_Cell_Data                                            data_t;
    typedef octree_cell_base<data_t, this_t>                       super_t;

    typedef octree_index32_t                                       index_t;
    typedef octree_index32_t::location_t                           location_t;

    typedef G_Vertex_Data                                          vertex_data_t;
    typedef cell_vertex_octree_vertex<vertex_data_t>               vertex_t;
#if 0
    typedef std::map<location_t,vertex_t>                          vertex_map_t;
#else
    typedef std::map<location_t,vertex_t, std::less<location_t>, fsb_allocator< std::pair<location_t,vertex_t> > >  vertex_map_t;
#endif
    typedef typename vertex_map_t::iterator                        vertex_iterator_t;
    typedef typename vertex_map_t::const_iterator                  const_vertex_iterator_t;
  protected:
    vertex_iterator_t  vertex_iterator_[8]; // Iterators to vertex data
  public:

    cell_vertex_octree_cell() {
    }

    vertex_iterator_t vertex_iterator(std::size_t i) {
      assert(i<8);
      return vertex_iterator_[i];
    }

    const_vertex_iterator_t vertex_iterator(std::size_t i) const {
      assert(i<8);
      return vertex_iterator_[i];
    }

    void set_vertex_iterator(std::size_t i, const vertex_iterator_t& x) {
      assert(i<8);
      vertex_iterator_[i] = x;
    }

  }; // cell_vertex_octree_cell

} // namespace sl

namespace sl {

  template <typename G_Cell_Data, typename G_Vertex_Data>
  class cell_vertex_octree: public octree_base< cell_vertex_octree_cell<G_Cell_Data, G_Vertex_Data> > {
  public:
    typedef octree_base< cell_vertex_octree_cell<G_Cell_Data,G_Vertex_Data> > super_t;
    typedef cell_vertex_octree<G_Cell_Data, G_Vertex_Data> this_t;
    typedef typename super_t::cell_t cell_t;
    typedef typename super_t::index_t index_t;
    typedef typename super_t::location_t location_t;
    typedef typename super_t::const_cell_pointer_t const_cell_pointer_t;
    typedef typename super_t::cell_pointer_t cell_pointer_t;

    typedef G_Vertex_Data                                          vertex_data_t;
    typedef cell_vertex_octree_vertex<vertex_data_t>               vertex_t;
    typedef typename cell_t::vertex_map_t                          vertex_map_t;
    typedef typename vertex_map_t::iterator                        vertex_iterator_t;
    typedef typename vertex_map_t::const_iterator                  const_vertex_iterator_t;

  protected:

    vertex_map_t vertex_map_;

  public:

    cell_vertex_octree(std::size_t nl) :
      super_t(nl) {
      assert(nl>0);				
    }

    virtual ~cell_vertex_octree() {

    }

  protected: // Callbacks

    virtual void initialize_cell(const cell_pointer_t& ptr) {
      super_t::initialize_cell(ptr);

      // Create vertices of new cell and register new ones in vertex map
      for (std::size_t i=0; i<8; ++i) {
	location_t vid = this->cell_vertex_id(ptr.index(), i);
#if 0
	std::cerr << "L:" << ptr.level() << " C:" << ptr.location()[0] << " " << ptr.location()[1] << " " << ptr.location()[2] << " V" << i << "->" << vid[0] << " " << vid[1] << " " << vid[2] << std::endl;
#endif
	vertex_iterator_t vit = vertex_map_.find(vid);
	if (vit == vertex_map_.end()) {
	  // New vertex
	  vit = vertex_map_.insert(std::make_pair(vid,vertex_t())).first;
	  vit->second.ref();
	  initialize_vertex(vit);
	} else {
	  // Already present
	  vit->second.ref();
	}

	// Store in cell table
	ptr.cell()->set_vertex_iterator(i, vit);
      }
    }

    virtual void finalize_cell(const cell_pointer_t& ptr) {
      // Dereference vertices of destroyed cell
      for (std::size_t i=0; i<8; ++i) {
	location_t vid = this->cell_vertex_id(ptr.index(), i);

	vertex_iterator_t vit = vertex_map_.find(vid);
	assert(vit != vertex_map_.end());
	vit->second.deref();
	if (vit->second.refcount() <=0) {
	  finalize_vertex(vit);
	  vertex_map_.erase(vit);
	}
      }

      super_t::finalize_cell(ptr);
    }

    virtual void initialize_vertex(const vertex_iterator_t& /*vit*/) {
    }

    virtual void finalize_vertex(const vertex_iterator_t& /*vit*/) {
    }

  public:

    std::size_t vertex_count() const {
      return vertex_map_.size();
    }

    vertex_iterator_t vertex_begin() {
      return vertex_map_.begin();
    }

    vertex_iterator_t vertex_end() {
      return vertex_map_.end();
    }

    vertex_iterator_t vertex_find(const location_t& vid) {
      return vertex_map_.find(vid);
    }

    const_vertex_iterator_t vertex_begin() const {
      return vertex_map_.begin();
    }

    const_vertex_iterator_t vertex_end() const {
      return vertex_map_.end();
    }

    const_vertex_iterator_t vertex_find(const location_t& vid) const {
      return vertex_map_.find(vid);
    }

    vertex_data_t* vertex_data(const location_t& vid) {
      vertex_data_t* result = 0;
      vertex_iterator_t vit = vertex_find(vid);
      if (vit != vertex_end()) {
	result = &(vit.second->data());
      }
      return result;
    }

    const vertex_data_t* vertex_data(const location_t& vid) const {
      const vertex_data_t* result = 0;
      const_vertex_iterator_t vit = vertex_find(vid);
      if (vit != vertex_end()) {
	result = &(vit.second->data());
      }
      return result;
    }


  };


} // namespace sl

#endif
