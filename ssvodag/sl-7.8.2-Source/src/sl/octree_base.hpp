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
#ifndef SL_OCTREE_BASE_HPP
#define SL_OCTREE_BASE_HPP

#include <sl/octree_indexing.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/utility.hpp> 
#include <cassert>

#include <sl/fsb_allocator.hpp>

namespace sl {
  
  template <typename G_Cell_Data, typename G_Derived>
  class octree_cell_base {
  public:
    typedef octree_cell_base<G_Cell_Data, G_Derived> this_t;
    typedef G_Derived                                derived_t;
    typedef G_Cell_Data                              data_t;
  protected:
    derived_t       *parent_;    // Pointer to parent cell
    derived_t       *children_;  // Pointer to first of 8 contiguous child cells
    data_t          data_;       // Application specific cell data
  public:
    SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);
  public:
    
    octree_cell_base() : 
      parent_(0), children_(0) {
    }
    
    const derived_t* parent() const {
      return parent_;
    }

    derived_t* parent() {
      return parent_;
    }

    void set_parent(derived_t* p) {
      parent_ = p;
    }

    const derived_t* children() const {
      return children_;
    }

    derived_t* children() {
      return children_;
    }

    void set_children(derived_t* c) {
      children_ = c;
    }

    const derived_t* child(std::size_t k) const {
      assert(children_);
      assert(k<8);
      return &children_[k];
    }

    derived_t* child(std::size_t k) {
      assert(children_);
      assert(k<8);
      return &children_[k];
    }

    void set_child(std::size_t k,
		   derived_t* c) {
      assert(children_);
      assert(k<8);
      children_[k] = c;
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
  }; // octree_cell_base
} // namespace sl

namespace sl {

  template <
    typename G_Cell
    >
  class octree_cell_pointer {
  public:
    typedef G_Cell                           cell_t;
    typedef octree_index32_t                 index_t;
    typedef octree_index32_t::location_t     location_t;
  protected:
    index_t index_;
    cell_t* cell_;
  public:
    octree_cell_pointer(const index_t& idx = index_t(),
			cell_t* ptr = 0) :
      index_(idx), cell_(ptr) {
    }

    const index_t& index() const {
      return index_;
    }

    std::size_t level() const {
      return index_.level();
    }

    location_t location() const {
      return index_.location();
    }

    cell_t* cell() const {
      return cell_;
    }

    bool is_null() const {
      return cell_ == 0;
    }

    bool is_leaf() const {
      return (!is_null()) && (cell()->children() == 0);
    }

    bool is_root() const {
      return (!is_null()) && (cell()->parent() == 0);
    }
      
  }; 

}  // namespace sl

namespace sl {

  template <typename G_Cell>
  class octree_base: public octree_indexing<32> {
  public:
    typedef octree_base<G_Cell>                this_t;
    typedef octree_indexing<32>                super_t;
    typedef G_Cell                             cell_t;
    typedef octree_index32_t                   index_t;
    typedef octree_index32_t::location_t       location_t;
    typedef octree_cell_pointer<const cell_t>  const_cell_pointer_t;
    typedef octree_cell_pointer<cell_t>        cell_pointer_t;
  protected:
    cell_t *                 root_;
    std::size_t              cell_count_;
    std::size_t              leaf_cell_count_;

  protected: // Callbacks

    virtual void initialize_cell(const cell_pointer_t& /*ptr*/) {
    }

    virtual void finalize_cell(const cell_pointer_t& /*ptr*/) {
    }

  protected: // Memory allocation

    typedef cell_t cell_children_t[8];

    fsb_allocator<cell_children_t> cell_children_allocator_;

    inline void root_new() {
      assert(!root_);
      root_ = new cell_t;
      this->initialize_cell(cell_pointer_t(index_t(0, location_t(0,0,0)),
						   root_));
      ++cell_count_;
      ++leaf_cell_count_;
    }
 
    inline void root_delete() {
      assert(root_);
      this->finalize_cell(cell_pointer_t(index_t(0, location_t(0,0,0)),
						   root_));
      delete root_; root_ = 0;
      --cell_count_;
      --leaf_cell_count_;
    }

    inline void cell_children_new(const index_t& parent_index, 
				  cell_t*  parent) {
#if 0
      cell_t* children = new cell_t[8];
#else
      cell_t* children = reinterpret_cast<cell_t*>(cell_children_allocator_.allocate(1));
      for (std::size_t i=0; i<8; ++i) {
	new(&(children[i])) cell_t();
      }
#endif
      for (std::size_t i=0; i<8; ++i) {
	children[i].set_parent(parent);
      }
      parent->set_children(children);

      for (std::size_t i=0; i<8; ++i) {
	this->initialize_cell(cell_pointer_t(this->cell_child_id(parent_index, i),
					     &(children[i])));
	++cell_count_;
      }
      leaf_cell_count_ += 7; // -1 + 8
    }

    inline void cell_children_delete(const index_t& parent_index,
				     cell_t* parent) {
      cell_t* children = parent->children();
      if (children) {
	for (std::size_t i=0; i<8; ++i) {
	  this->finalize_cell(cell_pointer_t(this->cell_child_id(parent_index, i),
					     &(children[i])));
	  --cell_count_;
	}
	leaf_cell_count_ -= 7; // -1 + 8

	parent->set_children(0);

#if 0
	delete[] children;
#else
	for (std::size_t i=0; i<8; ++i) {
	  children[i].~cell_t();
	}
	cell_children_allocator_.deallocate(reinterpret_cast<cell_children_t*>(children), 1);
#endif
      }
    }
    
  private:

    octree_base(const this_t& other);       // Masked: not implemented
    this_t& operator=(const this_t& other); // Masked: not implemented

  public:

    octree_base(std::size_t nl) :
      super_t(nl), root_(0), cell_count_(0), leaf_cell_count_(0) {
      assert(nl>0);				
    }

    virtual ~octree_base() {
      if (root_) {
	clear();
      }
    }

  public: // Cell count

    std::size_t cell_count() const {
      return cell_count_;
    }

    std::size_t leaf_cell_count() const {
      assert(leaf_cell_count_<=cell_count_);

      return leaf_cell_count_;
    }

    bool is_empty() const {
      return cell_count_ == 0;
    }

  public: // Root

    const_cell_pointer_t root() const {
      return const_cell_pointer_t(index_t(0,location_t(0,0,0)), root_);
    }

    cell_pointer_t root() {
      return cell_pointer_t(index_t(0,location_t(0,0,0)), root_);
    }

    const_cell_pointer_t cell_child_pointer(const const_cell_pointer_t& ptr,
					    std::size_t i) const {
      return const_cell_pointer_t(cell_child_id(ptr.index(), i),
				  ptr.cell()->child(i));
    }

    cell_pointer_t cell_child_pointer(const cell_pointer_t& ptr,
				      std::size_t i) {
      return cell_pointer_t(cell_child_id(ptr.index(), i),
			    ptr.cell()->child(i));
    }

    const_cell_pointer_t cell_parent_pointer(const const_cell_pointer_t& ptr) const {
      return const_cell_pointer_t(cell_parent_id(ptr.index()), ptr.cell()->parent());
    }

    cell_pointer_t cell_parent_pointer(const cell_pointer_t& ptr) {
      return cell_pointer_t(cell_parent_id(ptr.index()), ptr.cell()->parent());
    }

  public: // Point, cell or region location

    /**
     * Locate leaf of current octree containing x
     */
    cell_pointer_t locate(const location_t& x) {
      assert(is_good_location(x));

      cell_t* current_cell = root_;
      std::size_t current_level = 0; 
      if (current_cell) {
	while (current_cell->children()) {
	  std::size_t child_offset = this->child_offset(x, current_level);
	  current_cell = (current_cell->child(child_offset));
	  ++current_level;
	}
      }
      return cell_pointer_t(index_t(current_level, this->snap(x, current_level)), 
			    current_cell);
    }

    /**
     * Locate leaf of current octree containing x
     */
    const_cell_pointer_t locate(const location_t& x) const {
      assert(is_good_location(x));

      cell_pointer_t ptr = const_cast<this_t*>(this)->locate(x);
      return const_cell_pointer_t(ptr.index(), ptr.cell());
    }
 
    /**
     * Locate cell of current octree fully containing the
     * box xlo..xhi
     */
    cell_pointer_t locate(const location_t& xlo,
			  const location_t& xhi) {
      assert(is_good_location(xlo));
      assert(is_good_location(xhi));

      cell_t* current_cell = root_;
      std::size_t current_level = 0;
      for (;;) {
	if (!current_cell->children()) {
	  break;
	} else {
	  std::size_t child_offset_lo = this->child_offset(xlo, current_level);
	  std::size_t child_offset_hi = this->child_offset(xhi, current_level);
	  if (child_offset_lo!= child_offset_hi) {
	    break;
	  } else {
	    current_cell = (current_cell->child(child_offset_lo));
	    ++current_level;
	  }
	}
      }
      return cell_pointer_t(index_t(current_level, 
				    this->snap(xlo, current_level)), 
			    current_cell);
    }
    
   /**
     * Locate cell of current octree fully containing the
     * box xlo..xhi
     */
    const_cell_pointer_t locate(const location_t& xlo,
				const location_t& xhi) const {
      assert(is_good_location(xlo));
      assert(is_good_location(xhi));
      cell_pointer_t ptr = const_cast<this_t*>(this)->locate(xlo,xhi);
      return const_cell_pointer_t(ptr.index(), ptr.cell());
    }

    /**
     * Locate leaf of current octree containing cell lx 
     * (i.e., either lx or a parent)
     */
    cell_pointer_t locate(const index_t& lx) {
      assert(is_good_location(lx.location()));

      cell_t* current_cell = root_;
      std::size_t current_level = 0; 
      if (current_cell) {
	while (current_level<lx.level() && current_cell->children()) {
	  std::size_t child_offset = this->child_offset(lx.location(), current_level);
	  current_cell = (current_cell->child(child_offset));
	  ++current_level;
	}
      }
      return cell_pointer_t(index_t(current_level, 
				    this->snap(lx.location(), current_level)), 
			    current_cell);
    }

    /**
     * Locate leaf of current octree containing cell lx 
     * (i.e., either lx or a parent)
     */
    const_cell_pointer_t locate(const index_t& lx) const {
      assert(is_good_location(lx.location()));

      cell_pointer_t ptr = const_cast<this_t*>(this)->locate(lx);
      return const_cell_pointer_t(ptr.index(), ptr.cell());
    }

  public: // Refine and coarsen primitives

    /// Create children of pointed cell
    void refine(const cell_pointer_t& ptr) {
      assert(!ptr.is_null());
      assert(ptr.level()+1<max_level_count());

      if (ptr.is_leaf()) {
	this->cell_children_new(ptr.index(), ptr.cell());
      }
    }

    /// Remove children of pointed cell
    void coarsen(const cell_pointer_t& ptr) {
      assert(!ptr.is_null());
      if (ptr.is_leaf()) {
	// Do nothing 
      } else {
	// Recurse on children...
	for (std::size_t i=0; i<8; ++i) {
	  this->coarsen(cell_child_pointer(ptr, i));
	}
	// ... and then delete them
	this->cell_children_delete(ptr.index(), ptr.cell());
      }
    }

  public: // Clear

    virtual void clear() {
      if (root_) {
	this->coarsen(root());
	root_delete();
      }
    }
  public: // Region refinement operations
 
    /**
     *  Locate and eventually create all cells leading to lx
     */
    cell_pointer_t refine_locate(const index_t& lx) {
      assert(is_good_location(lx.location()));

      if (!root_) root_new();

      cell_t* current_cell = root_;
      std::size_t current_level = 0; 
      while (current_level<lx.level()) {
	if (!current_cell->children()) {
	  refine(cell_pointer_t(index_t(current_level, 
					this->snap(lx.location(), current_level)), 
				current_cell));
	}
	std::size_t child_offset = this->child_offset(lx.location(), current_level);
	current_cell = (current_cell->child(child_offset));
	++current_level;
      }
      return cell_pointer_t(index_t(current_level, 
				    this->snap(lx.location(), current_level)), 
			    current_cell);
    }

    /**
     *  Locate and eventually create all cells leading to x
     */
    cell_pointer_t refine_locate(const location_t& x) {
      assert(is_good_location(x));

      if (!root_) root_new();

      cell_t* current_cell = root_;
      std::size_t current_level = 0; 
      while (current_level+1<max_level_count_) {
	if (!current_cell->children()) {
	  refine(cell_pointer_t(index_t(current_level, 
					this->snap(x, current_level)), 
				current_cell));
	}
	std::size_t child_offset = this->child_offset(x, current_level);
	current_cell = (current_cell->child(child_offset));
	++current_level;
      }
      return cell_pointer_t(index_t(current_level, 
				    this->snap(x, current_level)), 
			    current_cell);
    }

    /**
     *  Refine all cells overlapping the box 
     *  xlo, xhi to level xlevel
     */
    void refine_overlapping(std::size_t xlevel,
			    const location_t& xlo,
			    const location_t& xhi) {
      if (!root_) { root_new(); }

      location_t cxlo = clamped_location(xlo);
      location_t cxhi = clamped_location(xhi);

      cell_pointer_t ptr = locate(cxlo, cxhi);
      assert(!ptr.is_null());
      
      this->refine_overlapping(ptr, 
			       xlevel, xlo, xhi);
    }

    /**
     *  Refine all cells overlapping the box 
     *  xlo, xhi to level xlevel
     */
    void refine_overlapping(const cell_pointer_t& ptr,
			    std::size_t xlevel,
			    const location_t& xlo,
			    const location_t& xhi) {
      assert(!ptr.is_null());
      if ((ptr.level()<xlevel) && is_overlapping(ptr.index(), xlo, xhi)) {
	if (ptr.is_leaf()) {
	  refine(ptr);
	}
	if (!ptr.is_leaf()) {
	  for (std::size_t i=0; i<8; ++i) {
	    refine_overlapping(cell_child_pointer(ptr, i),
			       xlevel,
			       xlo,
			       xhi);
	  }
	}
      }
    }

    /**
     *  Refine all cells fully contained in the box 
     *  xlo, xhi to level xlevel
     */
    void refine_contained(std::size_t xlevel,
			  const location_t& xlo,
			  const location_t& xhi) {
      if (!root_) { root_new(); }

      location_t cxlo = clamped_location(xlo);
      location_t cxhi = clamped_location(xhi);

      cell_pointer_t ptr = locate(cxlo, cxhi);
      assert(!ptr.is_null());
      
      this->refine_contained(ptr, 
			     xlevel, xlo, xhi);
    }

    /**
     *  Refine all cells fully contained  in the box 
     *  xlo, xhi to level xlevel
     */
    void refine_contained(const cell_pointer_t& ptr,
			  std::size_t xlevel,
			  const location_t& xlo,
			  const location_t& xhi) {
      assert(!ptr.is_null());
      if (ptr.level()<xlevel) {
	std::pair<location_t,location_t> lx_bounds = cell_bounds(ptr.index());
	
	if (is_overlapping(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	  // No overlap, stop
	} else {
	  if (ptr.is_leaf() && is_contained(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	    refine(ptr);
	  }
	  if (!ptr.is_leaf()) {
	    for (std::size_t i=0; i<8; ++i) {
	      refine_contained(cell_child_pointer(ptr,i),
			       xlevel,
			       xlo,
			       xhi);
	    }
	  }
	}
      }
    }

    /**
     *  Coarsen all cells fully contained in the box 
     *  xlo, xhi to level xlevel
     */
    void coarsen_contained(std::size_t xlevel,
			  const location_t& xlo,
			  const location_t& xhi) {
      if (!root_) { root_new(); }

      location_t cxlo = clamped_location(xlo);
      location_t cxhi = clamped_location(xhi);

      cell_pointer_t ptr = locate(cxlo, cxhi);
      assert(!ptr.is_null());
      
      this->coarsen_contained(ptr, 
			      xlevel, xlo, xhi);
    }

    /**
     *  Coarsen all cells fully contained  in the box 
     *  xlo, xhi to level xlevel
     */
    void coarsen_contained(const cell_pointer_t& ptr,
			  std::size_t xlevel,
			  const location_t& xlo,
			  const location_t& xhi) {
      assert(!ptr.is_null());

      if (is_disjoint(ptr.index(), xlo, xhi)) {
	// No overlap, stop
      } else if (ptr.level()>=xlevel && !ptr.is_leaf() && is_contained(ptr.index(), xlo, xhi)) {
	coarsen(ptr);
      } else {
	if (!ptr.is_leaf()) {
	  for (std::size_t i=0; i<8; ++i) {
	    coarsen_contained(cell_child_pointer(ptr,i),
			      xlevel, xlo, xhi);
	  }
	}
      }
    }
	   
  public: // Visit all leaf cells

    template <typename Visitor>
    void visit_all_leaf_cells(Visitor& visitor) const {
      if (root_) {
	visit_all_leaf_cells(root(), visitor);
      }
    }

    template <typename Visitor>
    void visit_all_leaf_cells(const const_cell_pointer_t& ptr,
			      Visitor& visitor) const {
      if (ptr.is_leaf()) {
	visitor(ptr);
      } else {
	for (std::size_t i=0; i<8; ++i) {
	  visit_all_leaf_cells(cell_child_pointer(ptr,i),
			       visitor);
	}
      }	  
    }

    template <typename Visitor>
    void visit_all_leaf_cells(Visitor& visitor) {
      if (root_) {
	visit_all_leaf_cells(root(), visitor);
      }
    }

    template <typename Visitor>
    void visit_all_leaf_cells(const cell_pointer_t& ptr,
			      Visitor& visitor) {
      if (ptr.is_leaf()) {
	visitor(ptr);
      } else {
	for (std::size_t i=0; i<8; ++i) {
	  visit_all_leaf_cells(cell_child_pointer(ptr,i),
			       visitor);
	}
      }	  
    }

  public: // Visit contained leaf cells

    template <typename Visitor>
    void visit_contained_leaf_cells(const location_t& xlo,
				    const location_t& xhi,
				    Visitor& visitor) const {
      if (root_) {
	location_t cxlo = clamped_location(xlo);
	location_t cxhi = clamped_location(xhi);

	const_cell_pointer_t ptr = locate(cxlo, cxhi);
	assert(!ptr.is_null());
	visit_contained_leaf_cells(ptr, xlo, xhi, visitor);
      }
    }

    template <typename Visitor>
    void visit_contained_leaf_cells(const const_cell_pointer_t& ptr,
				    const location_t& xlo,
				    const location_t& xhi,
				    Visitor& visitor) const {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(ptr.index());
      if (is_overlapping(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	if (ptr.is_leaf()) {
	  if (is_contained(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	    visitor(ptr);
	  }
	} else {
	  for (std::size_t i=0; i<8; ++i) {
	    visit_contained_leaf_cells(cell_child_pointer(ptr,i),
				       xlo, xhi, visitor);
	  }
	}
      }	  
    }

    template <typename Visitor>
    void visit_contained_leaf_cells(const location_t& xlo,
				    const location_t& xhi,
				    Visitor& visitor) {
      if (root_) {
	location_t cxlo = clamped_location(xlo);
	location_t cxhi = clamped_location(xhi);

	cell_pointer_t ptr = locate(cxlo, cxhi);
	assert(!ptr.is_null());
	visit_contained_leaf_cells(ptr, xlo, xhi, visitor);
      }
    }

    template <typename Visitor>
    void visit_contained_leaf_cells(const cell_pointer_t& ptr,
				    const location_t& xlo,
				    const location_t& xhi,
				    Visitor& visitor) {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(ptr.index());

      if (is_overlapping(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	if (ptr.is_leaf()) {
	  if (is_contained(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	    visitor(ptr);
	  }
	} else {
	  for (std::size_t i=0; i<8; ++i) {
	    visit_contained_leaf_cells(cell_child_pointer(ptr,i),
				       xlo, xhi, visitor);
	  }
	}
      }	  
    }

    std::size_t contained_leaf_cells_count(const location_t& xlo,
					   const location_t& xhi) const {
      counting_visitor cv;
      visit_contained_leaf_cells(xlo, xhi, cv);
      return cv.count();
    }

    std::size_t contained_leaf_cells_count(const const_cell_pointer_t& ptr,
					   const location_t& xlo,
					   const location_t& xhi) const {
      counting_visitor cv;
      visit_contained_leaf_cells(ptr, xlo, xhi, cv);
      return cv.count();
    }
 
  public: // Visit overlapping leaf cells

    template <typename Visitor>
    void visit_overlapping_leaf_cells(const location_t& xlo,
				      const location_t& xhi,
				      Visitor& visitor) const {
      if (root_) {
	location_t cxlo = clamped_location(xlo);
	location_t cxhi = clamped_location(xhi);

	const_cell_pointer_t ptr = locate(cxlo, cxhi);
	assert(!ptr.is_null());
	visit_overlapping_leaf_cells(ptr, xlo, xhi);
      }
    }

    template <typename Visitor>
    void visit_overlapping_leaf_cells(const location_t& xlo,
				      const location_t& xhi,
				      Visitor& visitor) {
      if (root_) {
	location_t cxlo = clamped_location(xlo);
	location_t cxhi = clamped_location(xhi);

	cell_pointer_t ptr = locate(cxlo, cxhi);
	assert(!ptr.is_null());
	visit_overlapping_leaf_cells(ptr, cxlo, cxhi, visitor);
      }
    }

    template <typename Visitor>
    void visit_overlapping_leaf_cells(const const_cell_pointer_t& ptr,
				      const location_t& xlo,
				      const location_t& xhi,
				      Visitor& visitor) const {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(ptr.index());
      if (is_overlapping(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	if (ptr.is_leaf()) {
	  visitor(ptr);
	} else if (is_contained(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	  for (std::size_t i=0; i<8; ++i) {
	    visit_all_leaf_cells(cell_child_pointer(ptr,i),
				 visitor);
	  }
	} else {
	  for (std::size_t i=0; i<8; ++i) {
	    visit_overlapping_leaf_cells(cell_child_pointer(ptr,i),
					 xlo, xhi, visitor);
	  }
	}
      }	  
    }

    template <typename Visitor>
    void visit_overlapping_leaf_cells(const cell_pointer_t& ptr,
				      const location_t& xlo,
				      const location_t& xhi,
				      Visitor& visitor) {
      std::pair<location_t,location_t> lx_bounds = cell_bounds(ptr.index());
      if (is_overlapping(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	if (ptr.is_leaf()) {
	  visitor(ptr);
	} else if (is_contained(lx_bounds.first, lx_bounds.second, xlo, xhi)) {
	  for (std::size_t i=0; i<8; ++i) {
	    visit_all_leaf_cells(cell_child_pointer(ptr,i),
				 visitor);
	  }
	} else {
	  for (std::size_t i=0; i<8; ++i) {
	    visit_overlapping_leaf_cells(cell_child_pointer(ptr,i),
					 xlo, xhi, visitor);
	  }
	}
      }	  
    }

    std::size_t overlapping_leaf_cells_count(const location_t& xlo,
					     const location_t& xhi) const {
      counting_visitor cv;
      visit_overlapping_leaf_cells(xlo, xhi, cv);
      return cv.count();
    }

    std::size_t overlapping_leaf_cells_count(const const_cell_pointer_t& ptr,
					     const location_t& xlo,
					     const location_t& xhi) const {
      counting_visitor cv;
      visit_overlapping_leaf_cells(ptr, xlo, xhi, cv);
      return cv.count();
    }

  }; // class octree_base
} // namespace sl

#endif
