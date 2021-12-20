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
#ifndef SL_CELL_OCTREE_HPP
#define SL_CELL_OCTREE_HPP

#include <sl/octree_base.hpp>
#include <cassert>


namespace sl {

  /// Cell of an octree with only cells
  template <typename G_Cell_Data>
  class cell_octree_cell: public octree_cell_base<G_Cell_Data, cell_octree_cell<G_Cell_Data> > {
  public:
    typedef cell_octree_cell<G_Cell_Data>                                 this_t;
    typedef G_Cell_Data                                                   data_t;
    typedef octree_cell_base<data_t, this_t> super_t;

  public:
    cell_octree_cell() {
    }
  }; // cell_octree_cell
} // namespace sl

namespace sl {

  template <typename G_Cell_Data>
  class cell_octree: public octree_base< cell_octree_cell<G_Cell_Data> > {
  public:
    typedef octree_base< cell_octree_cell<G_Cell_Data> >  super_t;
    typedef cell_octree_cell<G_Cell_Data>                 this_t;
    
    typedef typename super_t::cell_t cell_t;
    typedef typename super_t::index_t index_t;
    typedef typename super_t::location_t location_t;
    typedef typename super_t::const_cell_pointer_t const_cell_pointer_t;
    typedef typename super_t::cell_pointer_t cell_pointer_t;
  protected:

  public:

    cell_octree(std::size_t nl) :
      super_t(nl) {
      assert(nl>0);				
    }

    virtual ~cell_octree() {
    }

  };

} // namespace sl

#endif
