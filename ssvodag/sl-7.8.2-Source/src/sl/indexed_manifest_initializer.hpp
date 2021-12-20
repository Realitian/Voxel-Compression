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
#ifndef SL_INDEXED_MANIFEST_INITIALIZER_HPP
#define SL_INDEXED_MANIFEST_INITIALIZER_HPP

#include <sl/index.hpp>

namespace sl {

  /// Utility class for initializing indexed containers
  template <class  G_indexed_container>
  class indexed_manifest_initializer {
  public:
    typedef indexed_manifest_initializer<G_indexed_container> this_t;
    typedef G_indexed_container                               container_t;
    typedef typename container_t::value_t                     value_t;
    typedef typename container_t::subscript_t                 subscript_t;
  protected:
    container_t& a_;
    subscript_t  i_;
  public:

    inline indexed_manifest_initializer(container_t& a, value_t x)
      :
      a_(a), i_(subscript_t()) 
    {
      SL_REQUIRE("Good size", a.extent().element_count() > 0);
      a_.put(x, i_);
    };

    inline ~indexed_manifest_initializer() {
      SL_REQUIRE("Matching sizes", 
		 i_.element_count() == 0 || 
		 i_.extent_from_last() == a_.extent());
      if (i_.element_count() == 0) { 
	// Fill with constant
	a_ = a_[subscript_t()];
      }
    }
  
    inline this_t& operator, (value_t x) {
      i_.increment(a_.extent());
      SL_REQUIRE("Good index", a_.good_subscript(i_));
      a_.put(x, i_);
      return *this;
    }

  }; 

} // namespace sl 

#endif

