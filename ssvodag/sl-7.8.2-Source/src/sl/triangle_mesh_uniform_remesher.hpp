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
#ifndef SL_TRIANGLE_MESH_UNIFORM_REMESHER_HPP
#define SL_TRIANGLE_MESH_UNIFORM_REMESHER_HPP

#include <sl/triangle_mesh_refiner.hpp>
#include <sl/triangle_mesh_edge_collapser.hpp>
#include <sl/triangle_mesh_smoother.hpp>
#include <sl/triangle_mesh_utilities.hpp>

namespace sl {

  /// FIXME Just a first shot - needs better smoothing + valence optimization.
  template <class mesh_t>
  void triangle_mesh_remesh_to_target_edge_length(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    const value_t sl = value_t(4.0) * target_length / value_t(5.0);
    const value_t rl = value_t(4.0) * target_length / value_t(3.0); 

      SL_TRACE_OUT(1) << " - Initial" <<
        " - vcount = " << M.vertex_count() <<
        " - min = " << triangle_mesh_minimum_edge_length(M) <<
        " - avg = " << triangle_mesh_average_edge_length(M) <<
        " - max = " << triangle_mesh_maximum_edge_length(M) <<
        std::endl;
    // FIXME!!!
    for (std::size_t i=0; i<20; ++i) {
      triangle_mesh_refine_to_target_edge_len(M, rl);
      SL_TRACE_OUT(1) << i << " - after refine(" << rl << ")" <<
        " - vcount = " << M.vertex_count() <<
        " - min = " << triangle_mesh_minimum_edge_length(M) <<
        " - avg = " << triangle_mesh_average_edge_length(M) <<
        " - max = " << triangle_mesh_maximum_edge_length(M) <<
        std::endl;
      triangle_mesh_collapse_edges_to_target_edge_length(M,sl);
      SL_TRACE_OUT(1) << i << " - after collapse(" << sl << ")" <<
        " - vcount = " << M.vertex_count() <<
        " - min = " << triangle_mesh_minimum_edge_length(M) <<
        " - avg = " << triangle_mesh_average_edge_length(M) <<
        " - max = " << triangle_mesh_maximum_edge_length(M) <<
        std::endl;
      triangle_mesh_smooth(M,4);
      SL_TRACE_OUT(1) << i << " - after smoothing(" << 10 << ")" <<
        " - vcount = " << M.vertex_count() <<
        " - min = " << triangle_mesh_minimum_edge_length(M) <<
        " - avg = " << triangle_mesh_average_edge_length(M) <<
        " - max = " << triangle_mesh_maximum_edge_length(M) <<
        std::endl;
    }
  }
  
}

#endif
