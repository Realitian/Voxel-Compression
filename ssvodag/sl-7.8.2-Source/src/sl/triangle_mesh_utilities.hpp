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
#ifndef SL_TRIANGLE_MESH_UTILITIES_HPP
#define SL_TRIANGLE_MESH_UTILITIES_HPP

#include <sl/triangle_mesh.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/ball.hpp>

namespace sl {

  /// The average edge length of a triangle mesh. vertex_data_t must define point_t, value_t, position
  template <class mesh_t>
  typename mesh_t::vertex_data_t::value_t triangle_mesh_average_edge_length(const mesh_t& m) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;

    value_t result = sl::scalar_math<value_t>::zero();

    std::size_t N = m.edge_count();
    if (N) {
      for (typename mesh_t::const_edge_iterator_t e_it = m.edge_begin();
           e_it != m.edge_end();
           ++e_it) {
        result +=  m.vertex_data(e_it->first[0])->position().distance_to(m.vertex_data(e_it->first[1])->position());
      }
      result /= value_t(N);
    }
    return result;
  }
  
  /// The longest edge length of a triangle mesh. vertex_data_t must define point_t, value_t, position
  template <class mesh_t>
  typename mesh_t::vertex_data_t::value_t triangle_mesh_maximum_edge_length(const mesh_t& m) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;

    value_t result = sl::scalar_math<value_t>::zero();
    for (typename mesh_t::const_edge_iterator_t e_it = m.edge_begin();
         e_it != m.edge_end();
         ++e_it) {
      result = std::max(result,m.vertex_data(e_it->first[0])->position().distance_to(m.vertex_data(e_it->first[1])->position()));
    }
    return result;
  }

  /// The shortest edge length of a triangle mesh. vertex_data_t must define point_t, value_t, position
  template <class mesh_t>
  typename mesh_t::vertex_data_t::value_t triangle_mesh_minimum_edge_length(const mesh_t& m) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;

    value_t result = sl::scalar_math<value_t>::zero();
    if (m.edge_count() > 0) {
      result = sl::scalar_math<value_t>::finite_upper_bound();
      for (typename mesh_t::const_edge_iterator_t e_it = m.edge_begin();
           e_it != m.edge_end();
           ++e_it) {
        result = std::min(result,m.vertex_data(e_it->first[0])->position().distance_to(m.vertex_data(e_it->first[1])->position()));
      }
    }
    return result;
  }

  /// The area of a triangle mesh. vertex_data_t must define point_t, value_t, position
  template <class mesh_t>
  typename mesh_t::vertex_data_t::value_t triangle_mesh_area(const mesh_t& m) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;

    value_t result = sl::scalar_math<value_t>::zero();
    for (typename mesh_t::const_triangle_iterator_t t_it = m.triangle_begin();
         t_it != m.triangle_end();
         ++t_it) {
      result +=  std::sqrt(sl::triangle_area_squared(m.vertex_data(t_it->first[0])->position(),
                                                     m.vertex_data(t_it->first[1])->position(),
                                                     m.vertex_data(t_it->first[2])->position()));
    }

    return result;
  }
  
  /// The axis aligned box of a triangle mesh. vertex_data_t must define point_t, value_t, position
  template <class mesh_t>
  axis_aligned_box<mesh_t::vertex_data_t::point_t::dimension, typename mesh_t::vertex_data_t::point_t::value_t> triangle_mesh_aabox(const mesh_t& m) {
    typedef typename mesh_t::vertex_data_t::point_t point_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    enum {dimension = point_t::dimension };
    typedef axis_aligned_box<dimension,value_t> aabox_t;

    aabox_t result;
    result.to_empty();
    for (typename mesh_t::const_vertex_iterator_t v_it = m.vertex_begin();
         v_it!= m.vertex_end();
         ++v_it) {
      result.merge(v_it->second.data()->position());
    }
    return result;
  }

  /// The minumum enclosing ball of a triangle mesh. vertex_data_t must define point_t, value_t, position
  template <class mesh_t>
  ball< mesh_t::vertex_data_t::point_t::dimension,
               typename mesh_t::vertex_data_t::value_t > triangle_mesh_minimum_enclosing_ball(const mesh_t& M1) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef ball_builder<mesh_t::vertex_data_t::point_t::dimension, value_t> ball_builder_t;
    typedef typename mesh_t::const_vertex_iterator_t const_vertex_iterator_t;

    ball_builder_t mbb;
    mbb.begin_model();
    for (const_vertex_iterator_t v_it = M1.vertex_begin();
         v_it!= M1.vertex_end();
         ++v_it) {
      mbb.put_point(v_it->second.data()->position());
    }
    mbb.end_model();
    return mbb.last_bounding_volume();
  }

}

#endif
