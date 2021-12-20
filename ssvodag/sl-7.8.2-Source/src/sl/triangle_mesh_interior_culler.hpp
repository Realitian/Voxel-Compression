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
#ifndef SL_TRIANGLE_MESH_INTERIOR_CULLER_HPP
#define SL_TRIANGLE_MESH_INTERIOR_CULLER_HPP

#include <sl/uniform_grid_container.hpp>
#include <sl/triangle_mesh.hpp>
#include <sl/triangle_mesh_utilities.hpp>

namespace sl {

  template <class T_Mesh>
  class triangle_mesh_interior_culler {
  public:
    typedef triangle_mesh_interior_culler<T_Mesh> this_t;
    typedef T_Mesh   mesh_t;
    
    typedef typename mesh_t::vertex_t               vertex_t;
    typedef typename mesh_t::edge_t                 edge_t;
    typedef typename mesh_t::triangle_t             triangle_t;
    typedef typename mesh_t::vertex_data_t::point_t point_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    enum {dimension = point_t::dimension };

    typedef typename mesh_t::const_triangle_iterator_t const_triangle_iterator_t;
    typedef typename mesh_t::const_vertex_iterator_t const_vertex_iterator_t;
    typedef typename mesh_t::const_edge_iterator_t const_edge_iterator_t;

  protected:
    
    typedef fixed_size_vector<sl::column_orientation,dimension,value_t> vector_t;
    typedef axis_aligned_box<dimension,value_t>                         aabox_t;

    typedef uniform_grid_container<dimension, value_t, const_triangle_iterator_t> ugrid_t;

    typedef typename ugrid_t::subscript_t subscript_t;
    
  public:

    triangle_mesh_interior_culler() {
    }

    ~triangle_mesh_interior_culler() {
    }

    /**
     *  Remove from the given mesh all portions guaranteed to be
     *  invisible from the outside, closing gaps larger than eps.
     */
    void cull_interior(mesh_t& M,
                       value_t eps) {
      SL_TRACE_OUT(1) << "Before interior culling: " << M.triangle_count() << " triangles." << std::endl;

      aabox_t bbox = triangle_mesh_aabox(M);
      value_t area = triangle_mesh_area(M);
      if ( area < value_t(1e-6) && M.triangle_count() > 0) {
	SL_TRACE_OUT( -1 ) << "computed mesh area=" << area << " (tcount=" << M.triangle_count() << " changed to " << 1e-6 << std::endl;
	area = value_t(1e-6);	
      }
      value_t cell_sz = std::max(eps,
                                 value_t(std::sqrt(2.0f/sqrt(3.0f)*area/sl::max(M.triangle_count(), std::size_t(1)))));
      
      vector_t hsl = 1.01f*bbox.half_side_lengths();
      value_t hsl_max = hsl[hsl.iamax()];
      value_t hsl_eps = std::max(value_t(1e-5), value_t(hsl_max * 0.0001));
      for (std::size_t i=0; i<dimension; ++i) {
        hsl[i] = std::max(value_t(hsl[i]+4.0*cell_sz), value_t(hsl_eps+4.0*cell_sz));
      }
      bbox = aabox_t(bbox.center()-hsl, bbox.center()+hsl);
      
      ugrid_t         grid;
      grid.clear();
      grid.resize(bbox,cell_sz);
      
      SL_TRACE_OUT(1) << "Grid: " << grid.extent() << std::endl;
      SL_TRACE_OUT(1) << "VOX : " << grid.voxel_dimensions()[0] << " " << grid.voxel_dimensions()[1] << " " << grid.voxel_dimensions()[2] << std::endl;

      // Engrid
      for (const_triangle_iterator_t t_it = M.triangle_map().begin();
           t_it != M.triangle_map().end();
           ++t_it) {
        const point_t& p0 = M.vertex_data(t_it->first[0])->position();
        const point_t& p1 = M.vertex_data(t_it->first[1])->position();
        const point_t& p2 = M.vertex_data(t_it->first[2])->position();

        engrid(grid, t_it, p0, p1, p2);
      }

      // Mark visible only triangles reachable from the outside. 
      subscript_t seed = exterior_seed(grid);
      if (!grid[seed].empty()) {
        SL_TRACE_OUT(-1) << "Cannot find seed - should not occur!" << std::endl;
      } else {
        std::set<subscript_t>    visited_cells;
        std::vector<subscript_t> candidate_cells;
        std::set<triangle_t>     invisible_triangles;
        M.triangles_in(invisible_triangles);
        
        candidate_cells.push_back(seed);
        while (!candidate_cells.empty()) {
          subscript_t idx = candidate_cells.back();
          candidate_cells.pop_back();

          if (visited_cells.find(idx) == visited_cells.end()) {
            //std::cerr << "Visit: " << idx[0] << " " << idx[1] << " " << idx[2] << std::endl;
            visited_cells.insert(idx);
            const typename ugrid_t::object_list_t& tlist = grid[idx];
            if (tlist.empty()) {
              for (int di = -1; di<= 1; di+=1) {
                for (int dj = -1; dj<= 1; dj+=1) {
                  for (int dk = -1; dk<= 1; dk+=1) {
                    subscript_t idx2;
                    idx2[0] = idx[0] + di;
                    idx2[1] = idx[1] + dj;
                    idx2[2] = idx[2] + dk;
                    if (grid.good_subscript(idx2) && (visited_cells.find(idx2) == visited_cells.end())) {
                      candidate_cells.push_back(idx2);
                    }
                  }
                }
              }
            } else {
              // Hit boundary!
              for (typename ugrid_t::object_list_t::const_iterator tit = tlist.begin();
                   tit != tlist.end();
                   ++tit) {
                //std::cerr << "Erase: " << (*tit)->first[0] << " " << (*tit)->first[1] << " " << (*tit)->first[2] << std::endl;
                invisible_triangles.erase((*tit)->first);
              }
            }
          }
        }

        // Erase all triangles unreachable from outside the box
        for (typename std::set<triangle_t>::iterator it = invisible_triangles.begin();
             it != invisible_triangles.end();
             ++it) {
          M.erase_triangle(*it);
        }
      }

      SL_TRACE_OUT(1) << "After interior culling: " << M.triangle_count() << " triangles." << std::endl;
    }


  protected:

    subscript_t exterior_seed(const ugrid_t& /*g*/) const {
      subscript_t result = subscript_t();
      // FIXME CHECK CHECK CHECK
      return result;
    }
        
    void engrid(ugrid_t& g,
                const_triangle_iterator_t t_it,
                const point_t& p0,
                const point_t& p1,
                const point_t& p2) {
      value_t eps = 0.5f*g.voxel_dimensions().two_norm_squared();
      value_t l01 = p0.distance_squared_to(p1);
      value_t l12 = p1.distance_squared_to(p2);
      value_t l20 = p2.distance_squared_to(p0);
      if (l01>l12) {
        if (l01>l20) {
          if (l01>eps) {
            point_t pmid = p0.lerp(p1,0.5f);
            engrid(g, t_it, p0, pmid, p1);
            engrid(g, t_it, pmid, p1, p2);
          } else {
            aabox_t t_box;
            t_box.to(p0);
            t_box.merge(p1);
            t_box.merge(p2);
            g.insert(t_box, t_it);
          }
        } else if (l20>eps) {
          point_t pmid = p2.lerp(p0,0.5f);
          engrid(g, t_it, p2, pmid, p1);
          engrid(g, t_it, pmid, p0, p1);
        } else {
          aabox_t t_box;
          t_box.to(p0);
          t_box.merge(p1);
          t_box.merge(p2);
          g.insert(t_box, t_it);
         }
      } else if (l12>l20) {
        if (l12>eps) {
          point_t pmid = p0.lerp(p1,0.5f);
          engrid(g, t_it, p1, pmid, p0);
          engrid(g, t_it, pmid, p2, p0);
        } else {
          aabox_t t_box;
          t_box.to(p0);
          t_box.merge(p1);
          t_box.merge(p2);
          g.insert(t_box, t_it);
        }
      } else if (l20>eps) {
        point_t pmid = p2.lerp(p0,0.5f);
        engrid(g, t_it, p2, pmid, p1);
        engrid(g, t_it, pmid, p0, p1);
      } else {
        aabox_t t_box;
        t_box.to(p0);
        t_box.merge(p1);
        t_box.merge(p2);
        g.insert(t_box, t_it);
      }        
    }

  };

  /**
   *  Remove from the given mesh all portions guaranteed to be
   *  invisible from the outside, closing gaps larger than eps.
   */
  template <class mesh_t>
  void triangle_mesh_cull_interior(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& eps) {
    typedef triangle_mesh_interior_culler<mesh_t> interior_culler_t;
    interior_culler_t interior_culler;
    interior_culler.cull_interior(M, eps);
  }

  /**
   *  Remove from the given mesh all portions guaranteed to be
   *  invisible from the outside, closing gaps larger than eps.
   */
  template <class mesh_t>
  void triangle_mesh_cull_interior(mesh_t& M) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef triangle_mesh_interior_culler<mesh_t> interior_culler_t;
    interior_culler_t interior_culler;
    interior_culler.cull_interior(M, value_t(0));
  }
    

} // namespace sl
    
#endif
