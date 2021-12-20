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
#ifndef SL_TRIANGLE_MESH_VERTEX_CLUSTERER_HPP
#define SL_TRIANGLE_MESH_VERTEX_CLUSTERER_HPP

#include <sl/uniform_grid_container.hpp>
#include <sl/triangle_mesh.hpp>
#include <sl/triangle_mesh_utilities.hpp>

namespace sl {

  template <class T_Mesh>
  class triangle_mesh_vertex_clusterer {
  public:
    typedef triangle_mesh_vertex_clusterer<T_Mesh> this_t;
    typedef T_Mesh   mesh_t;
    
    typedef typename mesh_t::vertex_t               vertex_t;
    typedef typename mesh_t::edge_t                 edge_t;
    typedef typename mesh_t::triangle_t             triangle_t;
    typedef typename mesh_t::vertex_data_t          vertex_data_t;
    typedef typename mesh_t::vertex_data_t::point_t point_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    enum {dimension = point_t::dimension };

    typedef typename mesh_t::const_triangle_iterator_t const_triangle_iterator_t;
    typedef typename mesh_t::const_vertex_iterator_t   const_vertex_iterator_t;
    typedef typename mesh_t::const_edge_iterator_t     const_edge_iterator_t;

  protected:
    
    typedef fixed_size_vector<sl::column_orientation,dimension,value_t> vector_t;
    typedef axis_aligned_box<dimension,value_t>                         aabox_t;

    typedef uniform_grid_container<dimension, value_t, vertex_t>        ugrid_t;

    typedef typename ugrid_t::subscript_t subscript_t;

    bool is_non_boundary_vertex_clustering_enabled_;
    bool is_boundary_vertex_clustering_enabled_;
  public:

    triangle_mesh_vertex_clusterer() {
      is_non_boundary_vertex_clustering_enabled_ = true;
      is_boundary_vertex_clustering_enabled_ = true;
    }

    ~triangle_mesh_vertex_clusterer() {
    }

    bool is_non_boundary_vertex_clustering_enabled() const {
      return is_non_boundary_vertex_clustering_enabled_;
    }

    void set_non_boundary_vertex_clustering_enabled(bool x) {
      is_non_boundary_vertex_clustering_enabled_ = x;
    }

    bool is_boundary_vertex_clustering_enabled() const {
      return is_boundary_vertex_clustering_enabled_;
    }

    void set_boundary_vertex_clustering_enabled(bool x) {
      is_boundary_vertex_clustering_enabled_ = x;
    }

    bool is_collapsable(const mesh_t& M, const vertex_t& vidx) const {
      bool result = is_boundary_vertex_clustering_enabled() && is_non_boundary_vertex_clustering_enabled();
      if (!result) {
       bool b = M.is_boundary_vertex(vidx);
       result =
         (b && is_boundary_vertex_clustering_enabled()) ||
         (!b && is_non_boundary_vertex_clustering_enabled());
      }
      return result;
    }
    
    /**
     *  Cluster all vertices that are within eps
     */
    void cluster_vertices(mesh_t& M,
                          value_t eps) {
      SL_TRACE_OUT(1) << "Before vertex clustering: " << M.triangle_count() << " triangles / " << M.vertex_count() << " vertices." << std::endl;

      value_t eps2 = std::max(eps*eps, value_t(6)*sl::scalar_math<value_t>::epsilon());
      aabox_t bbox = triangle_mesh_aabox(M);
      value_t area = triangle_mesh_area(M);
      if ( area < value_t(1e-6) && M.triangle_count() > 0) {
	SL_TRACE_OUT( -1 ) << "computed mesh area=" << area << " (tcount=" << M.triangle_count() << " changed to " << 1e-6 << std::endl;
	area = value_t(1e-6);	
      }

      value_t elen = value_t(std::sqrt(2.0f/sqrt(3.0f)*area/sl::max(M.triangle_count(), std::size_t(1))));
      value_t cell_sz = std::max(eps,value_t(0.01f*elen));
      
      vector_t hsl = 1.01f*bbox.half_side_lengths();
      value_t hsl_max = hsl[hsl.iamax()];
      value_t hsl_eps = std::max(value_t(1e-5), value_t(hsl_max * 0.0001));
      for (std::size_t i=0; i<dimension; ++i) {
        hsl[i] = std::max(value_t(hsl[i]+4.0*cell_sz), value_t(hsl_eps+4.0*cell_sz));
      }
      
      ugrid_t         grid;
      aabox_t bbox0 = aabox_t(bbox.center()-hsl, bbox.center()+hsl);
      
      // Two passes, with a grid shift of half a cell size
      for (std::size_t i=0; i<2; ++i) {
        aabox_t shifted_bbox = aabox_t(bbox0[0]+(i*0.5f)*vector_t(cell_sz, cell_sz, cell_sz),
                                       bbox0[1]+(i*0.5f)*vector_t(cell_sz, cell_sz, cell_sz));
        grid.clear();
        grid.resize(shifted_bbox,cell_sz);
      
        SL_TRACE_OUT(1) << "Grid : " << grid.extent() << std::endl;
        SL_TRACE_OUT(1) << "ELEN : " << elen << std::endl;
        SL_TRACE_OUT(1) << "VOX  : " << grid.voxel_dimensions()[0] << " " << grid.voxel_dimensions()[1] << " " << grid.voxel_dimensions()[2] << std::endl;
        
        // Engrid
        SL_TRACE_OUT(1) << "Pass: " << i << " Engridding..." << std::endl;
        for (const_vertex_iterator_t v_it = M.vertex_map().begin();
             v_it != M.vertex_map().end();
             ++v_it) {
          vertex_t vidx = v_it->first;
          if (is_collapsable(M,vidx)) {
            const point_t&  p   = v_it->second.data()->position();
            grid.insert(p, vidx);
          }
        }
        
        SL_TRACE_OUT(1) << "Pass: " << i << " Clustering..." << std::endl;
        // Cluster all vertices that lie in the same cell
        subscript_t l = subscript_t();
        subscript_t h = grid.extent();
        for (subscript_t idx = l; idx.all_le(h); idx.increment(l,h)) {
          if (grid.good_subscript(idx)) {
            typedef typename ugrid_t::object_list_t vidxlist_t;
            vidxlist_t& vidxlist = grid[idx];
            std::size_t vcount = vidxlist.size();
            if (vcount >= 2) {
              // Cluster vertices - quick and dirty a la Rossignac/Borrel
              typename vidxlist_t::const_iterator vit = vidxlist.begin();
              vertex_t vidx0 = *vit;
              ++vit;
              while (vit != vidxlist.end() &&
                     ((!M.has_vertex(vidx0)) || !is_collapsable(M,vidx0))) {
                vidx0 = *vit;
                ++vit;
              }
              while (vit != vidxlist.end()) {
                vertex_t vidx1 = *vit;
                if (M.has_vertex(vidx0) && (is_collapsable(M,vidx0))) {
                  if (M.has_vertex(vidx1) && (is_collapsable(M,vidx1))) {
                    if (!M.has_edge(edge_t(vidx0,vidx1)) &&
                        M.vertex_data(vidx0)->position().distance_squared_to(M.vertex_data(vidx1)->position()) <= eps2) {
                      SL_TRACE_OUT(1) << "   Collapse: " << vidx1 << " -> " << vidx0 << std::endl;
                      vertex_data_t vdata_avg = M.vertex_data(vidx0)->lerp(*M.vertex_data(vidx1), 0.5f);
                      *(M.vertex_data(vidx0)) = vdata_avg;
                      M.collapse_edge(vidx1, vidx0);
                    }
                  }
                } else {
                  vidx0 = vidx1;
                }
                ++vit;
              }
            }
          }
        }
      }

      SL_TRACE_OUT(1) << "After vertex clustering: " << M.triangle_count() << " triangles / " << M.vertex_count() << " vertices." << std::endl;
    }
  };

  template <class mesh_t>
  void triangle_mesh_cluster_vertices(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& eps) {
    typedef triangle_mesh_vertex_clusterer<mesh_t> vertex_clusterer_t;
    vertex_clusterer_t vertex_clusterer;
    vertex_clusterer.cluster_vertices(M, eps);
  }

  template <class mesh_t>
  void triangle_mesh_cluster_vertices(mesh_t& M) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef triangle_mesh_vertex_clusterer<mesh_t> vertex_clusterer_t;
    vertex_clusterer_t vertex_clusterer;
    vertex_clusterer.cluster_vertices(M, value_t(0));
  }

  template <class mesh_t>
  void triangle_mesh_cluster_boundary_vertices(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& eps) {
    typedef triangle_mesh_vertex_clusterer<mesh_t> vertex_clusterer_t;
    vertex_clusterer_t vertex_clusterer;
    vertex_clusterer.set_non_boundary_vertex_clustering_enabled(false);
    vertex_clusterer.cluster_vertices(M, eps);
  }

  template <class mesh_t>
  void triangle_mesh_cluster_boundary_vertices(mesh_t& M) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef triangle_mesh_vertex_clusterer<mesh_t> vertex_clusterer_t;
    vertex_clusterer_t vertex_clusterer;
    vertex_clusterer.set_non_boundary_vertex_clustering_enabled(false);
    vertex_clusterer.cluster_vertices(M, value_t(0));
  }

  template <class mesh_t>
  void triangle_mesh_cluster_non_boundary_vertices(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& eps) {
    typedef triangle_mesh_vertex_clusterer<mesh_t> vertex_clusterer_t;
    vertex_clusterer_t vertex_clusterer;
    vertex_clusterer.set_boundary_vertex_clustering_enabled(false);
    vertex_clusterer.cluster_vertices(M, eps);
  }

  template <class mesh_t>
  void triangle_mesh_cluster_non_boundary_vertices(mesh_t& M) {
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef triangle_mesh_vertex_clusterer<mesh_t> vertex_clusterer_t;
    vertex_clusterer_t vertex_clusterer;
    vertex_clusterer.set_boundary_vertex_clustering_enabled(false);
    vertex_clusterer.cluster_vertices(M, value_t(0));
  }

} // namespace sl
    
#endif
