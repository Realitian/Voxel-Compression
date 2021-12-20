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
#ifndef SL_TRIANGLE_MESH_SMOOTHER_HPP
#define SL_TRIANGLE_MESH_SMOOTHER_HPP

#include <sl/triangle_mesh.hpp>
#include <sl/math.hpp>
#include <sl/keyed_heap.hpp>
#include <set>

namespace sl {

  /**
   *  A simple triangle mesh smoother based on Taubin's lambda-mu
   *  technique. See:
   *
   *  G. Taubin, Geometric Signal Processing on Polygonal Meshes,  
   *  Eurographics 2000 STAR.
   */
  template <class T_Mesh>
  class triangle_mesh_taubin_smoother {
  public:
    typedef triangle_mesh_edge_collapser<T_Mesh> this_t;
    typedef T_Mesh   mesh_t;
    
    typedef typename mesh_t::vertex_t               vertex_t;
    typedef typename mesh_t::edge_t                 edge_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef typename mesh_t::small_edge_set_t       small_edge_set_t;
    typedef typename mesh_t::vertex_data_t::point_t point_t;
    typedef typename point_t::vector_t              vector_t;
    
    enum {dimension = point_t::dimension};
  protected:
    bool is_non_boundary_smoothing_enabled_;
    bool is_boundary_smoothing_enabled_;
  public:

    triangle_mesh_taubin_smoother() {
      is_non_boundary_smoothing_enabled_ = true;
      is_boundary_smoothing_enabled_ = true;
    }

    ~triangle_mesh_taubin_smoother() {
    }

    bool is_non_boundary_smoothing_enabled() const {
      return is_non_boundary_smoothing_enabled_;
    }

    void set_non_boundary_smoothing_enabled(bool x) {
      is_non_boundary_smoothing_enabled_ = x;
    }

    bool is_boundary_smoothing_enabled() const {
      return is_boundary_smoothing_enabled_;
    }

    void set_boundary_smoothing_enabled(bool x) {
      is_boundary_smoothing_enabled_ = x;
    }

    bool is_movable(const mesh_t& M, const vertex_t& vidx) const {
      bool result = is_boundary_smoothing_enabled() && is_non_boundary_smoothing_enabled();
      if (!result) {
        bool b = M.is_boundary_vertex(vidx);
        result =
          (b && is_boundary_smoothing_enabled()) ||
          (!b && is_non_boundary_smoothing_enabled());
      }
      return result;
    }

    vector_t laplacian(const mesh_t& M, vertex_t vi) const {
      vector_t result;
      const small_edge_set_t* incident_edges = M.vertex_edges(vi);
      if (incident_edges) {
        value_t  wi = sl::scalar_math<value_t>::zero();
        const point_t& xi = M.vertex_data(vi)->position();
        for (typename small_edge_set_t::const_iterator it = incident_edges->begin();
             it != incident_edges->end();
             ++it) {
          vertex_t  vj = (*it).opposite(vi);
          const point_t& xj = M.vertex_data(vj)->position();
          vector_t delta_ij = xj-xi;
          value_t  c_ij= 1.0f; // Equal weights -- delta_ij.two_norm();
          if (c_ij) {
            value_t w_ij=reciprocal(c_ij);
            result += w_ij * delta_ij;
            wi += w_ij;
          }
        }
        if (wi) result/=wi;
      }
      return result;
    }
    
    void iterative_smooth(mesh_t& M, std::size_t N, value_t lambda, value_t mu) {
      SL_TRACE_OUT(1) << "lambda = " << lambda << " - " << " mu = " << mu << std::endl;
      
      assert(lambda >= value_t(0.0));
      assert(mu <= -lambda);
      for (std::size_t i=0; i<2*N; ++i) {
        std::map<vertex_t,vector_t> delta_x;
        for (typename mesh_t::const_vertex_iterator_t v_it = M.vertex_begin();
             v_it != M.vertex_end();
             ++v_it) {
          delta_x[v_it->first] = laplacian(M,v_it->first);
        }
        value_t weight = (i%2==0) ? lambda : mu;
        for (typename mesh_t::const_vertex_iterator_t v_it = M.vertex_begin();
             v_it != M.vertex_end();
             ++v_it) {
          M.vertex_data(v_it->first)->position() += weight * delta_x[v_it->first];
        }
      }
    }

    void iterative_smooth(mesh_t& M, std::size_t N=50, value_t k_pb = value_t(0.1)) {
      assert(k_pb > 0.0);
      const value_t mu = (-k_pb - std::sqrt(k_pb*k_pb+4*(5-3*k_pb)))/(2*(5-3*k_pb));
      const value_t lambda = mu/(mu*k_pb-1);
      iterative_smooth(M,N,lambda,mu);
    }

  };

  template <class mesh_t>
  void triangle_mesh_smooth(mesh_t& M, std::size_t N=50, const typename mesh_t::vertex_data_t::value_t k_pb=0.1) {
    typedef triangle_mesh_taubin_smoother<mesh_t> smoother_t;
    smoother_t smoother;
    smoother.iterative_smooth(M, N, k_pb);
  }
  
  template <class mesh_t>
  void triangle_mesh_smooth_non_boundary_vertices(mesh_t& M, std::size_t N=50, const typename mesh_t::vertex_data_t::value_t k_pb=0.1) {
    typedef triangle_mesh_taubin_smoother<mesh_t> smoother_t;
    smoother_t smoother;
    smoother.set_boundary_smoothing_enabled(false);
    smoother.iterative_smooth(M, N, k_pb);
  }

  template <class mesh_t>
  void triangle_mesh_smooth_boundary_vertices(mesh_t& M, std::size_t N=50, const typename mesh_t::vertex_data_t::value_t k_pb=0.1) {
    typedef triangle_mesh_taubin_smoother<mesh_t> smoother_t;
    smoother_t smoother;
    smoother.set_non_boundary_smoothing_enabled(false);
    smoother.iterative_smooth(M, N, k_pb);
  }
 
}

#endif
