//+++HDR+++
//======================================================================
//  This file is part of the SL software library.
//
//  Copyright (C) 1993-2010 by Enrico Gobbetti (gobbetti@crs4.it)
//  Copyright (C) 1996-2010 by CRS4 Visual Computing Group, Pula, Italy
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
#ifndef SL_TRIANGLE_MESH_EDGE_COLLAPSER_HPP
#define SL_TRIANGLE_MESH_EDGE_COLLAPSER_HPP

#include <sl/triangle_mesh.hpp>
#include <sl/math.hpp>
#include <sl/keyed_heap.hpp>
#include <set>

namespace sl {

  /// A simple purely distance based edge collapser
  /// FIXME - very rough/stupid implementation!
  template <class T_Mesh>
  class triangle_mesh_edge_collapser {
  public:
    typedef triangle_mesh_edge_collapser<T_Mesh> this_t;
    typedef T_Mesh   mesh_t;
    
    typedef typename mesh_t::vertex_t               vertex_t;
    typedef typename mesh_t::edge_t                 edge_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    typedef typename mesh_t::small_edge_set_t       small_edge_set_t;
    enum {dimension = mesh_t::vertex_data_t::point_t::dimension};
  protected:

    typedef std::pair<value_t, edge_t>  candidate_t;
    typedef keyed_heap<edge_t, value_t> candidate_set_t;
    
    bool is_non_boundary_edge_collapse_enabled_;
    bool is_boundary_edge_collapse_enabled_;
    bool is_single_pass_enabled_;
  public:

    triangle_mesh_edge_collapser() {
      is_non_boundary_edge_collapse_enabled_ = true;
      is_boundary_edge_collapse_enabled_ = true;
      is_single_pass_enabled_ = false;
    }

    ~triangle_mesh_edge_collapser() {
    }

    bool is_non_boundary_edge_collapse_enabled() const {
      return is_non_boundary_edge_collapse_enabled_;
    }

    void set_non_boundary_edge_collapse_enabled(bool x) {
      is_non_boundary_edge_collapse_enabled_ = x;
    }

    bool is_boundary_edge_collapse_enabled() const {
      return is_boundary_edge_collapse_enabled_;
    }

    void set_boundary_edge_collapse_enabled(bool x) {
      is_boundary_edge_collapse_enabled_ = x;
    }

    bool is_single_pass_enabled() const {
      return is_single_pass_enabled_;
    }

    void set_single_pass_enabled(bool x) {
      is_single_pass_enabled_ = x;
    }

    bool is_collapsable(const mesh_t& M, const edge_t& e) const {
      bool result = is_boundary_edge_collapse_enabled() && is_non_boundary_edge_collapse_enabled();
      if (!result) {
        bool b = M.is_boundary_edge(e);
        result =
          (b && is_boundary_edge_collapse_enabled()) ||
          (!b && is_non_boundary_edge_collapse_enabled());

        if (result && !is_boundary_edge_collapse_enabled()) {
          // Check if erasing the interior edge would remove a boundary vertex
          const typename mesh_t::small_triangle_set_t* e_tri_ptr = M.edge_triangles(e);
          for (typename mesh_t::small_triangle_set_t::const_iterator t_it = e_tri_ptr->begin();
               t_it != e_tri_ptr->end() && result;
               ++t_it) {
            typename mesh_t::vertex_t v = t_it->opposite(e);
            result = M.vertex_triangle_count(v) > 1;
          }
        }
      }
      return result;
    }
 
    void collapse_edges(mesh_t& M,
                        value_t target_edge_length,
			value_t max_new_edge_length,
                        std::size_t min_tricount) {
      const value_t target_l2 = target_edge_length*target_edge_length;
      const value_t max_new_l2 = max_new_edge_length*max_new_edge_length;
      const value_t eps_new_l2 = max_new_l2/(100.0f*100.0f);
     
      // insert all candidate edges into heap
      candidate_set_t candidate_set;
      for (typename mesh_t::const_edge_iterator_t e_it = M.edge_begin();
           e_it != M.edge_end(); ++e_it) {
        edge_t  e  = e_it->first;
        value_t l2 = len2(M, e);
        if ( l2 < target_l2 && is_collapsable(M,e) ) {
          candidate_set.push(std::make_pair(-l2, e));
        }
      }
        
      // collapse shortest edges
      while (!candidate_set.empty() && M.triangle_count() > min_tricount) {
        candidate_t c = candidate_set.top();
        candidate_set.pop();
        value_t l2 = -c.first;
	edge_t  e = c.second;
        
        // simplify
        if (M.has_edge(e)) {

          bool    collapse_feasible = false;
          value_t t = 0.5f;
          vertex_t src_idx = e[0];
          vertex_t dst_idx = e[1];
          bool     e0_boundary = M.is_boundary_vertex(e[0]);
          bool     e1_boundary = M.is_boundary_vertex(e[1]);

          // FIXME: Improve the following:
          //   * check corner vertices
          //   * check normal flipping
          //   * ...
          
          if (e0_boundary == e1_boundary) {
            collapse_feasible = is_collapsable(M,e);
            t = 0.5f;
          } else if (e0_boundary) {
            collapse_feasible = is_collapsable(M,e);
            src_idx = e[1];
            dst_idx = e[0];
            t = 1.0f;
          } else if (e1_boundary) {
            collapse_feasible = is_collapsable(M,e);
            src_idx = e[0];
            dst_idx = e[1];
            t = 1.0f;
          }

          if (collapse_feasible) {
            typename mesh_t::vertex_data_t vmid_data = (*M.vertex_data(src_idx)).lerp((*M.vertex_data(dst_idx)), t);
	    
	    const small_edge_set_t& incident_edges_src = *(M.vertex_edges(src_idx));
	    const small_edge_set_t& incident_edges_dst = *(M.vertex_edges(dst_idx));
	      
	    if (l2>eps_new_l2) {
	      // Check if we are creating too long edges if we are collapsing a non-null edge
	      for (typename small_edge_set_t::const_iterator it = incident_edges_dst.begin();
		   it != incident_edges_dst.end() && collapse_feasible;
		   ++it) {
		vertex_t  v = it->opposite(dst_idx);
		collapse_feasible = (M.vertex_data(v)->position().distance_squared_to(vmid_data.position()) <= max_new_l2);
	      }
	      for (typename small_edge_set_t::const_iterator it = incident_edges_src.begin();
		   it != incident_edges_src.end() && collapse_feasible;
		   ++it) {
		vertex_t  v = it->opposite(src_idx);
		collapse_feasible = (M.vertex_data(v)->position().distance_squared_to(vmid_data.position()) <= max_new_l2);
	      }
	    }

	    if (collapse_feasible) {
	      for (typename small_edge_set_t::const_iterator it = incident_edges_dst.begin();
		   it != incident_edges_dst.end();
		   ++it) {
		edge_t  e = *it;
		candidate_set.erase(e);
	      }
	      

	      M.collapse_edge(src_idx, dst_idx);
	      
	      if (M.has_vertex(dst_idx)) {
		*(M.vertex_data(dst_idx)) = vmid_data; // move to midpoint
		
		if (!is_single_pass_enabled_) {
		  const small_edge_set_t& incident_edges = *(M.vertex_edges(dst_idx));
		  for (typename small_edge_set_t::const_iterator it = incident_edges.begin();
		       it != incident_edges.end();
		       ++it) {
		    edge_t  e = *it;
		    value_t l2 = len2(M, e);
		    if (l2 < target_l2 && is_collapsable(M,e)) {
		    candidate_set.push(std::make_pair(-l2, e));
		    }
		  }
		}
	      }
	    }
          }
        }
      }
    }

  protected:
    
    value_t len2(const mesh_t& M,
                 const edge_t& e) const {
      return M.vertex_data(e[0])->position().distance_squared_to(M.vertex_data(e[1])->position());
    }

  };

  template <class mesh_t>
  void triangle_mesh_collapse_edges_to_target_tricount(mesh_t& M, std::size_t target_tricount) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.collapse_edges(M, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), target_tricount);
  }
 
  template <class mesh_t>
  void triangle_mesh_collapse_edges_to_target_edge_length(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length, const typename mesh_t::vertex_data_t::value_t max_new_length=scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound()) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.collapse_edges(M, target_length, max_new_length, std::size_t(0));
  }
  
  template <class mesh_t>
  void triangle_mesh_collapse_boundary_edges_to_target_tricount(mesh_t& M, std::size_t target_tricount) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_non_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), target_tricount);
  }
 
  template <class mesh_t>
  void triangle_mesh_collapse_boundary_edges_to_target_edge_length(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length, const typename mesh_t::vertex_data_t::value_t max_new_length=scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound()) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_non_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, target_length, max_new_length, std::size_t(0));
  }

  template <class mesh_t>
  void triangle_mesh_collapse_non_boundary_edges_to_target_tricount(mesh_t& M, std::size_t target_tricount) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), target_tricount);
  }
 
  template <class mesh_t>
  void triangle_mesh_collapse_non_boundary_edges_to_target_edge_length(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length, const typename mesh_t::vertex_data_t::value_t max_new_length=scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound()) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, target_length, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), std::size_t(0));
  }

  // Single pass version
  template <class mesh_t>
  void triangle_mesh_collapse_edges_to_target_tricount_single_pass(mesh_t& M, std::size_t target_tricount) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_single_pass_enabled(true);
    edge_collapser.collapse_edges(M, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), target_tricount);
  }
 
  template <class mesh_t>
  void triangle_mesh_collapse_edges_to_target_edge_length_single_pass(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length, const typename mesh_t::vertex_data_t::value_t max_new_length=scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound()) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_single_pass_enabled(true);
    edge_collapser.collapse_edges(M, target_length, max_new_length, std::size_t(0));
  }
  
  template <class mesh_t>
  void triangle_mesh_collapse_boundary_edges_to_target_tricount_single_pass(mesh_t& M, std::size_t target_tricount) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_single_pass_enabled(true);
    edge_collapser.set_non_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), target_tricount);
  }
 
  template <class mesh_t>
  void triangle_mesh_collapse_boundary_edges_to_target_edge_length_single_pass(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length, const typename mesh_t::vertex_data_t::value_t max_new_length=scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound()) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_single_pass_enabled(true);
    edge_collapser.set_non_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, target_length, max_new_length, std::size_t(0));
  }

  template <class mesh_t>
  void triangle_mesh_collapse_non_boundary_edges_to_target_tricount_single_pass(mesh_t& M, std::size_t target_tricount) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_single_pass_enabled(true);
    edge_collapser.set_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), target_tricount);
  }
 
  template <class mesh_t>
  void triangle_mesh_collapse_non_boundary_edges_to_target_edge_length_single_pass(mesh_t& M, const typename mesh_t::vertex_data_t::value_t target_length, const typename mesh_t::vertex_data_t::value_t max_new_length=scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound()) {
    typedef triangle_mesh_edge_collapser<mesh_t> edge_collapser_t;
    edge_collapser_t edge_collapser;
    edge_collapser.set_single_pass_enabled(true);
    edge_collapser.set_boundary_edge_collapse_enabled(false);
    edge_collapser.collapse_edges(M, target_length, scalar_math<typename mesh_t::vertex_data_t::value_t>::finite_upper_bound(), std::size_t(0));
  }

 
}

#endif
