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
#ifndef SL_TRIANGLE_MESH_REFINER_HPP
#define SL_TRIANGLE_MESH_REFINER_HPP

#include <sl/triangle_mesh.hpp>
#include <set>

namespace sl {

  template <class T_Mesh>
  class triangle_mesh_refiner {
  public:
    typedef triangle_mesh_refiner<T_Mesh> this_t;
    typedef T_Mesh   mesh_t;
    
    typedef typename mesh_t::vertex_t               vertex_t;
    typedef typename mesh_t::edge_t                 edge_t;
    typedef typename mesh_t::vertex_data_t::point_t point_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    enum {dimension = point_t::dimension };

  protected:

    typedef std::pair<value_t, edge_t> candidate_t;
    typedef std::set<candidate_t>    candidate_set_t;

    bool is_non_boundary_edge_refinement_enabled_;
    bool is_boundary_edge_refinement_enabled_;
    bool is_single_pass_enabled_;

  public:

    triangle_mesh_refiner() {
      is_non_boundary_edge_refinement_enabled_ = true;
      is_boundary_edge_refinement_enabled_ = true;
      is_single_pass_enabled_ = true;
    }

    ~triangle_mesh_refiner() {
    }

    bool is_non_boundary_edge_refinement_enabled() const {
      return is_non_boundary_edge_refinement_enabled_;
    }

    void set_non_boundary_edge_refinement_enabled(bool x) {
      is_non_boundary_edge_refinement_enabled_ = x;
    }

    bool is_boundary_edge_refinement_enabled() const {
      return is_boundary_edge_refinement_enabled_;
    }

    void set_boundary_edge_refinement_enabled(bool x) {
      is_boundary_edge_refinement_enabled_ = x;
    }

    bool is_single_pass_enabled() const {
      return is_single_pass_enabled_;
    }

    void set_single_pass_enabled(bool x) {
      is_single_pass_enabled_ = x;
    }

    bool is_splittable(const mesh_t& M, const edge_t& e) const {
      bool result = is_boundary_edge_refinement_enabled() && is_non_boundary_edge_refinement_enabled();
      if (!result) {
        bool b = M.is_boundary_edge(e);
        result =
          (b && is_boundary_edge_refinement_enabled()) ||
          (!b && is_non_boundary_edge_refinement_enabled());
      }
      return result;
    }
 
    /**
     *  Refine mesh M by longest edge bisection. Refinement
     *  stops when no edge is longer than target_edge_len or
     *  the mesh has exceded the given maximum triangle count
     */
    void refine(mesh_t& M,
                value_t target_edge_len,
                std::size_t max_tricount) {
      candidate_set_t candidate_set;
      for (typename mesh_t::const_edge_iterator_t e_it = M.edge_map().begin();
           e_it != M.edge_map().end(); ++e_it) {
        edge_t  e = e_it->first;
        value_t l2 = len2(M, e);
        if (l2>target_edge_len*target_edge_len && is_splittable(M,e)) {
          candidate_set.insert(std::make_pair(-l2, e));
        }
      }
      while (!candidate_set.empty() && M.triangle_count() < max_tricount) {
        candidate_t c = *(candidate_set.begin());
        candidate_set.erase(candidate_set.begin());
        edge_t                e = c.second;

        // Refine
        vertex_t              vmid = M.new_vertex_id();
        typename mesh_t::vertex_data_t vmid_data = (*M.vertex_data(e[0])).lerp((*M.vertex_data(e[1])), 0.5f);
        M.split_edge(e, vmid, vmid_data);

	if (!is_single_pass_enabled_) {
	  // Insert new edges
	  const typename mesh_t::small_edge_set_t* e_vmid = M.vertex_edges(vmid);
	  for (typename mesh_t::small_edge_set_t::const_iterator e_it = e_vmid->begin();
	       e_it != e_vmid->end();
	       ++e_it) {
	    edge_t  e = *e_it;
	    value_t l2 = len2(M, e);
	    if (l2>target_edge_len*target_edge_len && is_splittable(M,e)) {
	      candidate_set.insert(std::make_pair(-l2, e));
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

  // Multipass

  template <class mesh_t>
  void triangle_mesh_refine(mesh_t& M, typename mesh_t::vertex_data_t::value_t& target_edge_len, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.refine(M, target_edge_len, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_to_target_tricount(mesh_t& M, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.refine(M, 0.0f, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_to_target_edge_len(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& target_edge_len) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.refine(M, target_edge_len, std::size_t(-1));
  }
  
  template <class mesh_t>
  void triangle_mesh_refine_boundary(mesh_t& M, typename mesh_t::vertex_data_t::value_t& target_edge_len, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_non_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_boundary_to_target_tricount(mesh_t& M, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_non_boundary_edge_refinement_enabled(false);
    refiner.refine(M, 0.0f, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_boundary_to_target_edge_len(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& target_edge_len) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_non_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, std::size_t(-1));
  }
  
  template <class mesh_t>
  void triangle_mesh_refine_non_boundary(mesh_t& M, typename mesh_t::vertex_data_t::value_t& target_edge_len, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_non_boundary_to_target_tricount(mesh_t& M, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_boundary_edge_refinement_enabled(false);
    refiner.refine(M, 0.0f, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_non_boundary_to_target_edge_len(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& target_edge_len) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, std::size_t(-1));
  }
   
  // Single passs

  template <class mesh_t>
  void triangle_mesh_refine_single_pass(mesh_t& M, typename mesh_t::vertex_data_t::value_t& target_edge_len, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.refine(M, target_edge_len, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_to_target_tricount_single_pass(mesh_t& M, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.refine(M, 0.0f, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_to_target_edge_len_single_pass(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& target_edge_len) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.refine(M, target_edge_len, std::size_t(-1));
  }
  
  template <class mesh_t>
  void triangle_mesh_refine_boundary_single_pass(mesh_t& M, typename mesh_t::vertex_data_t::value_t& target_edge_len, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.set_non_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_boundary_to_target_tricount_single_pass(mesh_t& M, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.set_non_boundary_edge_refinement_enabled(false);
    refiner.refine(M, 0.0f, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_boundary_to_target_edge_len_single_pass(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& target_edge_len) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.set_non_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, std::size_t(-1));
  }
  
  template <class mesh_t>
  void triangle_mesh_refine_non_boundary_single_pass(mesh_t& M, typename mesh_t::vertex_data_t::value_t& target_edge_len, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.set_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_non_boundary_to_target_tricount_single_pass(mesh_t& M, std::size_t max_tricount) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
    refiner.set_single_pass_enabled(true);
    refiner.set_boundary_edge_refinement_enabled(false);
    refiner.refine(M, 0.0f, max_tricount);
  }

  template <class mesh_t>
  void triangle_mesh_refine_non_boundary_to_target_edge_len_single_pass(mesh_t& M, const typename mesh_t::vertex_data_t::value_t& target_edge_len) {
    typedef triangle_mesh_refiner<mesh_t> refiner_t;
    refiner_t refiner;
     refiner.set_single_pass_enabled(true);
     refiner.set_boundary_edge_refinement_enabled(false);
    refiner.refine(M, target_edge_len, std::size_t(-1));
  }

} // namespace sl
    
#endif
