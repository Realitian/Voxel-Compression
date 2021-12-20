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
#ifndef SL_TRIANGLE_MESH_STRIPIFIER_HPP
#define SL_TRIANGLE_MESH_STRIPIFIER_HPP

#include <sl/triangle_mesh.hpp>
#include <sl/fifo_cache_simulator.hpp>
#include <sl/stl_container_selector.hpp>
#include <list>

namespace sl {  

  // =============================================================================
  // triangle_mesh_stripifier
  // =============================================================================

  /**
   *  Objects that transform triangle lists into triangle strips
   */
  class triangle_mesh_stripifier {
  public:
    typedef std::vector<std::size_t>        strip_t;
  protected:
    std::size_t cache_size_;

    bool is_during_input_;
    bool is_during_output_;
    std::size_t output_strip_cursor_;
    std::size_t output_vertex_cursor_;
    
    std::vector< strip_t > tstrips_;

  public: // Creation

    triangle_mesh_stripifier();
    ~triangle_mesh_stripifier();
    
    // Set cache size: 0 will disable cache optimizer
    virtual void set_cache_size(std::size_t x = 32);

    inline std::size_t cache_size() const { return cache_size_; }
    
  public: // Input

    virtual void begin_input() = 0;
    bool is_during_input() const;
    virtual void insert_input_triangle(std::size_t i0,
				       std::size_t i1,
				       std::size_t i2) = 0;
    virtual void end_input() = 0;

  public: // Output

    std::size_t output_strip_count() const;
    void begin_output();
    bool is_during_output() const;
    bool has_output_strip() const;
    bool has_output_vertex() const;
    void begin_output_strip();
    std::size_t current_output_strip_index() const;
    std::size_t current_output_strip_vertex_index() const;
    std::size_t get_output_strip_vertex();
    void end_output_strip();
    void end_output();

  protected:

    void strip_grow_in(strip_t& strip,
		       std::size_t i0,
		       std::size_t i1,
		       std::size_t i2);
  };
  
  // =============================================================================
  // triangle_mesh_greedy_stripifier
  // =============================================================================

  namespace detail {
    struct tmg_vertex {
      std::size_t id_;
      int         lru_position_;
      float       score_;
      sl::vector_set<std::size_t>  active_incident_triangles_;

      inline tmg_vertex(std::size_t id): id_(id), lru_position_(-1), score_(-1.0f) {
      }

      inline std::size_t active_valence() const {
	return active_incident_triangles_.size();
      }
    };
    struct tmg_triangle {
      bool is_active_; // Still to be added to draw list
      triangle_connectivity connectivity_;
      float score_;
      
      inline tmg_triangle(): is_active_(false), score_(-3.0f) {
      }
    };
  } // namespace detail
      
  /**
   *  A fast, linear speed triangle mesh stripifier. Loosely based on
   *  Gang Lin, Thomas P.-Y. Yu, An Improved Vertex Caching Scheme for 3D Mesh Rendering,
   *  IEEE TVCG 12(4), 2006 and "Linear-Speed Vertex Cache Optimisation",
   *  Tom Forsyth, RAD Game Tools (http://home.comcast.net/~tom_forsyth/papers/fast_vert_cache_opt.html)
   */
  class triangle_mesh_greedy_stripifier: public triangle_mesh_stripifier {
  public: // Creation

    triangle_mesh_greedy_stripifier();
    ~triangle_mesh_greedy_stripifier();
    
  public: // Input

    virtual void begin_input();
    virtual void insert_input_triangle(std::size_t i0,
				       std::size_t i1,
				       std::size_t i2);
    virtual void end_input();
    
    std::size_t effective_cache_size() const;

  protected: // Stripification

    typedef stl_map_selector<std::size_t, std::size_t>::unordered_t vertex_id_to_vertex_index_map_t;

    vertex_id_to_vertex_index_map_t vertex_id_to_vertex_index_;
    std::vector<detail::tmg_vertex>   vertices_;
    std::vector<detail::tmg_triangle> triangles_;
    std::list<std::size_t> lru_cache_;

    std::size_t active_triangle_count_;
    float best_score_;
    std::size_t best_tri_;

    void stripify();
    void stripify_cleanup();
    void stripify_init();
    float stripify_vertex_score(std::size_t vidx) const;
    void stripify_commit_triangle(std::size_t tidx);
    void stripify_output_strip_grow(std::size_t vid0, std::size_t vid1, std::size_t vid2);
    void stripify_update_scores(std::size_t tidx);
  };
    

  // =============================================================================
  // triangle_mesh_backtracking_stripifier
  // =============================================================================
  
  /**
   *  A high quality, but slow, triangle mesh stripifier.
   */
  class triangle_mesh_backtracking_stripifier: public triangle_mesh_stripifier {
  public: // Creation

    triangle_mesh_backtracking_stripifier();
    ~triangle_mesh_backtracking_stripifier();
    
  public: // Input

    virtual void begin_input();
    virtual void insert_input_triangle(std::size_t i0,
				       std::size_t i1,
				       std::size_t i2);
    virtual void end_input();
    
  protected: // Stripification: types

    struct edata_t {
    };
    struct vdata_t {
    };

    static const std::size_t Tdata_mark_not_visited     = 0;
    static const std::size_t Tdata_mark_committed       = 1;
    static const std::size_t Tdata_mark_experiment_none = 2;
    static const std::size_t Tdata_mark_experiment_base = 3;
    
    struct tdata_t {
      std::size_t mark_;
      inline tdata_t(): mark_(Tdata_mark_not_visited) {}
      inline std::size_t mark() const { return mark_; }
      inline void set_mark(std::size_t x) { mark_ = x; }
    };
    
    typedef triangle_mesh<vdata_t, edata_t, tdata_t> tmesh_t;
    typedef tmesh_t::vertex_t            vertex_t;
    typedef tmesh_t::half_edge_t         half_edge_t;
    typedef tmesh_t::edge_t              edge_t;
    typedef tmesh_t::triangle_t          triangle_t;
    typedef std::set<triangle_t> triangle_set_t;

    typedef fifo_cache_simulator<vertex_t> vertex_cache_t;
    
  protected: // Data
    tmesh_t                tmesh_;
    
    triangle_set_t triangles_by_degree_[5]; // 4: >3 neighboring triangles
    vertex_cache_t vertex_cache_;

    std::size_t stripify_current_experiment_;
    
  protected: // Algorithm
    
    void stripify();
    
    bool stripify_is_not_visited(const triangle_t& tri,
                                 std::size_t mark = Tdata_mark_experiment_none) const;
    
    std::size_t stripify_degree(const triangle_t& tri,
                                std::size_t mark = Tdata_mark_experiment_none) const;
      
    void stripify_triangle_buckets_init();
    void stripify_triangle_buckets_insert(const triangle_t& tri);
    void stripify_triangle_buckets_erase(const triangle_t& tri);
    void stripify_cleanup();
    void stripify_candidate_roots_in(std::vector<triangle_t>& roots) const;
    void stripify_grow();
    void stripify_commit_unconnected_triangles();
    void stripify_commit(const triangle_t& tri);
    void stripify_commit(const strip_t& best_strip);
    void stripify_grow_in(strip_t& best_strip,
                          const std::vector<triangle_t>& roots);
    std::size_t stripify_strip_grow_cost(std::size_t i0,
                                         std::size_t i1,
                                         std::size_t i2,
                                         std::size_t experiment_mark,
                                         const strip_t& strip,
                                         const vertex_cache_t& vertex_cache) const;
    void stripify_strip_grow(std::size_t i0,
                             std::size_t i1,
                             std::size_t i2,
                             std::size_t experiment_mark,
                             strip_t& strip,
                             vertex_cache_t& vertex_cache);
    void stripify_best_next_tri_in(triangle_t&  best_tri,
                                   std::size_t& best_tri_cost,
                                   std::size_t experiment_mark,
                                   const strip_t& strip,
                                   const vertex_cache_t& vertex_cache) const;
    
    void stripify_grow_in(strip_t& best_strip,
                          double& cost,
                          const vertex_t& i0,
                          const vertex_t& i1,
                          const vertex_t& i2);
    
  };

}

#endif
