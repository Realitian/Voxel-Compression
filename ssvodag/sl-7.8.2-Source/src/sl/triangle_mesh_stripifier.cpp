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
#include <sl/triangle_mesh_stripifier.hpp>
#include <algorithm>
#include <cassert>

namespace sl {
  
  // =============================================================================
  // triangle_mesh_stripifier
  // =============================================================================

  // -----------------------------------------------------------------------
  // Creation/destruction
  // -----------------------------------------------------------------------

  triangle_mesh_stripifier::triangle_mesh_stripifier() {
    cache_size_ = 32;
    is_during_input_ = false;
    is_during_output_ = false;
    output_strip_cursor_ = 0;
    output_vertex_cursor_ = 0;
  }
  
  triangle_mesh_stripifier::~triangle_mesh_stripifier() {
  }

  // -----------------------------------------------------------------------
  // Parameterization
  // -----------------------------------------------------------------------

  void triangle_mesh_stripifier::set_cache_size(std::size_t x) {
    cache_size_ = x;
  }

  // -----------------------------------------------------------------------
  // Input
  // -----------------------------------------------------------------------
  
  bool triangle_mesh_stripifier::is_during_input() const {
    return is_during_input_;
  }


  // -----------------------------------------------------------------------
  // Output
  // -----------------------------------------------------------------------

  std::size_t triangle_mesh_stripifier::output_strip_count() const {
    return tstrips_.size();
  }
  
  void triangle_mesh_stripifier::begin_output() {
    assert(!is_during_input());
    assert(!is_during_output());
    is_during_output_ = true;
    
    output_strip_cursor_ = 0;
    output_vertex_cursor_ = 0;
    
    assert(is_during_output());
  }
    
  bool triangle_mesh_stripifier::is_during_output() const {
    return is_during_output_;
  }
  
  bool triangle_mesh_stripifier::has_output_strip() const {
    assert(is_during_output());
    return output_strip_cursor_ < tstrips_.size();
  }

  bool triangle_mesh_stripifier::has_output_vertex() const {
    assert(is_during_output());
    return
      output_strip_cursor_ < tstrips_.size() &&
      output_vertex_cursor_ < tstrips_[output_strip_cursor_].size();
  }

  void  triangle_mesh_stripifier::begin_output_strip() {
    output_vertex_cursor_ = 0;
  }

  void  triangle_mesh_stripifier::end_output_strip() {
    ++output_strip_cursor_;
  }
  
  std::size_t triangle_mesh_stripifier::get_output_strip_vertex() {
    assert(is_during_output());
    assert(has_output_vertex());
    std::size_t result = tstrips_[output_strip_cursor_][output_vertex_cursor_];
    ++output_vertex_cursor_;
    return result;
  }
  
  void triangle_mesh_stripifier::end_output() {
    assert(is_during_output());
    tstrips_.clear();
  }
  
  std::size_t triangle_mesh_stripifier::current_output_strip_index() const {
    return output_strip_cursor_;
  }
  
  std::size_t triangle_mesh_stripifier::current_output_strip_vertex_index() const {
    return output_vertex_cursor_;
  }

  // -----------------------------------------------------------------------
  // Strip growing
  // -----------------------------------------------------------------------
  
  void triangle_mesh_stripifier::strip_grow_in(strip_t& strip,
					       std::size_t i0,
					       std::size_t i1,
					       std::size_t i2) {
    triangle_connectivity tri(i0,i1,i2);
    
    std::size_t old_strip_size = strip.size();
    if (old_strip_size == 0) {
      // Strip start
      strip.push_back(i0);
      strip.push_back(i1);
      strip.push_back(i2);
    } else {
      // Strip grow
      bool is_odd = (old_strip_size%2) == 1;
      std::size_t v0 = strip[old_strip_size-3];
      std::size_t v1 = strip[old_strip_size-2];
      std::size_t v2 = strip[old_strip_size-1];
      
      std::size_t v2_idx = tri.index_of(v2);
      if (v2_idx == tri.index_end()) {
	std::size_t v1_idx = tri.index_of(v1);
	if (v1_idx == tri.index_end()) {
	  // Strip restart
	  strip.push_back(v2);
	  strip.push_back(tri[0]);
	  if (is_odd) {
	    strip.push_back(tri[0]);
	    strip.push_back(tri[2]);
	    strip.push_back(tri[1]);
	  } else {
	    strip.push_back(tri[0]);
	    strip.push_back(tri[1]);
	    strip.push_back(tri[2]);
	  }
	} else {
	  // Rename vertices in canonical order
	  std::size_t va = tri[(v1_idx+0)%3]; // = v1
	  std::size_t vb = tri[(v1_idx+1)%3];
	  std::size_t vc = tri[(v1_idx+2)%3];
	  strip.push_back(va);
	  if (is_odd) {
	    strip.push_back(va);
	    strip.push_back(vb);
	    strip.push_back(vc);
	  } else {
	    strip.push_back(va);
	    strip.push_back(vc);
	    strip.push_back(vb);
	  }
	}
      } else {
        // Rename vertices in canonical order
        std::size_t va = tri[(v2_idx+0)%3]; // = v2
        std::size_t vb = tri[(v2_idx+1)%3];
        std::size_t vc = tri[(v2_idx+2)%3];
#if 0
        std::cerr << "012 = " << v0 << " " << v1 << " " << v2 << std::endl;
        std::cerr << "ABC = " << va << " " << vb << " " << vc << std::endl;
#endif
        if (v1 == vb) {
          if (is_odd) {
            // ?BA->BAC
            //std::cerr << "BAC" << std::endl;
            strip.push_back(vc);
          } else {
            // ?BA->BABC
            //std::cerr << "BAABC" << std::endl;
            strip.push_back(va);
            strip.push_back(vb);
            strip.push_back(vc);
          }
        } else if (v1 == vc) {
          if (is_odd) {
            // ?CA->CACB
            //std::cerr << "CAACB" << std::endl;
            strip.push_back(va);
            strip.push_back(vc);
            strip.push_back(vb);
          } else {
            // ?CA->CAB
            strip.push_back(vb);
          }
        } else if (v0 == vc) {
          if (is_odd) {
            // C?A -> C?AACB
            // FIXME: Swap earlier tri
            //std::cerr << "C?AACB" << std::endl;
            strip.push_back(va);
            strip.push_back(vc);
            strip.push_back(vb);          
          } else {
            // C?B -> C?AABC
            // FIXME: Swap earlier tri
            //std::cerr << "C?AABC" << std::endl;
            strip.push_back(va);
            strip.push_back(vb);
            strip.push_back(vc);          
          }
        } else {
          if (is_odd) {
            // Strip restart
            //std::cerr << "??AACB" << std::endl;
            strip.push_back(va);
            strip.push_back(vc);
            strip.push_back(vb);
          } else {
            //std::cerr << "??AABC" << std::endl;
            strip.push_back(va);
            strip.push_back(vb);
            strip.push_back(vc);
          }
        }
      }
    }
  }

  // =============================================================================
  // triangle_mesh_greedy_stripifier
  // =============================================================================

  // -----------------------------------------------------------------------
  // Creation/destruction
  // -----------------------------------------------------------------------

  triangle_mesh_greedy_stripifier::triangle_mesh_greedy_stripifier() {
    best_score_ = -9999.0f;
    best_tri_ = std::size_t(-1);
    active_triangle_count_ = 0;
    cache_size_ = 32;
  }
  
  triangle_mesh_greedy_stripifier::~triangle_mesh_greedy_stripifier() {
  }

  // -----------------------------------------------------------------------
  // Parameterization
  // -----------------------------------------------------------------------

  std::size_t triangle_mesh_greedy_stripifier::effective_cache_size() const {
    return sl::median(cache_size_, std::size_t(3), std::size_t(32));
  }

  // -----------------------------------------------------------------------
  // Input
  // -----------------------------------------------------------------------

  void triangle_mesh_greedy_stripifier::begin_input() {
    assert(!is_during_input());
    assert(!is_during_output());

    tstrips_.clear();
    stripify_cleanup();

    is_during_input_ = true;
    assert(is_during_input());
  }
  
  void triangle_mesh_greedy_stripifier::insert_input_triangle(std::size_t i0,
							      std::size_t i1,
							      std::size_t i2) {
    SL_TRACE_OUT(1) << i0 << " " << i1 << " " << i2 << std::endl;

    assert(is_during_input());
    if (i0 != i1 && i0 != i2 && i1 != i2) {
      // Non degenerate triangle, insert

      // Insert tri
      std::size_t tri_idx = triangles_.size();
      triangles_.push_back(detail::tmg_triangle());
      detail::tmg_triangle& tmg_tri = triangles_.back();

      // Insert vertices
      std::size_t vtx_idx[3];
      std::size_t vtx_id[3]; vtx_id[0] = i0; vtx_id[1] = i1; vtx_id[2] = i2;
      for (std::size_t k=0; k<3; ++k) {
	vertex_id_to_vertex_index_map_t::iterator vtx_it = vertex_id_to_vertex_index_.find(vtx_id[k]);
	if (vtx_it==vertex_id_to_vertex_index_.end()) {
	  // New vertex
	  vtx_idx[k] = vertices_.size();
	  vertices_.push_back(detail::tmg_vertex(vtx_id[k]));
	  vertex_id_to_vertex_index_[vtx_id[k]] = vtx_idx[k];
	} else {
	  vtx_idx[k] = vtx_it->second;
	}
	vertices_[vtx_idx[k]].active_incident_triangles_.insert(tri_idx);
      }

      tmg_tri.connectivity_ = triangle_connectivity(vtx_idx[0], vtx_idx[1], vtx_idx[2]);
    }
  }
    
  void triangle_mesh_greedy_stripifier::end_input() {
    assert(is_during_input());
    is_during_input_ = false;
    stripify();
    assert(!is_during_input());
  }

  // -----------------------------------------------------------------------
  // Stripification
  // -----------------------------------------------------------------------

  void triangle_mesh_greedy_stripifier::stripify() {
    stripify_init();
    while (best_tri_!=std::size_t(-1)) {
      assert(best_tri_<triangles_.size());
      stripify_commit_triangle(best_tri_);
    }
    stripify_cleanup();
    assert(active_triangle_count_==0);
  }
  
  void triangle_mesh_greedy_stripifier::stripify_cleanup() {
    triangles_.clear();
    vertices_.clear();
    vertex_id_to_vertex_index_.clear();
    lru_cache_.clear();
    best_score_ = -9999.0f;
    best_tri_ = std::size_t(-1);
    active_triangle_count_ = 0;
  }
  
  void triangle_mesh_greedy_stripifier::stripify_init() {
    SL_TRACE_OUT(1) << "Vertices: " << vertices_.size() << " Triangles: " << triangles_.size() << std::endl;
    
    // Init vertex scores
    for (std::size_t i=0; i<vertices_.size(); ++i) {
      vertices_[i].score_ = stripify_vertex_score(i);
    }
    best_score_ = -9999.0f;
    best_tri_ = std::size_t(-1);
    active_triangle_count_ = 0;
    for (std::size_t i=0; i<triangles_.size(); ++i) {
      float s_i =
	vertices_[triangles_[i].connectivity_[0]].score_ +
	vertices_[triangles_[i].connectivity_[1]].score_ +
	vertices_[triangles_[i].connectivity_[2]].score_;
      ++active_triangle_count_;
      triangles_[i].is_active_ = true;
      triangles_[i].score_ = s_i;
      if (s_i > best_score_) {
	best_tri_ = i;
	best_score_ = s_i;
      }
    }
    SL_TRACE_OUT(1) << "Best: " << best_tri_ << " with score " << best_score_ << std::endl;
  }
  
  float triangle_mesh_greedy_stripifier::stripify_vertex_score(std::size_t vidx) const {
    const float Cache_decay_power   = 1.5f;
    //const float Last_tri_score      = 0.75f; // Unused
    const float Valence_boost_scale = 2.0f;
    const float Valence_boost_power = -0.5f;

    const detail::tmg_vertex& vtx = vertices_[vidx];
    const std::size_t incident_tricount = vtx.active_valence();
      
    float result = -1.0f;
    if (incident_tricount > 0) {
 
      result = 0.0f;
      if (vtx.lru_position_ < 0) {
	// Not in cache, no score
      } else if (effective_cache_size()<=3) {
	// We are simulating using a "cacheless" simulation - give a high cache score only to
	// the last edge
	if (vtx.lru_position_<2) {
	  result = 10.0f * Valence_boost_scale;
	}
      } else {
	assert(vtx.lru_position_<int(effective_cache_size()));
	assert(effective_cache_size()>0);
	result = std::pow(1.0f - float(vtx.lru_position_)/float(effective_cache_size()), Cache_decay_power);
      }

      // Bonus points for having a low number of tris still to
      // use the vertex, so we get rid of lone verts quickly.

      result += Valence_boost_scale * std::pow(float(incident_tricount), Valence_boost_power);
    }

    SL_TRACE_OUT(1) << vtx.id_ << "(tricount:" << incident_tricount << " lru: " << vtx.lru_position_ << ") -> " << result << std::endl;
    return result;
  }

  void triangle_mesh_greedy_stripifier::stripify_commit_triangle(std::size_t tidx) {
    assert(tidx<triangles_.size());
    
    detail::tmg_triangle& tmg_tri = triangles_[tidx];

    // Reorder so that vid0 is of valence 0 or has highest score
    const triangle_connectivity& vidx = tmg_tri.connectivity_;
    std::size_t  vid[3];
    
    std::size_t best = 0;
    if (vertices_[vidx[1]].active_valence()==1) best=1;
    if (vertices_[vidx[2]].active_valence()==1) best=2;
    if (vertices_[vidx[best]].active_valence()>1) {
      best = 0;
      if (vertices_[vidx[1]].score_>vertices_[vidx[best]].score_) best=1;
      if (vertices_[vidx[2]].score_>vertices_[vidx[best]].score_) best=2;
    }
    vid[0] = vertices_[vidx[(best+0)%3]].id_; 
    vid[1] = vertices_[vidx[(best+1)%3]].id_; 
    vid[2] = vertices_[vidx[(best+2)%3]].id_; 

    stripify_output_strip_grow(vid[0], vid[1], vid[2]);
    stripify_update_scores(tidx);
  }

  void triangle_mesh_greedy_stripifier::stripify_output_strip_grow(std::size_t vid0,
								   std::size_t vid1,
								   std::size_t vid2) {

    if (tstrips_.empty()) {
      tstrips_.push_back(strip_t());
    }
    strip_t& strip = tstrips_.back();

    strip_grow_in(strip, vid0, vid1, vid2);
    
    SL_TRACE_OUT(1) << vid0 << " " << vid1 << " " << vid2 << " => sz=" << strip.size() << std::endl;
  }
  
  void triangle_mesh_greedy_stripifier::stripify_update_scores(std::size_t tidx) {
    detail::tmg_triangle& tmg_tri = triangles_[tidx];

    // Reorder triangle in the same order added to tri strip

    triangle_connectivity vidx = tmg_tri.connectivity_;
#if 0
    triangle_connectivity  vid = triangle_connectivity(vertices_[vidx[0]].id_,
						       vertices_[vidx[1]].id_,
						       vertices_[vidx[2]].id_);
#endif
    
    strip_t& strip = tstrips_.back();
    std::size_t N = strip.size();
    std::size_t vidx_strip[3];

    vidx_strip[0] = vertex_id_to_vertex_index_[strip[N-3]];
    vidx_strip[1] = vertex_id_to_vertex_index_[strip[N-2]];
    vidx_strip[2] = vertex_id_to_vertex_index_[strip[N-1]];
    // FIXME Check that the indices correspond to the triangle...
    
    // Update lru cache
    for (std::size_t k=0; k<3; ++k) {
      std::list<std::size_t>::iterator it = std::find(lru_cache_.begin(),
						      lru_cache_.end(),
						      vidx_strip[k]);
      if (it != lru_cache_.end()) {
	lru_cache_.erase(it);
      }
      lru_cache_.push_front(vidx_strip[k]);
    }

    // Update valence
    tmg_tri.is_active_ = false;
    vertices_[vidx[0]].active_incident_triangles_.erase(tidx);
    vertices_[vidx[1]].active_incident_triangles_.erase(tidx);
    vertices_[vidx[2]].active_incident_triangles_.erase(tidx);
    --active_triangle_count_;
 
    // Update scores of vertices in cache and vertices soon to be removed
    int lru_pos = 0;
    for (std::list<std::size_t>::iterator vit = lru_cache_.begin();
	 vit != lru_cache_.end();
	 ++vit, ++lru_pos) {
      std::size_t vidx = *vit;
      vertices_[vidx].lru_position_ = ((lru_pos >= int(effective_cache_size())) ? -1 : lru_pos);
      vertices_[vidx].score_ = stripify_vertex_score(vidx);
    }

    // Update scores of triangles references by vertices in cache and vertices soon to be removed
    best_score_ = -9999.0f;
    best_tri_ = std::size_t(-1);
    for (std::list<std::size_t>::iterator vit = lru_cache_.begin();
	 vit != lru_cache_.end();
	 ++vit) {
      std::size_t vidx = *vit;
      for (sl::vector_set<std::size_t>::iterator tit = vertices_[vidx].active_incident_triangles_.begin();
	   tit != vertices_[vidx].active_incident_triangles_.end();
	   ++tit) {
	std::size_t tidx = *tit;
	if (triangles_[tidx].is_active_) {
	  float s_i =
	    vertices_[triangles_[tidx].connectivity_[0]].score_ +
	    vertices_[triangles_[tidx].connectivity_[1]].score_ +
	    vertices_[triangles_[tidx].connectivity_[2]].score_;
	  triangles_[tidx].score_ = s_i;

	  if (s_i > best_score_) {
	    best_tri_ = tidx;
	    best_score_ = s_i;
	  } // if s_i > best score
	} else {
	  best_score_ = -9999.0f;
	} // if triangle is active
      } // for each incident triangle
    } // for each triangle

    // Shrink back cache to maximum size
    while (lru_cache_.size() > effective_cache_size()) {
      lru_cache_.pop_back();
    }
    
    // If unusually low score, also check triangles not in cache
    if ((active_triangle_count_>0) && ((best_tri_!=size_t(-1)) || (best_score_ < 0.1f))) {
      SL_TRACE_OUT(1) << "Unusually low score: " << best_score_ << std::endl;
      best_score_ = -9999.0f;
      best_tri_ = std::size_t(-1);
      for (std::size_t i=0; i<triangles_.size(); ++i) {
	if (triangles_[i].is_active_) {
	  float s_i = triangles_[i].score_;
	  if (s_i > best_score_) {
	    best_tri_ = i;
	    best_score_ = s_i;
	  }
	}
      }
      SL_TRACE_OUT(1) << "New low score: " << best_score_ << std::endl;
    } 
  }
  
  // =============================================================================
  // triangle_mesh_backtracking_stripifier
  // =============================================================================

  // -----------------------------------------------------------------------
  // Creation/destruction
  // -----------------------------------------------------------------------

  triangle_mesh_backtracking_stripifier::triangle_mesh_backtracking_stripifier() {
  }
  
  triangle_mesh_backtracking_stripifier::~triangle_mesh_backtracking_stripifier() {
  }

  // -----------------------------------------------------------------------
  // Input
  // -----------------------------------------------------------------------

  void triangle_mesh_backtracking_stripifier::begin_input() {
    assert(!is_during_input());
    assert(!is_during_output());

    tmesh_.clear();
    tstrips_.clear();
    
    is_during_input_ = true;
    assert(is_during_input());
  }
  
  void triangle_mesh_backtracking_stripifier::insert_input_triangle(std::size_t i0,
								    std::size_t i1,
								    std::size_t i2) {
    assert(is_during_input());
    if (i0 != i1 && i0 != i2 && i1 != i2) {
      tmesh_.insert_triangle(tmesh_t::triangle_t(i0,i1,i2));
    }
  }
    
  void triangle_mesh_backtracking_stripifier::end_input() {
    assert(is_during_input());
    is_during_input_ = false;
    stripify();
    assert(!is_during_input());
  }

  // -----------------------------------------------------------------------
  // Stripification
  // -----------------------------------------------------------------------

  void triangle_mesh_backtracking_stripifier::stripify() {
    stripify_current_experiment_ = Tdata_mark_experiment_none;
    stripify_triangle_buckets_init();
    stripify_grow();
    stripify_cleanup();
  }

  bool triangle_mesh_backtracking_stripifier::stripify_is_not_visited(const triangle_t& tri,
                                                  std::size_t mark) const {
    const tmesh_t::triangle_data_t* tridata_ptr = tmesh_.triangle_data(tri);
    assert(tridata_ptr);
    return
      tridata_ptr->mark() != Tdata_mark_committed &&
      tridata_ptr->mark() != mark;
  }

  std::size_t triangle_mesh_backtracking_stripifier::stripify_degree(const triangle_t& tri,
                                                        std::size_t mark) const {
    std::size_t result = 0;
    for (std::size_t i=0; i<3; ++i) {
      const tmesh_t::small_triangle_set_t* e_tri_ptr = tmesh_.edge_triangles(tri.undirected_edge(i));
      if (e_tri_ptr) {
        for (tmesh_t::small_triangle_set_t::const_iterator t_it = e_tri_ptr->begin();
             t_it != e_tri_ptr->end();
             ++t_it) {
          if (*t_it != tri && stripify_is_not_visited(*t_it, mark)) {
            ++result;
          }
        }
      }
    }
    return result;
  }

  void triangle_mesh_backtracking_stripifier::stripify_triangle_buckets_init() {
    for (tmesh_t::triangle_map_t::const_iterator it = tmesh_.triangle_map().begin();
         it != tmesh_.triangle_map().end();
         ++it) {
      stripify_triangle_buckets_insert(it->first);
    }
  }

  void triangle_mesh_backtracking_stripifier::stripify_triangle_buckets_insert(const triangle_t& tri) {
    std::size_t triangle_degree = stripify_degree(tri);
    triangles_by_degree_[std::min(std::size_t(4), triangle_degree)].insert(tri);
  }

  void triangle_mesh_backtracking_stripifier::stripify_triangle_buckets_erase(const triangle_t& tri) {
    for (std::size_t i=0; i<5; ++i) {
      triangles_by_degree_[i].erase(tri);
    }
  }
  
  void triangle_mesh_backtracking_stripifier::stripify_cleanup() {
    tmesh_.clear();
    for (std::size_t i=0; i<5; ++i) {
      triangles_by_degree_[i].clear();
    }
  }

  void triangle_mesh_backtracking_stripifier::stripify_candidate_roots_in(std::vector<triangle_t>& roots) const {
    roots.clear();
    // Roots from cache
    {
      bool found = false;
      std::size_t i =0;
      for (vertex_cache_t::const_iterator vit = vertex_cache_.begin();
           !found && vit != vertex_cache_.end();
           ++vit) {
        const tmesh_t::small_triangle_set_t* v_tri_ptr = tmesh_.vertex_triangles(*vit);
        if (v_tri_ptr) {
          for (tmesh_t::small_triangle_set_t::const_iterator t_it = v_tri_ptr->begin();
               !found && t_it != v_tri_ptr->end();
               ++t_it) {
            const triangle_t& tri = *t_it;
            if ((std::find(roots.rbegin(), roots.rend(), tri) == roots.rend()) &&
                stripify_is_not_visited(tri)) {
              if (vertex_cache_.hit_count(tri[0])+
                  vertex_cache_.hit_count(tri[1])+
                  vertex_cache_.hit_count(tri[2]) >= 2) {
                roots.push_back(tri);
                //found = true;
              }
            }
          }
        }
        ++i;
        //found = !roots.empty() && 3*i>vertex_cache_.capacity();
      }
    }

    // Roots from mesh
    {
      bool found = false;
      for (std::size_t i=0; i<5 && !found; ++i) {
        for (triangle_set_t::const_iterator t_it = triangles_by_degree_[i].begin();
             !found && t_it != triangles_by_degree_[i].end();
             ++t_it) {
          const triangle_t& tri = *t_it;
          if ((std::find(roots.rbegin(), roots.rend(), tri) == roots.rend()) && stripify_is_not_visited(tri)) {
            //std::cerr << "RD_" << i << ": " << tri[0]<< " " << tri[1] << " " << tri[2] << std::endl;
            roots.push_back(tri);
            found = true;
          }
        }
      }
    }
  }
  
  void triangle_mesh_backtracking_stripifier::stripify_grow() {
    std::vector<triangle_t> roots;

    vertex_cache_.clear();
    vertex_cache_.resize(cache_size_);

    stripify_commit_unconnected_triangles();
    stripify_candidate_roots_in(roots);
    while (!roots.empty()) {
      strip_t strip_i;
      // Grow from selected roots
      stripify_grow_in(strip_i, roots);

      // If stuck, grow from mesh triangles until we find a strip 
      for (std::size_t i=0; i<5 && strip_i.empty(); ++i) {
        for (triangle_set_t::iterator t_it = triangles_by_degree_[i].begin();
             strip_i.empty() && t_it != triangles_by_degree_[i].end();
             ++t_it) {
          const triangle_t& tri = *t_it;
          roots.clear();
          roots.push_back(tri);
          stripify_grow_in(strip_i, roots);
        }
      }
      
      // If stuck, commit generate a single tri root to advance
      for (std::size_t i=0; i<5 && strip_i.empty(); ++i) {
        for (triangle_set_t::iterator t_it = triangles_by_degree_[i].begin();
             strip_i.empty() && t_it != triangles_by_degree_[i].end();
             ++t_it) {
          const triangle_t& tri = *t_it;
          strip_i.push_back(tri[0]);
          strip_i.push_back(tri[1]);
          strip_i.push_back(tri[2]);
          //std::cerr << "STUCK!" << std::endl;
        }
      }
      assert(!strip_i.empty());
      
      stripify_commit(strip_i);
      stripify_commit_unconnected_triangles();
      stripify_candidate_roots_in(roots);
    }
  }

  void triangle_mesh_backtracking_stripifier::stripify_commit(const triangle_t& tri) {
    tmesh_t::triangle_data_t* tridata_ptr = tmesh_.triangle_data(tri);
    std::size_t tri_mark = tridata_ptr->mark();
    tridata_ptr->set_mark(Tdata_mark_committed);
    stripify_triangle_buckets_erase(tri);
    for (std::size_t i=0; i<3; ++i) {
      const tmesh_t::small_triangle_set_t* e_tri_ptr = tmesh_.edge_triangles(tri.undirected_edge(i));
      if (e_tri_ptr) {
        for (tmesh_t::small_triangle_set_t::const_iterator t_it = e_tri_ptr->begin();
             t_it != e_tri_ptr->end();
             ++t_it) {
          const triangle_t& tri2 = *t_it;
          if (tri2 != tri && stripify_is_not_visited(tri2, tri_mark)) {
            stripify_triangle_buckets_erase(tri2);
            stripify_triangle_buckets_insert(tri2);
          }
        }
      }
    }
  }

  void triangle_mesh_backtracking_stripifier::stripify_commit(const strip_t& best_strip) {
    // Remove strip from mesh/buckets and update vertex cache
    //std::cerr << "COMMIT: " << best_strip.size() << std::endl;
    for (std::size_t i=2; i<best_strip.size(); ++i) {
      vertex_t v0= best_strip[i-2];
      vertex_t v1= best_strip[i-1];
      vertex_t v2= best_strip[i-0];
      triangle_t t = (i%2==0) ? triangle_t(v0,v1,v2) : triangle_t(v0,v2,v1);
      if (t.is_valid()) {
        stripify_commit(t);
        if (i>2) {
          vertex_cache_.insert(v0);
          vertex_cache_.insert(v1);
        }
        vertex_cache_.insert(v2);
      }
    }
    // Record strip
    tstrips_.push_back(best_strip);
  }

  void triangle_mesh_backtracking_stripifier::stripify_commit_unconnected_triangles() {
    while (!triangles_by_degree_[0].empty()) {
      const triangle_t& tri = *(triangles_by_degree_[0].begin());
      strip_t strip;
      strip.push_back(tri[0]);
      strip.push_back(tri[1]);
      strip.push_back(tri[2]);
      stripify_commit(strip);
    }
  }

  void triangle_mesh_backtracking_stripifier::stripify_grow_in(strip_t& best_strip,
							       const std::vector<triangle_t>& roots) {
    double best_strip_cost = 1e30;
    best_strip.clear();

    std::size_t old_experiments_mark = stripify_current_experiment_;
    for (std::vector<triangle_t>::const_iterator root_it = roots.begin();
         root_it != roots.end();
         ++root_it) {
      const triangle_t& root_tri = *root_it;
      if (tmesh_.triangle_data(root_tri)->mark() <= old_experiments_mark) {
        for (std::size_t i=0; i<3; ++i) {
          strip_t strip_i;
          double strip_i_cost;
          stripify_grow_in(strip_i,
                           strip_i_cost,
                           root_tri[(i+0)%3],
                           root_tri[(i+1)%3],
                           root_tri[(i+2)%3]);
          if (best_strip.empty() || (strip_i_cost < best_strip_cost)) {
            if (strip_i.size() > 3) {
              std::swap(best_strip, strip_i);
              std::swap(best_strip_cost, strip_i_cost);
            }
          }
        }
      }
    }
  }
  
  std::size_t triangle_mesh_backtracking_stripifier::stripify_strip_grow_cost(std::size_t i0,
									      std::size_t i1,
									      std::size_t i2,
									      std::size_t experiment_mark,
									      const strip_t& strip,
									      const vertex_cache_t& vertex_cache) const {
    std::size_t result = 0;

    // Count neighbors
    triangle_t candidate_tri(i0,i1,i2);
    std::size_t triangle_neighbors = stripify_degree(candidate_tri, experiment_mark);
    result += triangle_neighbors;
    
    // If we are incorporating a lonely tri, we are done, othewise we check
    // other criteria
    if (result > 0) {
      // Cost due to the increase in the number of indices
      std::size_t index_cost = 0;
      if (strip.empty()) {
        index_cost = 10+3;
      } else {
        bool is_odd = (strip.size()%2) == 1;
      
        vertex_t v0 = strip[strip.size()-3];
        vertex_t v1 = strip[strip.size()-2];
        vertex_t v2 = strip[strip.size()-1];

        tmesh_t::vertex_t v2_idx = candidate_tri.index_of(v2);
        if (v2_idx == candidate_tri.index_end()) {
          // Strip restart
          index_cost+= is_odd ? 6 : 5;
        } else {
          // Rename vertices in canonical order
          //vertex_t va = candidate_tri[v2_idx];
          vertex_t vb = candidate_tri[(v2_idx+1)%3];
          vertex_t vc = candidate_tri[(v2_idx+2)%3];
          
          if (v1 == vb) {
            // Natural direction
            index_cost+= is_odd ? 1 : 3;
          } else if (v1 == vc) {
            // Reverse direction
            index_cost+= is_odd ? 3 : 1;
          } else if (v0 == vc) {
            // Edge Swap
            index_cost+= is_odd ? 3 : 3;
          } else {
            // Strip restart
            index_cost+= is_odd ? 3 : 3;
          }
        }
      }
      result += index_cost;
#if 1
      // Cost due to cache misses
      std::size_t miss_count = 0;
      if (!vertex_cache.has(candidate_tri[0])) ++miss_count;
      if (!vertex_cache.has(candidate_tri[1])) ++miss_count;
      if (!vertex_cache.has(candidate_tri[2])) ++miss_count;
      result += 6*miss_count; // FIXME: parameterize weight
#endif
    }
    return result;
  }

  void triangle_mesh_backtracking_stripifier::stripify_strip_grow(std::size_t i0,
								  std::size_t i1,
								  std::size_t i2,
								  std::size_t experiment_mark,
								  strip_t& strip,
								  vertex_cache_t& vertex_cache) {
    std::size_t old_strip_size = strip.size();
    
    strip_grow_in(strip, i0, i1, i2);
    
    tmesh_.triangle_data(triangle_t(i0,i1,i2))->set_mark(experiment_mark);

    // Update cache
    for (std::size_t i=old_strip_size; i< strip.size(); ++i) {
      vertex_cache.insert(strip[i]);
    }
  }

  void triangle_mesh_backtracking_stripifier::stripify_best_next_tri_in(triangle_t&  best_tri,
									std::size_t& best_tri_cost,
									std::size_t experiment_mark,
									const strip_t& strip,
									const vertex_cache_t& vertex_cache) const {
    edge_t e12(strip[strip.size()-2], strip[strip.size()-1]);
    edge_t e20(strip[strip.size()-1], strip[strip.size()-3]);

    best_tri = triangle_t(0,0,0);
    best_tri_cost = std::size_t(-1);
    
    // Choose forward direction
    const tmesh_t::small_triangle_set_t* e_tri_ptr = tmesh_.edge_triangles(e12);
    if (e_tri_ptr) {
      for (tmesh_t::small_triangle_set_t::const_iterator t_it = e_tri_ptr->begin();
           t_it != e_tri_ptr->end();
           ++t_it) {
        const triangle_t& tri = *t_it;
        if (stripify_is_not_visited(tri, experiment_mark)) {
          std::size_t tri_cost = stripify_strip_grow_cost(tri[0], tri[1], tri[2],
                                                          experiment_mark,
                                                          strip,
                                                          vertex_cache);
          if (tri_cost<best_tri_cost) {
            best_tri = tri;
            best_tri_cost = tri_cost;
          }            
        }
      }
    }

    // If not root tri, evaluate swap direction
    //    if (best_tri_cost >0 && strip.size()>3) {
    if (!best_tri.is_valid() && strip.size()>3) {
      const tmesh_t::small_triangle_set_t* e_tri_ptr = tmesh_.edge_triangles(e20);
      if (e_tri_ptr) {
        for (tmesh_t::small_triangle_set_t::const_iterator t_it = e_tri_ptr->begin();
             t_it != e_tri_ptr->end();
             ++t_it) {
          const triangle_t& tri = *t_it;
          if (stripify_is_not_visited(tri, experiment_mark)) {
            std::size_t tri_cost = stripify_strip_grow_cost(tri[0], tri[1], tri[2],
                                                            experiment_mark,
                                                            strip,
                                                            vertex_cache);
            //if (tri_cost<best_tri_cost) {
            if (tri_cost == 0) {
              best_tri = tri;
              best_tri_cost = tri_cost;
            }            
          }
        }
      }
    }
  }
  
  void triangle_mesh_backtracking_stripifier::stripify_grow_in(strip_t& strip,
							       double& cost,
							       const vertex_t& i0,
							       const vertex_t& i1,
							       const vertex_t& i2) {
    
    ++stripify_current_experiment_;
    std::size_t experiment_mark = stripify_current_experiment_;
    
    // Init cache
    vertex_cache_t vertex_cache = vertex_cache_;
    vertex_cache.reset_counts();
    std::size_t Max_miss_count =
      (vertex_cache.capacity()<3)
      ?
      std::size_t(-1)
      :
      std::max(vertex_cache.capacity(),std::size_t(8))-std::size_t(6);

    // Init strip
    std::size_t strip_tricount = 0;
    strip.clear();
    std::size_t incremental_cost = stripify_strip_grow_cost(i0,i1,i2,
                                                            experiment_mark,
                                                            strip,
                                                            vertex_cache);
    stripify_strip_grow(i0,i1,i2, experiment_mark, strip, vertex_cache);
    ++strip_tricount;
    
    bool done = false;
    while (!done) {
      triangle_t  best_tri(0,0,0);
      std::size_t best_tri_cost = std::size_t(-1);
      stripify_best_next_tri_in(best_tri, best_tri_cost,
                                experiment_mark,
                                strip,
                                vertex_cache);                             
      if (best_tri.is_valid()) {
        // Extend strip with best tri
        ++strip_tricount;
        stripify_strip_grow(best_tri[0], best_tri[1], best_tri[2], experiment_mark, strip, vertex_cache);
        incremental_cost += best_tri_cost;
        
        if (vertex_cache.miss_count() >= Max_miss_count) {
          // We laid down another full cache, stop growing unless we can gobble a lonely tri
          stripify_best_next_tri_in(best_tri, best_tri_cost,
                                    experiment_mark,
                                    strip,
                                    vertex_cache);
          if (best_tri_cost == 0) {
            ++strip_tricount;
            stripify_strip_grow(best_tri[0], best_tri[1], best_tri[2], experiment_mark, strip, vertex_cache);
            //incremental_cost += best_tri_cost;
          }
          done = true;            
        }
      } else {
        // Cannot find a good tri
        //std::cerr << "No tri" << std::endl;
        done = true;
      }
    }

    // Normalize cost
    cost = double(incremental_cost)/strip_tricount;
    //std::cerr << "C[" << i0 << " " << i1 << " " << i2 << "]: " << incremental_cost << "/" << strip_tricount << " = " << cost << std::endl;
  }
  
}


