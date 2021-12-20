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
#ifndef SL_TRIANGLE_MESH_HPP
#define SL_TRIANGLE_MESH_HPP

#include <sl/connectivity.hpp>
#include <sl/sorted_vector_set.hpp>
#include <sl/vector_set.hpp>
#include <sl/operators.hpp>
#include <fstream>
#include <stack>
#include <map>
#include <set>
#include <vector>
#include <sl/utility.hpp> // hashing functions
#include <sl/random.hpp> 
#include <sl/stl_container_selector.hpp>

namespace sl {
  // -----------------------------------------------------------------------
  /// Vertex with attributes
  template <class T_Data>
  class triangle_mesh_vertex_attribute {
  public:
    typedef T_Data                            data_t;
    typedef directed_edge_connectivity        half_edge_t;
    typedef undirected_edge_connectivity      edge_t;
    typedef triangle_connectivity             triangle_t;
    typedef sorted_vector_set<edge_t>         small_edge_set_t;
    typedef sorted_vector_set<triangle_t>     small_triangle_set_t;
  protected:
    small_edge_set_t      edge_set_;
    small_triangle_set_t  triangle_set_;
    data_t                data_;
  public:
    triangle_mesh_vertex_attribute() {
      edge_set_.reserve(6);
      triangle_set_.reserve(6);
    }
    
    const small_edge_set_t& edge_set() const {
      return edge_set_;
    }
    small_edge_set_t& edge_set() {
      return edge_set_;
    }

    const  small_triangle_set_t& triangle_set() const {
      return triangle_set_;
    }
    
    small_triangle_set_t& triangle_set() {
      return triangle_set_;
    }

    const data_t* data() const {
      return &data_;
    }

    data_t* data() {
      return &data_;
    }
  };
}

namespace sl {
  // -----------------------------------------------------------------------
  /// Edge with attributes
  template <class T_Data>
  class triangle_mesh_edge_attribute {
  public:
    typedef T_Data                            data_t;
    typedef directed_edge_connectivity        half_edge_t;
    typedef undirected_edge_connectivity      edge_t;
    typedef triangle_connectivity             triangle_t;
    typedef sorted_vector_set<edge_t>         small_edge_set_t;
    typedef sorted_vector_set<triangle_t>     small_triangle_set_t;
  protected:
    small_triangle_set_t  triangle_set_;
    data_t                data_;
  public:
    triangle_mesh_edge_attribute() {
      triangle_set_.reserve(2);
    }

    const  small_triangle_set_t& triangle_set() const {
      return triangle_set_;
    }
    
    small_triangle_set_t& triangle_set() {
      return triangle_set_;
    }

    const data_t* data() const {
      return &data_;
    }

    data_t* data() {
      return &data_;
    }
  };

  // -----------------------------------------------------------------------
  /// Triangle with attributes
  template <class T_Data>
  class triangle_mesh_triangle_attribute {
  public:
    typedef T_Data                            data_t;
    typedef directed_edge_connectivity        half_edge_t;
    typedef undirected_edge_connectivity      edge_t;
    typedef triangle_connectivity             triangle_t;
  protected:
    data_t                         data_;
  public:
    triangle_mesh_triangle_attribute() {
    }

    const data_t* data() const {
      return &data_;
    }

    data_t* data() {
      return &data_;
    }
  };
}

namespace sl {
  // -----------------------------------------------------------------------
  /**
   *  A general triangle mesh with attributes on triangles, vertices, or edges.
   *  Memory-hungry representation, but O(logN) topology update operations.
   */
  template <class T_Vertex_Data, class T_Edge_Data, class T_Triangle_Data>
  class triangle_mesh {
  public:
    typedef T_Vertex_Data   vertex_data_t;
    typedef T_Edge_Data     edge_data_t;
    typedef T_Triangle_Data triangle_data_t;
    
    typedef triangle_mesh<vertex_data_t, edge_data_t, triangle_data_t> this_t;

    typedef directed_edge_connectivity::vertex_t	vertex_t;
    typedef directed_edge_connectivity			half_edge_t;
    typedef undirected_edge_connectivity		edge_t;
    typedef triangle_connectivity			triangle_t;

    typedef triangle_mesh_vertex_attribute<vertex_data_t>      vertex_attribute_t;
    typedef triangle_mesh_edge_attribute<edge_data_t>          edge_attribute_t;
    typedef triangle_mesh_triangle_attribute<triangle_data_t>  triangle_attribute_t;

    typedef typename stl_map_selector<vertex_t,vertex_attribute_t>::unordered_t vertex_map_t;
    typedef typename vertex_map_t::iterator                    vertex_iterator_t;
    typedef typename vertex_map_t::const_iterator              const_vertex_iterator_t;

    typedef typename stl_map_selector<edge_t, edge_attribute_t>::unordered_t    edge_map_t;
    typedef typename edge_map_t::iterator                      edge_iterator_t;
    typedef typename edge_map_t::const_iterator                const_edge_iterator_t;
    
    typedef typename stl_map_selector<triangle_t, triangle_attribute_t>::unordered_t triangle_map_t;
    typedef typename triangle_map_t::iterator                  triangle_iterator_t;
    typedef typename triangle_map_t::const_iterator            const_triangle_iterator_t;

    typedef sorted_vector_set<vertex_t>                        small_vertex_set_t;
    typedef sorted_vector_set<edge_t>                          small_edge_set_t;
    typedef sorted_vector_set<triangle_t>                      small_triangle_set_t;
   
  protected:
    mutable random::std_irng_t irng_;

    vertex_map_t    vertex_map_;
    edge_map_t      edge_map_;
    triangle_map_t  triangle_map_;

  protected:

    const vertex_attribute_t   null_vertex_attribute_;
    const edge_attribute_t     null_edge_attribute_;
    const triangle_attribute_t null_triangle_attribute_;
    
  public:  // construction and destruction

    triangle_mesh () {
    }
    
    virtual ~triangle_mesh () {
      clear();
    }

    /** Used for operations that create new meshes from the current one.  This
     * allows derived class construction within the base class operations.
     */
    virtual this_t* new_mesh() const {
      return new this_t;
    }
    
  public: // accessors for sizes

    std::size_t vertex_count () const {
      return vertex_map_.size();
    }
    
    std::size_t edge_count () const {
      return edge_map_.size();
    }
    
    std::size_t triangle_count () const {
      return triangle_map_.size();
    }
    
  public: // insert/remove triangles
    
    /// Insert triangle, also inserting new edges and vertices if not already present
    virtual void insert_triangle(const triangle_t& tri) {
      // insert triangle
      std::pair<triangle_iterator_t,bool> t_it = triangle_map_.insert(std::make_pair(tri,null_triangle_attribute_));
      if (t_it.second) {
        edge_t e01(tri[0],tri[1]), e12(tri[1],tri[2]), e20(tri[2],tri[0]);
      
        // insert vertices
        std::pair<vertex_iterator_t,bool> v_it0 = vertex_map_.insert(std::make_pair(tri[0],null_vertex_attribute_));
        v_it0.first->second.edge_set().insert(e01);
        v_it0.first->second.edge_set().insert(e20);
        v_it0.first->second.triangle_set().insert(tri);
        
        std::pair<vertex_iterator_t,bool> v_it1 = vertex_map_.insert(std::make_pair(tri[1],null_vertex_attribute_));
        v_it1.first->second.edge_set().insert(e01);
        v_it1.first->second.edge_set().insert(e12);
        v_it1.first->second.triangle_set().insert(tri);
        
        std::pair<vertex_iterator_t,bool> v_it2 = vertex_map_.insert(std::make_pair(tri[2],null_vertex_attribute_));
        v_it2.first->second.edge_set().insert(e12);
        v_it2.first->second.edge_set().insert(e20);
        v_it2.first->second.triangle_set().insert(tri);

        // insert edges
        std::pair<edge_iterator_t,bool> e_it0 = edge_map_.insert(std::make_pair(e01,null_edge_attribute_));
        e_it0.first->second.triangle_set().insert(tri);
        
        std::pair<edge_iterator_t,bool> e_it1 = edge_map_.insert(std::make_pair(e12,null_edge_attribute_));
        e_it1.first->second.triangle_set().insert(tri);
        
        std::pair<edge_iterator_t,bool> e_it2 = edge_map_.insert(std::make_pair(e20,null_edge_attribute_));
        e_it2.first->second.triangle_set().insert(tri);
      }
      
      SL_ENSURE("Inserted", triangle_map_.find(tri) != triangle_map_.end());
      SL_ENSURE("Inserted", vertex_map_.find(tri[0]) != vertex_map_.end());
      SL_ENSURE("Inserted", vertex_map_.find(tri[1]) != vertex_map_.end());
      SL_ENSURE("Inserted", vertex_map_.find(tri[2]) != vertex_map_.end());
      SL_ENSURE("Inserted", edge_map_.find(tri.undirected_edge(0)) != edge_map_.end());
      SL_ENSURE("Inserted", edge_map_.find(tri.undirected_edge(1)) != edge_map_.end());
      SL_ENSURE("Inserted", edge_map_.find(tri.undirected_edge(2)) != edge_map_.end());
    }
    
    /// Insert triangle, also inserting new edges and vertices if not already present
    void insert_triangle (vertex_t i0, vertex_t i1, vertex_t i2) {
      insert_triangle(triangle_t(i0,i1,i2));      
    }

    /// Erase triangle, also erasing isolated edges and vertices
    virtual void erase_triangle (const triangle_t& tri) {
      triangle_iterator_t t_it = triangle_map_.find(tri);
      if (t_it!=triangle_map_.end()) {
        edge_t e01(tri[0],tri[1]), e12(tri[1],tri[2]), e20(tri[2],tri[0]);

        vertex_iterator_t v_it0 = vertex_map_.find(tri[0]);
        SL_CHECK("Found",  v_it0 != vertex_map_.end() );
        vertex_iterator_t v_it1 = vertex_map_.find(tri[1]);
        SL_CHECK("Found",  v_it1 != vertex_map_.end() );
        vertex_iterator_t v_it2 = vertex_map_.find(tri[2]);
        SL_CHECK("Found",  v_it2 != vertex_map_.end() );

        edge_iterator_t e_it01 = edge_map_.find(e01);
        SL_CHECK("Found",  e_it01 != edge_map_.end() );
        edge_iterator_t e_it12 = edge_map_.find(e12);
        SL_CHECK("Found",  e_it12 != edge_map_.end() );
        edge_iterator_t e_it20 = edge_map_.find(e20);
        SL_CHECK("Found",  e_it20 != edge_map_.end() );

        // Update edges
        e_it01->second.triangle_set().erase(tri);
        e_it12->second.triangle_set().erase(tri);
        e_it20->second.triangle_set().erase(tri);

        // update vertices
        v_it0->second.triangle_set().erase(tri);
        v_it1->second.triangle_set().erase(tri);
        v_it2->second.triangle_set().erase(tri);

        // Erase unreferenced edges
	bool erase_e01 = e_it01->second.triangle_set().empty();
	bool erase_e12 = e_it12->second.triangle_set().empty();
	bool erase_e20 = e_it20->second.triangle_set().empty();
	
        if (erase_e01) {
          v_it0->second.edge_set().erase(e01);
          v_it1->second.edge_set().erase(e01);
          edge_map_.erase(e01); // do not use iterators as they might be invalidated
	}

        if (erase_e12) {
          v_it1->second.edge_set().erase(e12);
          v_it2->second.edge_set().erase(e12);
          edge_map_.erase(e12); // do not use iterators as they might be invalidated
	}

        if (erase_e20) {
          v_it2->second.edge_set().erase(e20);
          v_it0->second.edge_set().erase(e20);
          edge_map_.erase(e20); // do not use iterators as they might be invalidated
        }

        // Erase unreferenced vertices
	bool erase_v0 = v_it0->second.edge_set().empty();
	bool erase_v1 = v_it1->second.edge_set().empty();
	bool erase_v2 = v_it2->second.edge_set().empty();
	  
        if (erase_v0) {
          // SL_CHECK("Consistent", v_it0->second.triangle_set().empty());
          vertex_map_.erase(tri[0]); // do not use iterators as they might be invalidated
        }
        if (erase_v1) {
          // SL_CHECK("Consistent", v_it1->second.triangle_set().empty());
          vertex_map_.erase(tri[1]); // do not use iterators as they might be invalidated
        }
        if (erase_v2) {
          // SL_CHECK("Consistent", v_it2->second.triangle_set().empty());
          vertex_map_.erase(tri[2]); // do not use iterators as they might be invalidated
        }

        // Erase triangle
        triangle_map_.erase(t_it);
      }
      SL_CHECK("Erased", !has_triangle(tri));
    }

    /// Erase triangle, also erasing isolated edges and vertices
    void erase_triangle (vertex_t i0, vertex_t i1, vertex_t i2) {
      erase_triangle(triangle_t(i0,i1,i2));
    }
    
    /// Erase all triangles
    virtual void clear() {
      vertex_map_.clear();
      edge_map_.clear();
      triangle_map_.clear();
    }

  public: // Edge removal

    /// Remove all triangles sharing the given edge
    void erase_edge(const edge_t& e) {
      const small_triangle_set_t* e_tri_ptr = edge_triangles(e);
      if (e_tri_ptr) {
        small_triangle_set_t e_tri = (*e_tri_ptr); // note local copy
        for (small_triangle_set_t::iterator t_it = e_tri.begin();
             t_it != e_tri.end();
             ++t_it) {
          erase_triangle(*t_it);
        }
      }
      SL_ENSURE("Erased", !has_edge(e));
    }

    /// Remove all triangles sharing the given edge
    void erase_edge(vertex_t i0, vertex_t i1) {
      erase_edge(edge_t(i0,i1));
    }

  public: // Vertex removal

    /// Remove all triangles sharing the given vertex
    void erase_vertex(const vertex_t i0) {
      const small_triangle_set_t* v_tri_ptr = vertex_triangles(i0);
      if (v_tri_ptr) {
        small_triangle_set_t v_tri = (*v_tri_ptr); // note local copy
        for (small_triangle_set_t::iterator t_it = v_tri.begin();
             t_it != v_tri.end();
             ++t_it) {
          erase_triangle(*t_it);
        }
      }
      SL_ENSURE("Erased", !has_vertex(i0));
    }
    
  public: // Collapse operations

    /// Replace src_vidx with dst_vidx, removing all degenerate triangles
    void collapse_edge(vertex_t src_vidx, vertex_t dst_vidx) {
      SL_REQUIRE("Exists", has_vertex(src_vidx));
      SL_REQUIRE("Exists", has_vertex(dst_vidx));

      if (src_vidx != dst_vidx) {
        // Get neighborhood of removed vertex
        small_triangle_set_t src_vidx_tri = *vertex_triangles(src_vidx);
        
        // Insert new triangles
        for (small_triangle_set_t::iterator t_it = src_vidx_tri.begin();
             t_it != src_vidx_tri.end();
             ++t_it) {
          const triangle_t& old_tri = *t_it;

          triangle_t new_tri = old_tri.remapped(src_vidx,dst_vidx);
          if (new_tri.is_valid()) {
            SL_CHECK("Exists", has_vertex(new_tri[0]));
            SL_CHECK("Exists", has_vertex(new_tri[1]));
            SL_CHECK("Exists", has_vertex(new_tri[2]));
            
            // Insert new triangle
            insert_triangle(new_tri);

            // Copy triangle data
            const triangle_data_t* old_tridata_ptr = triangle_data(old_tri);
            triangle_data_t* new_tridata_ptr = triangle_data(new_tri);
            SL_CHECK("Exists", old_tridata_ptr);
            SL_CHECK("Exists", new_tridata_ptr);
            *new_tridata_ptr = *old_tridata_ptr;

            // Copy edge data
            for (std::size_t i=0; i<3; ++i) {
              edge_t old_edge = old_tri.undirected_edge(i);
              const edge_data_t* old_edgedata_ptr = edge_data(old_edge);
              edge_t new_edge = old_edge.remapped(src_vidx,dst_vidx);
              edge_data_t* new_edgedata_ptr = edge_data(new_edge);
              SL_CHECK("Exists", old_edgedata_ptr);
              SL_CHECK("Exists", new_edgedata_ptr);
              *new_edgedata_ptr = *old_edgedata_ptr;
            }

            // Vertex data remains the same...
          }
        }

        // Erase old triangles
        for (small_triangle_set_t::iterator t_it = src_vidx_tri.begin();
             t_it != src_vidx_tri.end();
             ++t_it) {
          const triangle_t& old_tri = *t_it;
          erase_triangle(old_tri);
        }
      }
      
      SL_ENSURE("Erased", (src_vidx == dst_vidx) || !has_vertex(src_vidx));
    }

    /// Replace v1 and v2 with v3, removing all degenerate triangles
    void collapse_triangle(vertex_t v1, vertex_t v2, vertex_t v3) {
      collapse_edge(v1, v3);
      if (has_vertex(v2)) collapse_edge(v2, v3);
    }

  public: // split operations
    
    void split_edge(edge_t e,
                    vertex_t vmid,
                    const vertex_data_t& vmid_data) {
      SL_REQUIRE("Edge exists", has_edge(e));

      // Refine triangles incident on e
      small_triangle_set_t e_triangles = *(edge_triangles(e));
      for (small_triangle_set_t::iterator it = e_triangles.begin();
           it != e_triangles.end();
           ++it) {
        const triangle_t& old_tri = *it;
        if        (old_tri.undirected_edge(0) == e) {
          triangle_t t1 = triangle_t(old_tri[0], vmid, old_tri[2]);
          triangle_t t2 = triangle_t(vmid, old_tri[1], old_tri[2]);
          insert_triangle(t1);
          insert_triangle(t2);
          *triangle_data(t1) = *triangle_data(old_tri); 
          *triangle_data(t2) = *triangle_data(old_tri);
          // FIXME ??? Edge data
        } else if (old_tri.undirected_edge(1) == e) {
          triangle_t t1 = triangle_t(old_tri[1], vmid, old_tri[0]);
          triangle_t t2 = triangle_t(vmid, old_tri[2], old_tri[0]);
          insert_triangle(t1);
          insert_triangle(t2);
          *triangle_data(t1) = *triangle_data(old_tri); 
          *triangle_data(t2) = *triangle_data(old_tri);
          // FIXME ??? Edge data
        } else if (old_tri.undirected_edge(2) == e) {
          triangle_t t1 = triangle_t(old_tri[2], vmid, old_tri[1]);
          triangle_t t2 = triangle_t(vmid, old_tri[0], old_tri[1]);
          insert_triangle(t1);
          insert_triangle(t2);
          *triangle_data(t1) = *triangle_data(old_tri); 
          *triangle_data(t2) = *triangle_data(old_tri);
          // FIXME ??? Edge data
        }
      }

      // Update vdata
      vertex_data_t* vdata = vertex_data(vmid);
      if (vdata) {
        *vdata = vmid_data;
      }
      
      // Erase triangles incident on e
      for (small_triangle_set_t::iterator t_it = e_triangles.begin();
           t_it != e_triangles.end();
           ++t_it) {
        const triangle_t& old_tri = *t_it;
        erase_triangle(old_tri);
      }

      // FIXME edge data???
      // FIXME triangle data???
    }
    
  public: // vertex attributes
    
    const vertex_map_t& vertex_map() const {
      return vertex_map_;
    }
    
    const_vertex_iterator_t vertex_begin() const {
      return vertex_map_.begin();
    }

    const_vertex_iterator_t vertex_end() const {
      return vertex_map_.end();
    }

    vertex_iterator_t vertex_begin() {
      return vertex_map_.begin();
    }

    vertex_iterator_t vertex_end() {
      return vertex_map_.end();
    }


    bool has_vertex(vertex_t i) const {
      return vertex_map_.find(i) != vertex_map_.end();
    }

    vertex_t new_vertex_id() const {
      // FIXME: Pretty baroque
      vertex_t result = vertex_t(irng_.value());
      while (has_vertex(result)) {
        result = vertex_t(irng_.value());
      }
      return result;
    }
    
    void vertices_in(std::set<vertex_t>& s) const {
      for (const_vertex_iterator_t v_it = vertex_map_.begin(); v_it != vertex_map_.end(); ++v_it) {
        s.insert(v_it->first);
      }
    }

    const vertex_data_t* vertex_data(vertex_t i) const {
      const_vertex_iterator_t v_it = vertex_map_.find(i);
      return (v_it != vertex_map_.end() ? v_it->second.data() : 0);
    }

    vertex_data_t* vertex_data(vertex_t i) {
      vertex_iterator_t v_it = vertex_map_.find(i);
      return (v_it != vertex_map_.end() ? v_it->second.data() : 0);
    }
    
    const small_edge_set_t* vertex_edges(vertex_t i) const {
      const_vertex_iterator_t v_it = vertex_map_.find(i);
      return ( v_it != vertex_map_.end() ? &v_it->second.edge_set() : 0 );
    }

    std::size_t vertex_edge_count(vertex_t i) const {
      const_vertex_iterator_t v_it = vertex_map_.find(i);
      return ( v_it != vertex_map_.end() ? v_it->second.edge_set().size() : 0 );
    }
    
    const small_triangle_set_t* vertex_triangles(vertex_t i) const {
      const_vertex_iterator_t v_it = vertex_map_.find(i);
      return ( v_it != vertex_map_.end() ? &v_it->second.triangle_set() : 0 );
    }

    std::size_t vertex_triangle_count(vertex_t i) const {
      const_vertex_iterator_t v_it = vertex_map_.find(i);
      return ( v_it != vertex_map_.end() ? v_it->second.triangle_set().size() : 0 );
    }
    
  public: // edge attributes
    
    const edge_map_t& edge_map() const {
      return edge_map_;
    }

    const_edge_iterator_t edge_begin() const {
      return edge_map_.begin();
    }

    const_edge_iterator_t edge_end() const {
      return edge_map_.end();
    }

    edge_iterator_t edge_begin() {
      return edge_map_.begin();
    }

    edge_iterator_t edge_end() {
      return edge_map_.end();
    }

    
    bool has_edge(edge_t i) const {
      return edge_map_.find(i) != edge_map_.end();
    }

    bool has_edge(vertex_t i0, vertex_t i1) const {
      return has_edge(edge_t(i0,i1));
    }

    void edges_in(std::set<edge_t>& s) const {
      for (const_edge_iterator_t v_it = edge_map_.begin(); v_it != edge_map_.end(); ++v_it) {
        s.insert(v_it->first);
      }
    }

#if HAVE_STLEXT_UNORDERED_CONTAINERS
    void edges_in(stlext::unordered_set<edge_t, sl::hash<edge_t> >& s) const {
      for (const_edge_iterator_t v_it = edge_map_.begin(); v_it != edge_map_.end(); ++v_it) {
        s.insert(v_it->first);
      }
    }
#endif
    
    const edge_data_t* edge_data(edge_t i) const {
      const_edge_iterator_t v_it = edge_map_.find(i);
      return (v_it != edge_map_.end() ? v_it->second.data() : 0);
    }

    edge_data_t* edge_data(edge_t i) {
      edge_iterator_t v_it = edge_map_.find(i);
      return (v_it != edge_map_.end() ? v_it->second.data() : 0);
    }
        
    const small_triangle_set_t* edge_triangles(edge_t i) const {
      const_edge_iterator_t v_it = edge_map_.find(i);
      return ( v_it != edge_map_.end() ? &v_it->second.triangle_set() : 0 );
    }

    std::size_t edge_triangle_count(edge_t i) const {
      const_edge_iterator_t v_it = edge_map_.find(i);
      return ( v_it != edge_map_.end() ? v_it->second.triangle_set().size() : 0 );
    }

  public: // triangle attributes

    const triangle_map_t& triangle_map() const {
      return triangle_map_;
    }

    const_triangle_iterator_t triangle_begin() const {
      return triangle_map_.begin();
    }

    const_triangle_iterator_t triangle_end() const {
      return triangle_map_.end();
    }

    triangle_iterator_t triangle_begin() {
      return triangle_map_.begin();
    }

    triangle_iterator_t triangle_end() {
      return triangle_map_.end();
    }

    
    bool has_triangle(triangle_t i) const {
      return triangle_map_.find(i) != triangle_map_.end();
    }

    bool has_triangle(vertex_t i0, vertex_t i1, vertex_t i2) const {
      return has_triangle(triangle_t(i0,i1,i2));
    }

    void triangles_in(std::set<triangle_t>& s) const {
      for (const_triangle_iterator_t t_it = triangle_map_.begin(); t_it != triangle_map_.end(); ++t_it) {
        s.insert(t_it->first);
      }
    }
    
#if HAVE_STLEXT_UNORDERED_CONTAINERS
    void triangles_in(stlext::unordered_set<triangle_t, sl::hash<triangle_t> >& s) const {
      for (const_triangle_iterator_t t_it = triangle_map_.begin(); t_it != triangle_map_.end(); ++t_it) {
        s.insert(t_it->first);
      }
    }
#endif
    
    const triangle_data_t* triangle_data(triangle_t i) const {
      const_triangle_iterator_t t_it = triangle_map_.find(i);
      return (t_it != triangle_map_.end() ? t_it->second.data() : 0);
    }

    triangle_data_t* triangle_data(triangle_t i) {
      triangle_iterator_t t_it = triangle_map_.find(i);
      return (t_it != triangle_map_.end() ? t_it->second.data() : 0);
    }

  public: // Queries

    /**
     *  The mesh is manifold if each edge has at most two adjacent triangles.
     *  It is possible that the mesh has multiple connected components.
     */
    bool is_manifold () const {
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        if (e_it->second.triangle_set().size()> 2) {
          return false;
        }
      }
      return true;
    }

    /**
     *  The mesh is a consistently oriented manifold if each edge has at most two adjacent triangles,
     *  that traverse it in two different directions.
     *  It is possible that the mesh has multiple connected components.
     */
    bool is_consistently_oriented_manifold () const {
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        const edge_t& e = e_it->first;
        const small_triangle_set_t& e_tri = e_it->second.triangle_set();
        std::size_t N = e_tri.size();
        if (N>2) {
          return false;
        } if (N==2) {
          // Check orientation
          triangle_t t0 = *(e_tri.begin());
          triangle_t t1 = *(++(e_tri.begin()));
          std::size_t i0 = 3;
          if (e == t0.undirected_edge(0)) i0 = 0;
          if (e == t0.undirected_edge(1)) i0 = 1;
          if (e == t0.undirected_edge(2)) i0 = 2;
          SL_CHECK("Edge exists", i0 < 3);
          half_edge_t he = t0.directed_edge(i0);
          if (he == t1.directed_edge(0)) return false;
          if (he == t1.directed_edge(1)) return false;
          if (he == t1.directed_edge(2)) return false;
        }
      }
      return true;
    }

    /**
     *  The mesh is closed if each edge has exactly two adjacent triangles.
     *  It is possible that the mesh has multiple connected components.
     */
    bool is_closed () const {
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        if (e_it->second.triangle_set().size() != 2) {
          return false;
        }
      }
      return true;
    }

    /**
     *  The mesh is connected if each triangle can be reached from any other
     *  triangle by a single traversal.
     */
    bool is_connected () const {
      std::size_t tcount = triangle_map_.size();
      if (tcount == 0) {
        return true;
      } else {
        typename stl_set_selector<triangle_t>::unordered_t visited_triangles;
        std::stack<triangle_t> triangle_stack;
        triangle_stack.push(triangle_map_.begin()->first);
        visited_triangles.insert(triangle_map_.begin()->first);
        while (!triangle_stack.empty()) {
          triangle_t tri = triangle_stack.top();
          triangle_stack.pop();
          // Visit unvisited adjacent triangles
          for (std::size_t i=0; i<3; i++) {
            const_edge_iterator_t e_it = edge_map_.find(edge_t(tri[i],tri[(i+1)%3]));
            SL_CHECK("Found", e_it != edge_map_.end());
            const small_triangle_set_t& triset = e_it->second.triangle_set();
            for (small_triangle_set_t::const_iterator t_it = triset.begin();
                 t_it != triset.end();
                 ++t_it) {
              const triangle_t& tri_adj = *(t_it);
              if (visited_triangles.find(tri_adj) == visited_triangles.end()) {
                triangle_stack.push(tri_adj);
                visited_triangles.insert(tri_adj);
              }
            }
          }
        }
        return tcount == visited_triangles.size();
      }
    }

  public: // Boundaries

    bool is_boundary_edge(const edge_t& e) const {
      const small_triangle_set_t* triset_ptr = edge_triangles(e);
      return triset_ptr && (triset_ptr->size() < 2);
    }
    
    std::size_t boundary_edge_count() const {
      std::size_t result = 0;
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        if (e_it->second.triangle_set().size() < 2) {
          ++result;
        }
      }
      return result;
    }

    bool is_boundary_vertex(const vertex_t& i) const {
      bool result = false;
      const_vertex_iterator_t v_it = vertex_map_.find(i);
      if (v_it != vertex_map_.end()) {
        const small_edge_set_t& eset = v_it->second.edge_set();
        for (small_edge_set_t::const_iterator e_it = eset.begin();
             !result && e_it != eset.end();
             ++e_it) {
          result = is_boundary_edge(*e_it);
        }
      }
      return result;
    }

    void boundary_edges_in(std::set<edge_t>& s) const {
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        if (e_it->second.triangle_set().size() < 2) {
          s.insert(e_it->first);
        }
      }
    }

    std::size_t boundary_vertex_count() const {
      std::size_t result = 0;
      for (const_vertex_iterator_t v_it = vertex_map_.begin();
           v_it != vertex_map_.end();
           ++v_it) {
        bool is_part_of_boundary_edge = false;
        const small_edge_set_t& eset = *vertex_edges(v_it->first);
        for (small_edge_set_t::const_iterator e_it = eset.begin();
             e_it != eset.end() && !is_part_of_boundary_edge;
             ++e_it) {
          is_part_of_boundary_edge = is_boundary_edge(*e_it);
        }
        if (is_part_of_boundary_edge) ++result;
      }
      return result;
    }

    void boundary_vertices_in(std::set<vertex_t>& s) const {
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        if (e_it->second.triangle_set().size() < 2) {
          s.insert(e_it->first[0]);
          s.insert(e_it->first[1]);
        }
      }
    }

#if HAVE_STLEXT_UNORDERED_CONTAINERS
    void boundary_vertices_in(stlext::unordered_set<vertex_t, sl::hash<vertex_t> >& s) const {
      for (const_edge_iterator_t e_it = edge_map_.begin(); e_it != edge_map_.end(); ++e_it) {
        if (e_it->second.triangle_set().size() < 2) {
          s.insert(e_it->first[0]);
          s.insert(e_it->first[1]);
        }
      }
    }
#endif
    
  public: // Related meshes
    
    /// The mesh with all triangles in reversed order
    this_t* reversed_order_mesh(bool copy_vertex_attributes = true,
                                bool copy_edge_attributes = true,
                                bool copy_triangle_attributes = true) const {
      this_t* result = new_mesh();
      for (const_triangle_iterator_t t_it = triangle_map_.begin();
           t_it != triangle_map_.end();
           ++t_it) {
        result->insert_triangle(t_it->first);
        if (copy_triangle_attributes) {
          triangle_data_t* ptr = result->triangle_data(t_it->first);
          if (ptr && t_it->second.data()) {
            *ptr = *(t_it->second.data());
          }
        }
      }
      if (copy_edge_attributes) {
        for (const_edge_iterator_t e_it = edge_map_.begin();
             e_it != edge_map_.end();
             ++e_it) {
          edge_data_t* ptr = result->edge_data(e_it->first);
          if (ptr && e_it->second.data()) {
            *ptr = *(e_it->second.data());
          }
        }
      }
      if (copy_vertex_attributes) {
        for (const_vertex_iterator_t v_it = vertex_map_.begin();
             v_it != vertex_map_.end();
             ++v_it) {
          vertex_data_t* ptr = result->vertex_data(v_it->first);
          if (ptr && v_it->second.data()) {
            *ptr = *(v_it->second.data());
          }
        }
      }
      
      return result;
    }
  };
} // namespace sl

#endif
