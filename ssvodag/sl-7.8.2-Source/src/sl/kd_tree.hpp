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
#ifndef SL_KD_TREE_HPP
#define SL_KD_TREE_HPP

/////////////////////////////// DEPRECATED CLASS ////////////////////////
#warning \
Class sl::kd_tree from sl/kd_tree.hpp is obsolete and will be \
removed from future sl releases. Please use class sl::kdtree \
defined in sl/kdtree.hpp.
/////////////////////////////// DEPRECATED CLASS ////////////////////////

#include "sl/assert.hpp"
#include "sl/fixed_size_point.hpp"
#include "sl/utility.hpp"
#include <sl/memory_pool.hpp>
#include <map>

namespace sl {

  namespace detail {

    /**
     *  A simple object storing K coordinates
     */
    template <std::size_t G_dimension, 
      class G_value>
    class kd_coords {
    public:
      enum { dimension = G_dimension };

      typedef kd_coords<G_dimension, G_value> this_t;
      typedef G_value value_t;
    protected:
      value_t c_[dimension];
    public:
      inline kd_coords() {
      }

      inline const value_t &operator[](std::size_t i) const {
        SL_REQUIRE("Good index", i<dimension);
        return c_[i];
      }
      inline value_t &operator[](std::size_t i) {
        SL_REQUIRE("Good index", i<dimension);
        return c_[i];
      }
    };
       
    /**
     *  A simple predicate always returning true
     */
    template <class T>
    class kd_true {
    public:
      inline kd_true() {}
      inline bool operator()(const T&) const { return true; }
    };
    
    /**
     *  A node in a kd_tree
     */
    template <std::size_t G_dimension, 
      class G_value,
      class G_data>
    class kd_node {
    public:
      enum { dimension = G_dimension };

      typedef kd_node<G_dimension, G_value, G_data> this_t;
      typedef G_value value_t;
      typedef G_data  data_t;
    protected:
      kd_node* lo_child_;
      kd_node* hi_child_;
      std::size_t discriminant_;
      std::size_t tree_size_;
      data_t data_;
    public:

      /// Init an empty kd_node
      inline kd_node()
          :
          lo_child_(0), hi_child_(0), discriminant_(0), tree_size_(1), data_(data_t()) {
      }
	
      /// Init a kd_node containing x
      inline kd_node(const data_t& x) 
          :
          lo_child_(0), hi_child_(0), discriminant_(0), tree_size_(1), data_(x) {
      }

      /// The number of nodes in the tree rooted at this (including this)
      inline std::size_t tree_size() const {
        return tree_size_;
      }

      /// Set to x the number of nodes in the tree rooted at this (including this)
      inline void set_tree_size(std::size_t x) {
        SL_REQUIRE("Good tree size", x >= 1);
        tree_size_ = x;
      }
      
      /// The data associated to this node
      inline const data_t& data() const {
        return data_;
      }

      /// The data associated to this node
      inline data_t& data() {
        return data_;
      }

      /// Is this a leaf node?
      inline bool is_leaf() const {
        return !lo_child_ && !hi_child_;
      }

      /// The child on the negative side of the discriminating discriminant
      inline const kd_node* lo_child() const {
        return lo_child_;
      }

      /// The child on the negative side of the discriminating discriminant
      inline kd_node* lo_child() {
        return lo_child_;
      }

      /// Set to x the child on the negative side of the discriminating discriminant
      inline void set_lo_child(kd_node* x) {
        lo_child_ = x;
      }

      /// The child on the positive side of the discriminating discriminant
      inline const kd_node* hi_child() const {
        return hi_child_;
      }

      /// The child on the positive side of the discriminating discriminant
      inline kd_node* hi_child() {
        return hi_child_;
      }

      /// Set to x the child on the positive side of the discriminating discriminant
      inline void set_hi_child(kd_node* x) {
        hi_child_ = x;
      }

      /// The split discriminant associated to this node
      inline std::size_t discriminant() const {
        return discriminant_;
      }

      // The constant of the discriminant
      inline value_t discriminant_plane_constant() const {
        return data_[discriminant()];
      }
      
      /// Set the spltting discriminant associated to this node to x
      inline void set_discriminant(std::size_t x) {
        SL_REQUIRE("Good discriminant", x<dimension);
        SL_REQUIRE("Is leaf", is_leaf());
        discriminant_ = x;
      }

      /// Apply action a to the tree rooted at this
      template <class action_t>
      inline void apply(action_t& a) const {
        a.apply(data());
        if (lo_child()) lo_child()->apply(a);
        if (hi_child()) hi_child()->apply(a);
      }
 
    };

  }

  /**
   *  kd-trees of dimension G_dimension containing
   *  entities of type G_data with coordinates of
   *  type G_value. The class G_data
   *  is constrained to export G_value operator[](std::size_t i)
   *  to provide read access to coordinates 0 to dimension.
   */
  template <std::size_t G_dimension, 
    class G_value,
    class G_data>
  class kd_tree {
  public:
    enum { dimension = G_dimension };

    typedef kd_tree<G_dimension, G_value, G_data> this_t;
    typedef G_value value_t;
    typedef G_data  data_t;

  protected:

    typedef detail::kd_coords<G_dimension, G_value> kd_coords_t;
    typedef detail::kd_node<G_dimension, G_value, G_data> node_t;

    // Root of the tree
    node_t* root_;

    // Number of points
    std::size_t count_;

    // Bounding box
    kd_coords_t bbox_lo_;
    kd_coords_t bbox_hi_;

    /// The pool of nodes
    memory_pool* node_pool_;

    /// Allocate a new node from the pool associated to this tree
    inline node_t* pool_new_node(const data_t& d) {
      SL_REQUIRE("Pool exists", node_pool_);
      node_t *const result = static_cast<node_t*>(node_pool_->allocate());
      if (result) { 
        try { new (result) node_t(d); }
        catch (...) { node_pool_->release(result); throw; }
      }
      return result;
    }
      
    /// Delete the node from the pool associated to this tree
    inline void pool_delete_node(node_t* const node_ptr) {
      SL_REQUIRE("Node exists", node_ptr);
      SL_REQUIRE("Pool exists", node_pool_);
      node_ptr->~node_t();
      node_pool_->release(node_ptr);
    }
      
  public:

    /// Create an empty kd_tree
    inline kd_tree(const std::size_t memory_pool_chunk_first_count = 0,
                   const std::size_t memory_pool_chunk_grow_factor = 1) 
        : 
        root_(0), count_(0), node_pool_(new memory_pool("kd_node", sizeof(node_t),memory_pool_chunk_first_count,memory_pool_chunk_grow_factor))
    {
      for (std::size_t i=0; i<dimension; ++i) {
        bbox_lo_[i] = static_cast<value_t>(0.0);
        bbox_hi_[i] = static_cast<value_t>(0.0);
      }
    }

    /// Destroy a kd_tree
    virtual inline ~kd_tree() {
      clear();
    }

    /// The number of data points
    inline std::size_t count() const {
      return count_;
    }

    /// Is the tree empty
    inline bool is_empty() const {
      return count() == 0;
    }

    /// Clear the tree
    virtual void clear();

    /// Insert a new data point in the tree. The point is copied.
    virtual void insert(const data_t& p);

    /// Apply action a to each datum in arbitrary order
    template <class action_t>
    inline void apply(action_t& a) const {
      if (root_) root_->apply(a);
    }
      
    /**
     *  Initialize knn with a maximum of k points that are closest to the 
     *  query point p and are in the sphere centered at p with radius sqrt(cut_r2).
     *  The key of the map knn is the square of the Euclidean distance to the
     *  query point, the value is a pointer to the data point.
     */
    template <class point_t>
    inline void nearest_neighbors_in(std::multimap<value_t, data_t*>& knn,
                                     const point_t& p,
                                     std::size_t k,
                                     value_t cut_r2) const {
      detail::kd_true<data_t> q_true;
      nearest_neighbors_in(knn,p,k,cut_r2,q_true);
    }

    /**
     *  The closest neighbor of p within a distance sqrt(cut_r2) from it.
     *  NULL if not existent.
     */
    template <class point_t>
    const data_t *nearest_neighbor(const point_t& p,
                                   const value_t& cut_r2) const {
      data_t* result = 0;
      std::multimap<value_t, data_t*> knn;
      nearest_neighbors_in(knn, p, 1, cut_r2);
      if (!knn.empty()) {
        result = knn.begin()->second;
      }
      return result;
    }

    /**
     *  The closest neighbor of p within a distance sqrt(cut_r2) from it.
     *  NULL if not existent.
     */
    template <class point_t>
    data_t *nearest_neighbor(const point_t& p,
                             const value_t& cut_r2) {
      data_t* result = 0;
      std::multimap<value_t, data_t*> knn;
      nearest_neighbors_in(knn, p, 1, cut_r2);
      if (!knn.empty()) {
        result = knn.begin()->second;
      }
      return result;
    }

    /**
     *  Initialize knn with a maximum of k points that are closest to the 
     *  query point p, are in the sphere centered at p with radius sqrt(cut_r2),
     *  and match the additional criterion q.
     *  The key of the map knn is the square of the Euclidean distance to the
     *  query point, the value is a pointer to the data point.
     */
    template <class point_t, class predicate_t>
    void nearest_neighbors_in(std::multimap<value_t, data_t*>& knn,
                              const point_t& p,
                              std::size_t k,
                              value_t cut_r2,
                              const predicate_t& q) const;

    /**
     *  The closest neighbor of p that meets criterion q and lies
     *  within a distance sqrt(cut_r2) from p.
     *  NULL if not existent.
     */
    template <class point_t, class predicate_t>
    const data_t *nearest_neighbor(const point_t& p,
                                   const value_t& cut_r2,
                                   const predicate_t& q) const {
      data_t* result = 0;
      std::multimap<value_t, data_t*> knn;
      nearest_neighbors_in(knn, p, 1, cut_r2, q);
      if (!knn.empty()) {
        result = knn.begin()->second;
      }
      return result;
    }

    /**
     *  The closest neighbor of p that meets criterion q and lies
     *  within a distance sqrt(cut_r2) from p.
     *  NULL if not existent.
     */
    template <class point_t, class predicate_t>
    data_t *nearest_neighbor(const point_t& p,
                             const value_t& cut_r2,
                             const predicate_t& q) {
      data_t* result = 0;
      std::multimap<value_t, data_t*> knn;
      nearest_neighbors_in(knn, p, 1, cut_r2, q);
      if (!knn.empty()) {
        result = knn.begin()->second;
      }
      return result;
    }

    /**
     *  Set min_depth and max_depth to the shortest and longest distance from
     *  the root to the leafs
     */
    void tree_depth_in(std::size_t* min_depth,
                       std::size_t* max_depth) const;

    /**
     *  True iff the invariant is respected. *Very* slow, use it for
     *  debugging purposes.
     */
    virtual bool invariant() const;

    /**
     *  The sum of the depths of all internal nodes
     */
    std::size_t internal_path_length() const;

    /**
     *  The sum of the depths of all leaf nodes
     */
    std::size_t external_path_length() const;
    

  protected:
    
    void delete_nodes(node_t* n);

    template <class point_t, class predicate_t>
    void nearest_neighbors_update_in(std::multimap<value_t, data_t*>& knn,
                                     const node_t* node,
                                     const point_t& p,
                                     std::size_t k,
                                     value_t box_dist2,
                                     kd_coords_t& box_offset,
                                     value_t& cut_r2,
                                     const predicate_t& q) const;

    void tree_depth_update_in(const node_t* node,
                              std::size_t node_depth,
                              std::size_t* min_depth,
                              std::size_t* max_depth) const;

    virtual bool node_invariant(const node_t* node) const;

    bool is_lo(const node_t* node,
               std::size_t i,
               value_t x_i) const;

    bool is_hi(const node_t* node,
               std::size_t i,
               value_t x_i) const;

    std::size_t external_path_length(const node_t* node,
                                     std::size_t depth) const;

    std::size_t internal_path_length(const node_t* node,
                                     std::size_t depth) const;

    
    void split(node_t* t,
               std::size_t i,
               value_t x_i,
               node_t** t_lo,
               node_t** t_hi);

    node_t *join(node_t* a,
                 node_t* b,
                 std::size_t i);

  }; // class kd_tree
    			 
} // namespace sl

// -- Inline implementation

namespace sl {

  template <std::size_t G_dimension, class G_value,class G_data>
  void kd_tree<G_dimension, G_value, G_data>::delete_nodes(node_t* n) {
    if (n) {
      delete_nodes(n->lo_child());
      n->set_lo_child(0);
      delete_nodes(n->hi_child());
      n->set_hi_child(0);
      pool_delete_node(n);
      --count_;
    }
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  void kd_tree<G_dimension, G_value, G_data>::clear() {
    delete_nodes(root_);
    root_ = 0;
    SL_ENSURE("Empty", is_empty());
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  void kd_tree<G_dimension, G_value, G_data>::insert(const data_t& p) {
    node_t * new_node = pool_new_node(p);
    ++count_;

    if (!root_) {
      // First node
      root_ = new_node;
      // Init bbox
      for (std::size_t i = 0; i<dimension; ++i) {
        bbox_lo_[i] = p[i];
        bbox_hi_[i] = p[i];
      }
    } else {
      // Update bbox
      for (std::size_t i = 0; i<dimension; ++i) {
        if (p[i]<bbox_lo_[i]) bbox_lo_[i] = p[i];
        if (p[i]>bbox_hi_[i]) bbox_hi_[i] = p[i];
      }
#define KD_SQUARISH
#ifdef KD_SQUARISH
      kd_coords_t parent_box_lo = bbox_lo_;
      kd_coords_t parent_box_hi = bbox_hi_;
#endif

      // Find insertion point and increment tree size along the path
      node_t* parent = 0;
      node_t* leaf   = root_;
      while (leaf) {
        parent = leaf;

        parent->set_tree_size(parent->tree_size() + 1);

        std::size_t d = parent->discriminant();
        if (p[d] <= parent->data()[d]) {
#ifdef KD_SQUARISH
          parent_box_hi[d] = parent->data()[d];
#endif
          leaf = parent->lo_child();
        } else {
#ifdef KD_SQUARISH
          parent_box_lo[d] = parent->data()[d];
#endif
          leaf = parent->hi_child();
        }
      }
      SL_CHECK("Parent exists", parent);
      
      // If parent has a child, the new is the other one
      if (parent->lo_child()) {
        SL_CHECK("Coherent", !parent->hi_child());
        SL_CHECK("Ordered", new_node->data()[parent->discriminant()] > parent->data()[parent->discriminant()]);
        parent->set_hi_child(new_node);
      } else if (parent->hi_child()) {
        SL_CHECK("Coherent", !parent->lo_child());
        SL_CHECK("Ordered", new_node->data()[parent->discriminant()] <= parent->data()[parent->discriminant()]);
        parent->set_lo_child(new_node);
      } else {
        // Choose an appropriate discriminator for the parent
        std::size_t d_best_discriminant = 0;
#ifdef KD_SQUARISH
        value_t d_best_delta = parent_box_hi[0] - parent_box_lo[0];
        for (std::size_t i=1; i<dimension; ++i) {
          value_t d_ith_delta = parent_box_hi[i] - parent_box_lo[i];
          if (d_ith_delta > d_best_delta) {
            d_best_discriminant = i;
            d_best_delta = d_ith_delta;
          }
        }
#else
        value_t d_best_delta  = p[0] - parent->data()[0];
        value_t d_best_dist = (d_best_delta >= value_t(0.0)) ? (d_best_delta) : (-d_best_delta);
        for (std::size_t i=1; i<dimension; ++i) {
          value_t d_ith_delta  = p[i] - parent->data()[i];
          value_t d_ith_dist = (d_ith_delta >= value_t(0.0)) ? (d_ith_delta) : (-d_ith_delta);
          if (d_ith_dist > d_best_dist) {
            d_best_discriminant = i;
            d_best_delta = d_ith_delta;
            d_best_dist = d_ith_dist;
          }
        }
#endif
        // Set discriminator and insert
        parent->set_discriminant(d_best_discriminant);
        if (p[d_best_discriminant] <= parent->data()[d_best_discriminant]) {
          SL_CHECK("Ordered", new_node->data()[parent->discriminant()] <= parent->data()[parent->discriminant()]);
          parent->set_lo_child(new_node);
        } else {
          SL_CHECK("Ordered", new_node->data()[parent->discriminant()] > parent->data()[parent->discriminant()]);
          parent->set_hi_child(new_node);
        }
      }
    }
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  template <class point_t, class predicate_t>
  void kd_tree<G_dimension, G_value, G_data>::
  nearest_neighbors_in(std::multimap<value_t, data_t*>& knn,
                       const point_t& p,
                       std::size_t k,
                       value_t cut_r2,
                       const predicate_t& q) const {
    knn.clear();
    if (root_ && k>0 && cut_r2 >= 0.0f) {
      value_t cut_r2_updated = cut_r2;
      
      // Init squared distance to root bounding box
      value_t box_dist2 = value_t(0.0);
      kd_coords_t box_offset;
      for (std::size_t d=0; d<dimension; ++d) {
        box_offset[d] = value_t(0.0);
        value_t delta = p[d] - bbox_lo_[d];
        if (delta < 0.0) {
          box_offset[d]  = delta;
          box_dist2 += delta * delta;
        } else {
          delta = p[d] - bbox_hi_[d];
          if (delta > 0.0) {
            box_offset[d]  = delta;
            box_dist2 += delta * delta;
          }
        }
      }

      nearest_neighbors_update_in(knn, 
                                  root_, 
                                  p, 
                                  k, 
                                  box_dist2, 
                                  box_offset, 
                                  cut_r2_updated,
                                  q);
    }
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  template <class point_t, class predicate_t>
  void kd_tree<G_dimension, G_value, G_data>::
  nearest_neighbors_update_in(std::multimap<value_t, data_t*>& knn,
                              const node_t* node,
                              const point_t& p,
                              std::size_t k,
                              value_t box_dist2,
                              kd_coords_t &box_offset,
                              value_t& cut_r2,
                              const predicate_t& q) const {
    SL_REQUIRE("Node exists", node);
    SL_REQUIRE("Good k", k>0);
    SL_REQUIRE("Good radius squared", cut_r2>=0.0f);
    SL_REQUIRE("Good knn size", knn.size() <= k);

    // Sort children
    const node_t* near_node;
    const node_t* far_node;
    
    const std::size_t discriminant = node->discriminant();
    const value_t cut_diff = p[discriminant] - node->data()[discriminant];
    const value_t cut_diff2 = cut_diff*cut_diff;
    if (cut_diff < value_t(0.0)) {
      near_node = node->lo_child();
      far_node = node->hi_child();
    } else {
      near_node = node->hi_child();
      far_node = node->lo_child();
    }	

    // Always visit closer child first
    if (near_node) {
      nearest_neighbors_update_in(knn,
                                  near_node,
                                  p,
                                  k,
                                  box_dist2,
                                  box_offset,
                                  cut_r2,
                                  q);
    }

    // If cutting discriminant is closer than (updated)
    // cut distance, test data point and possibly far nodes
    if (cut_diff2 < cut_r2) {
      // Compute distance to query point
      value_t dist2 = cut_diff2;
      for (std::size_t i=1; i<dimension; ++i) {
        const std::size_t d = (discriminant+i)%dimension;
        const value_t delta = (node->data()[d] - p[d]);
        dist2 += delta * delta;
      }
      
      // Add to knn if closer than cut distance and matching criterion
      if (dist2 < cut_r2 && q(((node_t*)node)->data())) {
        knn.insert(std::pair<value_t, data_t*>(dist2, &(((node_t*)node)->data())));
        // Update distance bound if needed
        if (knn.size() == k) {
          cut_r2 = (*knn.rbegin()).first;
        } else if (knn.size() > k) {
          knn.erase(--knn.end());
          cut_r2 = (*knn.rbegin()).first;
        }
      }

      // Visit further child if close enough
      if (far_node) {
        // Incrementally update distance to bounding box
        const value_t old_off = box_offset[discriminant];
        box_dist2 += cut_diff2 - old_off*old_off;
        if (box_dist2 < cut_r2) {
          box_offset[discriminant] = cut_diff;
          nearest_neighbors_update_in(knn,
                                      far_node,
                                      p,
                                      k,
                                      box_dist2,
                                      box_offset,
                                      cut_r2,
                                      q);
          box_offset[discriminant] = old_off;
        }
      }
    }
  }	

  template <std::size_t G_dimension, class G_value,class G_data>
  void kd_tree<G_dimension, G_value, G_data>::tree_depth_update_in(const node_t* node,
                                                                   std::size_t node_depth,
                                                                   std::size_t* min_depth,
                                                                   std::size_t* max_depth) const {
    SL_REQUIRE("Node exists", node);
    SL_REQUIRE("Min depth exists", min_depth);
    SL_REQUIRE("Max depth exists", max_depth);
    if (node->is_leaf()) {
      if (node_depth < *min_depth) *min_depth = node_depth;
      if (node_depth > *max_depth) *max_depth = node_depth;
    } else {
      if (node->lo_child()) tree_depth_update_in(node->lo_child(), node_depth+1, min_depth, max_depth);
      if (node->hi_child()) tree_depth_update_in(node->hi_child(), node_depth+1, min_depth, max_depth);
    }
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  void kd_tree<G_dimension, G_value, G_data>::tree_depth_in(std::size_t* min_depth,
                                                            std::size_t* max_depth) const {
    SL_REQUIRE("Min depth exists", min_depth);
    SL_REQUIRE("Max depth exists", max_depth);
    if (root_) {
      *min_depth = count()+1;
      *max_depth = 0;
      tree_depth_update_in(root_, 1, min_depth, max_depth);
    } else {
      *min_depth = 0;
      *max_depth = 0;
    }      
  }  
  
  template <std::size_t G_dimension, class G_value,class G_data>
  bool kd_tree<G_dimension, G_value, G_data>::is_lo(const node_t* node,
                                                    std::size_t i,
                                                    value_t x_i) const {
    bool result = true;
    if (node) {
      result = result && node->data()[i] <= x_i;
      result = result && is_lo(node->lo_child(), i, x_i);
      result = result && is_lo(node->hi_child(), i, x_i);
    }
    return result;
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  bool kd_tree<G_dimension, G_value, G_data>::is_hi(const node_t* node,
                                                    std::size_t i,
                                                    value_t x_i) const {
    bool result = true;
    if (node) {
      result = result && node->data()[i] > x_i;
      result = result && is_hi(node->lo_child(), i, x_i);
      result = result && is_hi(node->hi_child(), i, x_i);
    }
    return result;
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  bool kd_tree<G_dimension, G_value, G_data>::invariant() const {
    bool result = true;

    const std::size_t tree_sz = (root_ == 0) ? 0 : root_->tree_size();

    if (result) {
      result = (count() == tree_sz);

      if (!result) SL_TRACE_OUT(0) << "kd_tree::invariant() : count() != tree_sz !" << std::endl;
    }

    if (result) {
      result = node_invariant(root_);
    }

    return result;
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  bool kd_tree<G_dimension, G_value, G_data>::node_invariant(const node_t* node) const {
    bool result = true;
    if (node) {
      const node_t* l = node->lo_child();
      const node_t* h = node->hi_child();
      
      if (result) {
        const std::size_t l_sz = (l==0) ? 0 : l->tree_size();
        const std::size_t h_sz = (h==0) ? 0 : h->tree_size();
      
        result = (node->tree_size() == (1 + l_sz + h_sz));

        if (!result) SL_TRACE_OUT(0) << "kd_tree::invariant(node) : inconsistent tree size !" << std::endl;
      }

      if (result) {
        result = is_lo(node->lo_child(), node->discriminant(), node->discriminant_plane_constant());

        if (!result) SL_TRACE_OUT(0) << "kd_tree::invariant(node) : lo split discriminant test failed !" << std::endl;
      }

      if (result) {
        result = is_hi(node->hi_child(), node->discriminant(), node->discriminant_plane_constant());

        if (!result) SL_TRACE_OUT(0) << "kd_tree::invariant(node) : hi split discriminant test failed !" << std::endl;
      }
	
      if (result) {
        // Check children invariant
        result = result && node_invariant(l);
        result = result && node_invariant(h);
      }
    }
    return result;
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  std::size_t kd_tree<G_dimension, G_value, G_data>::external_path_length() const {
    return external_path_length(root_, 1);
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  std::size_t kd_tree<G_dimension, G_value, G_data>::internal_path_length() const {
    return internal_path_length(root_, 1);
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  std::size_t kd_tree<G_dimension, G_value, G_data>::external_path_length(const node_t* node,
                                                                          std::size_t depth) const {
    std::size_t result = 0;
    if (node) {
      if (node->is_leaf()) {
        result = depth;
      } else {
        result = 
          external_path_length(node->lo_child(), depth+1) +
          external_path_length(node->hi_child(), depth+1);
      }	  
    }
    return result;
  }


  template <std::size_t G_dimension, class G_value,class G_data>
  std::size_t kd_tree<G_dimension, G_value, G_data>::internal_path_length(const node_t* node,
                                                                          std::size_t depth) const {
    std::size_t result = 0;
    if (node) {
      if (node->is_leaf()) {
        result = 0;
      } else {
        result = depth +
          internal_path_length(node->lo_child(), depth+1) +
          internal_path_length(node->hi_child(), depth+1);
      }	  
    }
    return result;
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  void kd_tree<G_dimension, G_value, G_data>::split(node_t* t,
                                                    std::size_t i,
                                                    value_t x_i,
                                                    node_t** t_lo,
                                                    node_t** t_hi) {
    SL_REQUIRE("Good i", i<dimension);
    SL_REQUIRE("t_lo exists", t_lo);
    SL_REQUIRE("t_hi exists", t_hi);
    
    if (t == 0) {
      *t_lo = 0;
      *t_hi = 0;
    } else {
      std::size_t j = t->discriminant();
      value_t y_i = t->data()[i];
      if (i==j) {
        if (y_i <= x_i) {
          node_t *t_hi_lo;
          node_t *t_hi_hi;
          split(t->hi_child(), i, x_i, &t_hi_lo, &t_hi_hi);
          t->set_hi_child(t_hi_lo);
          *t_lo = t;
          *t_hi = t_hi_hi;
        } else {
          SL_CHECK("Ordered", y_i > x_i);
          node_t *t_lo_lo;
          node_t *t_lo_hi;
          split(t->lo_child(), i, x_i, &t_lo_lo, &t_lo_hi);
          t->set_lo_child(t_lo_hi);
          *t_lo = t_lo_lo;
          *t_hi = t;
        }
      } else {
        SL_CHECK("Different", i!=j);
        node_t *t_lo_lo;
        node_t *t_lo_hi;
        split(t->lo_child(), i, x_i, &t_lo_lo, &t_lo_hi);
        node_t *t_hi_lo;
        node_t *t_hi_hi;
        split(t->hi_child(), i, x_i, &t_hi_lo, &t_hi_hi);
        if (y_i <= x_i) {
          t->set_lo_child(t_lo_lo);
          t->set_hi_child(t_hi_lo);
          *t_lo = t;
          *t_hi = join(t_lo_hi,t_hi_hi,j);
        } else {
          *t_lo = join(t_lo_lo,t_hi_lo,j);
          t->set_lo_child(t_lo_hi);
          t->set_hi_child(t_hi_hi);
          *t_hi = t;
        }
      }

      // Update tree size
      t->set_tree_size(1);
      if (t->lo_child()) t->set_tree_size(t->tree_size() + t->lo_child()->tree_size());
      if (t->hi_child()) t->set_tree_size(t->tree_size() + t->hi_child()->tree_size());
    }
  }

  template <std::size_t G_dimension, class G_value,class G_data>
  typename kd_tree<G_dimension, G_value, G_data>::node_t *kd_tree<G_dimension, G_value, G_data>::join(typename kd_tree<G_dimension, G_value, G_data>::node_t * a,
                                                                                                      typename kd_tree<G_dimension, G_value, G_data>::node_t * b,
                                                                                                      std::size_t i) {
    SL_REQUIRE("Good dimension", i<dimension);
    
    if (a == 0) {
      return b; 
    } else if (b == 0) {
      return a;
    } else {
      std::size_t n = a->tree_size();
      std::size_t m = b->tree_size();

      if (n>m) {
        // select a as root of the join
        std::size_t j_a = a->discriminant();
        if (i==j_a) {
          a->set_hi_child(join(a->hi_child(), b, i));
        } else {
          node_t* b_lo;
          node_t* b_hi;
          split(b, j_a, a->discriminant_plane_constant(), &b_lo, &b_hi);
          a->set_lo_child(join(a->lo_child(),b_lo, i));
          a->set_hi_child(join(a->hi_child(),b_hi, i));
        }
        a->set_tree_size(n+m);
        return a;
      } else {
        // select b as root of the join
        std::size_t j_b = b->discriminant();
        if (i==j_b) {
          b->set_lo_child(join(a, b->lo_child(), i));
        } else {
          node_t* a_lo;
          node_t* a_hi;
          split(a, j_b, b->discriminant_plane_constant(), &a_lo, &a_hi);
          b->set_lo_child(join(a_lo, b->lo_child(), i));
          b->set_hi_child(join(a_hi, b->hi_child(), i));
        }
        b->set_tree_size(n+m);
        return b;
      }  
    }
  }


}

#endif
