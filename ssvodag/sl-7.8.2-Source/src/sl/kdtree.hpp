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
#ifndef SL_KDTREE_HPP
#define SL_KDTREE_HPP

#include <sl/fixed_size_point.hpp>
#include <sl/utility.hpp>
#include <sl/random.hpp>
#include <sl/math.hpp>
#include <sl/numeric_traits.hpp>
#include <sl/reservable_priority_queue.hpp>
#include <sl/arena_allocator.hpp>
#include <cassert>
#include <utility>
#include <vector>

#include <sl/fsb_allocator.hpp>

namespace sl {
  
  template <std::size_t G_dimension, class G_subvalue, class G_value>
  class kdtree;
    
  /**
   *  A simple object storing K coordinates
   */
  template <std::size_t G_dimension, class G_subvalue>
  class kdcoords {
  public:
    enum { dimension = G_dimension };
    
    typedef kdcoords<G_dimension, G_subvalue> this_t;
    typedef G_subvalue subvalue_t;
  protected:
    subvalue_t c_[dimension];
  public:
    inline kdcoords() {
      for (std::size_t i=0; i<dimension; ++i) c_[i] = 0;
    }
    
    inline subvalue_t operator[](std::size_t i) const {
      assert(i<dimension);
      return c_[i];
    }
    inline subvalue_t &operator[](std::size_t i) {
      assert(i<dimension);
      return c_[i];
    }
  };
       
  /**
   *  A simple predicate always returning true
   */
  class kdalways {
  public:
    inline kdalways() {}

    template <class T>
    inline bool operator()(const T&) const { return true; }
  };

  /**
   *  A simple predicate returning true for values equal to
   *  a reference value
   */
  template <class T>
  class kdsame {
  protected:
    T x_;
  public:
    inline kdsame(const T& x): x_(x) {}
    inline bool operator()(const T& x) const { return x == x_; }
  };

  /**
   *  A simple predicate returning true for values whose
   *  coords are equal to a reference value
   */
  template <std::size_t N, class T>
  class kdsamecoords {
  protected:
    T x_;
  public:
    inline kdsamecoords(const T& x): x_(x) {}
    inline bool operator()(const T& x) const { 
      bool result = true;
      for (std::size_t i=0; i<N && result; ++i) {
	result = result && (x[i]==x_[i]);
      }
      return result;
    }
  };
     
  /**
   *  A node in a kdtree
   */
  template <std::size_t G_dimension, class G_subvalue, class G_value>
  class kdnode {
  public:
    enum { dimension = G_dimension };
    
    typedef kdnode<G_dimension, G_subvalue, G_value> this_t;
    typedef G_subvalue subvalue_t;
    typedef G_value  value_t;
  protected:
    kdnode* parent_;
    kdnode* lo_child_;
    kdnode* hi_child_;
    std::size_t discriminant_;
    std::size_t tree_size_;
    value_t value_;
  public:
    
    /// Init an empty kdnode
    inline kdnode()
      :
      parent_(0), lo_child_(0), hi_child_(0), discriminant_(0), tree_size_(1), value_(value_t()) {
    }
    
    /// Init a kdnode containing x
    inline kdnode(const value_t& x) 
      :
      parent_(0), lo_child_(0), hi_child_(0), discriminant_(0), tree_size_(1), value_(x) {
    }
    
    /// The number of nodes in the tree rooted at this (including this)
    inline std::size_t tree_size() const {
      return tree_size_;
    }
    
    /// Set to x the number of nodes in the tree rooted at this (including this)
    inline void set_tree_size(std::size_t x) {
      assert(x >= 1);
      tree_size_ = x;
    }
    
    /// The value associated to this node
    inline const value_t& value() const {
      return value_;
    }
    
    /// The value associated to this node
    inline value_t& value() {
      return value_;
    }
    
    /// Is this a leaf node?
    inline bool is_leaf() const {
      return !lo_child_ && !hi_child_;
    }
    
    /// Is this the root node?
    inline bool is_root() const {
      return parent_ == 0;
    }

    // The parent of this node
    inline const kdnode* parent() const {
      return parent_;
    }

    // The parent of this node
    inline kdnode* parent() {
      return parent_;
    }

    /// The child on the negative side of the discriminating discriminant
    inline const kdnode* lo_child() const {
      return lo_child_;
    }

    /// The child on the negative side of the discriminating discriminant
    inline kdnode* lo_child() {
      return lo_child_;
    }

    /// Set to x the child on the negative side of the discriminating discriminant
    inline void set_lo_child(kdnode* x) {
      if (lo_child_ && lo_child_->parent_ == this) lo_child_->parent_ = 0;
      lo_child_ = x;
      if (lo_child_) lo_child_->parent_ = this;
    }

    /// The child on the positive side of the discriminating discriminant
    inline const kdnode* hi_child() const {
      return hi_child_;
    }
    
    /// The child on the positive side of the discriminating discriminant
    inline kdnode* hi_child() {
      return hi_child_;
    }

    /// Set to x the child on the positive side of the discriminating discriminant
    inline void set_hi_child(kdnode* x) {
      if (hi_child_ && hi_child_->parent_ == this) hi_child_->parent_ = 0;
      hi_child_ = x;
      if (hi_child_) hi_child_->parent_ = this;
    }

    /// The split discriminant associated to this node
    inline std::size_t discriminant() const {
      return discriminant_;
    }

    // The constant of the discriminant
    inline subvalue_t discriminant_plane_constant() const {
      return value_[discriminant()];
    }
      
    /// Set the spltting discriminant associated to this node to x
    inline void set_discriminant(std::size_t x) {
      assert(x<dimension);
      assert(is_leaf());
      discriminant_ = x;
    }

    /// Apply action a to the tree rooted at this
    template <class action_t>
    inline void apply(action_t& a) const {
      a.apply(value());
      if (lo_child()) lo_child()->apply(a);
      if (hi_child()) hi_child()->apply(a);
    }
    
  };

  template <class G_node, class G_value>
  class kditerator {
  public:
    typedef kditerator<G_node, G_value> this_t;
    typedef G_node                      node_t;
    typedef G_value                     value_t;

  public: // std::iterator types

    typedef value_t        value_type;
    typedef value_t&       reference;
    typedef const value_t& const_reference; 
    typedef value_t*       pointer;
    typedef const value_t* const_pointer;

    typedef kditerator<node_t, value_t> iterator;
    typedef kditerator<const node_t, const value_t> const_iterator;

    // Iterator could be bidirectional, but it's not implemented yet
    // because we currently use a null pointer as end() marker
    typedef std::forward_iterator_tag iterator_category; 
    typedef ptrdiff_t difference_type;

  protected:
    node_t* node_;

  public:

    kditerator(node_t* n = 0): node_(n) {
    }
     
    kditerator(const iterator& it): node_(it.node_) {
    }

    operator const_iterator() {
      return const_iterator(node_);
    }

    bool operator == (const iterator& it) {
      return node_ == it.node_;
    }
    bool operator != (const iterator& it) {
      return node_ != it.node_;
    }

    reference operator*() const {
      assert(node_);
      return node_->value();
    }

    pointer operator->() const {
      assert(node_);
      return &(operator*());
    }
    
    this_t& operator++() {
      increment(); 
      return *this;
    }
    
    this_t operator++(int) {
      this_t old_this = *this;
      increment();
      return old_this;
    }

  public:

    node_t* node_ptr() const {
      return node_;
    }

  protected:
    inline void increment() {
      assert(node_);

      // Preorder iteration, starts with root
      if (node_->lo_child()) {
	node_ = node_->lo_child();
      } else if (node_->hi_child()) {
	node_ = node_->hi_child();
      } else {
	// We've just visited a leaf node.
	// Go back up the tree until we find a node
	// with a hi child that we haven't seen yet.
	node_t* parent= node_->parent();
	node_t* child = node_;
	while (parent != 0 &&
	       ((parent->hi_child() == child) || 
		(parent->hi_child() == 0))) {
	  child = parent;
	  parent = parent->parent();
	}
	if (parent) {
	  node_ = parent->hi_child();
	} else {
	  node_ = 0;
	}
      }
    }

  };
    
  namespace detail {
    
    template <std::size_t G_dimension, class G_subvalue, class G_value>
    struct kdt_priority_search_queue_element {
      typedef G_subvalue subvalue_t;
      typedef G_value  value_t;
      typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances

      typedef kdcoords<G_dimension, G_subvalue> kdcoords_t;
      typedef kdnode<G_dimension, G_subvalue, G_value> node_t;

      node_t*    node_;
      kdcoords_t box_offset_;
      float_t    box_dist2_;

      inline kdt_priority_search_queue_element(node_t* node,
					       const kdcoords_t& box_offset,
					       float_t box_dist2) :
	node_(node), box_offset_(box_offset), box_dist2_(box_dist2) {
      }

      inline bool operator<(const kdt_priority_search_queue_element& other) const {
	// Priority is inverse of distance -- so that small distance nodes are visited first
	return box_dist2_ > other.box_dist2_;
      }
    };

    /**
     * Squared Euclidean distance functor, optimized version
     * This is highly optimised, with loop unrolling, as it is one
     * of the most expensive inner loops.
     */
    template <std::size_t G_dimension, class G_subvalue, class G_value>
    struct kdt_L2 {
      typedef G_subvalue subvalue_t;
      typedef G_value  value_t;
      typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances
            
      template <class G_querypoint>
      static float_t value(const value_t& node_value, const G_querypoint& p, const float_t knn_maxdist2) {
	float_t result = float_t(0);
	
	std::size_t d =0;
	// Process 4 items with each loop for efficiency
	while (d+4<G_dimension) {
	  const float_t diff0 = float_t(node_value[d+0] - p[d+0]);
	  const float_t diff1 = float_t(node_value[d+1] - p[d+1]);
	  const float_t diff2 = float_t(node_value[d+2] - p[d+2]);
	  const float_t diff3 = float_t(node_value[d+3] - p[d+3]);
	  
	  result += diff0*diff0 + diff1*diff1 + diff2*diff2 + diff3*diff3;
	
	  d+= 4;
	  if (result>knn_maxdist2) return result;
	}
	// Process last 0-3 items
	while (d < G_dimension) {
	  const float_t diff0 = float_t(node_value[d] - p[d]);
	  result += diff0 * diff0;
	  ++d;
	}
	return result;
      }
    };
    
    /**
     * Squared Euclidean distance functor, optimized version for dimension=1
     */
    template <class G_subvalue, class G_value>
    struct kdt_L2<1, G_subvalue, G_value> {
      typedef G_subvalue subvalue_t;
      typedef G_value  value_t;
      typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances
            
      template <class G_querypoint>
      static inline float_t value(const value_t& node_value, const G_querypoint& p, const float_t knn_maxdist2) {
	SL_USEVAR(knn_maxdist2);
	const float_t diff0 = float_t(node_value[0] - p[0]);
	return diff0*diff0;
      }
    };
    
    /**
     * Squared Euclidean distance functor, optimized version for dimension=2
     */
    template <class G_subvalue, class G_value>
    struct kdt_L2<2, G_subvalue, G_value> {
      typedef G_subvalue subvalue_t;
      typedef G_value  value_t;
      typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances
      
      template <class G_querypoint>
      static inline float_t value(const value_t& node_value, const G_querypoint& p, const float_t knn_maxdist2) {
	SL_USEVAR(knn_maxdist2);
	const float_t diff0 = float_t(node_value[0] - p[0]);
	const float_t diff1 = float_t(node_value[1] - p[1]);
	return diff0*diff0 + diff1*diff1;
      }
    };

    /**
     * Squared Euclidean distance functor, optimized version for dimension=3
     */
    template <class G_subvalue, class G_value>
    struct kdt_L2<3, G_subvalue, G_value> {
      typedef G_subvalue subvalue_t;
      typedef G_value  value_t;
      typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances
      
      template <class G_querypoint>
      static inline float_t value(const value_t& node_value, const G_querypoint& p, const float_t knn_maxdist2) {
	SL_USEVAR(knn_maxdist2);
	const float_t diff0 = float_t(node_value[0] - p[0]); 
	const float_t diff1 = float_t(node_value[1] - p[1]); 
	const float_t diff2 = float_t(node_value[2] - p[2]);
	
	return diff0*diff0 + diff1*diff1 + diff2*diff2;
      }
    };

    /**
     * Squared Euclidean distance functor, optimized version for dimension=4
     */
    template <class G_subvalue, class G_value>
    struct kdt_L2<4, G_subvalue, G_value> {
      typedef G_subvalue subvalue_t;
      typedef G_value  value_t;
      typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances
      
      template <class G_querypoint>
      static inline float_t value(const value_t& node_value, const G_querypoint& p, const float_t knn_maxdist2) {
	SL_USEVAR(knn_maxdist2);
	const float_t diff0 = float_t(node_value[0] - p[0]);
	const float_t diff1 = float_t(node_value[1] - p[1]);
	const float_t diff2 = float_t(node_value[2] - p[2]);
	const float_t diff3 = float_t(node_value[3] - p[3]);
	return diff0*diff0 + diff1*diff1 + diff2*diff2 + diff3*diff3;
      }
    };

  } // namespace detail
  
  /**
   *  Randomized kd-trees of dimension G_dimension containing
   *  entities of type G_value with coordinates of
   *  type G_subvalue. The class G_value
   *  is constrained to export G_subvalue operator[](std::size_t i)
   *  to provide read access to coordinates 0 to dimension.
   *  The insertion procedure uses a randomized approach in
   *  order to keep the tree (probabilistically) balanced
   *  independently from the value.
   *
   *  Randomized KD tree based on:
   *
   *  @inproceedings{ duch98randomized,
   *     author = {Duch and Estivill-Castro and Martinez},
   *     title = {Randomized {K}-Dimensional Binary Search Trees},
   *     booktitle = {{ISAAC}: 9th International Symposium on Algorithms and Computation 
   *              (formerly {SIGAL} International Symposium on Algorithms), 
   *		 Organized by Special Interest Group on Algorithms ({SIGAL}) 
   *		 of the Information Processing Society of Japan ({IPSJ}) 
   *		 and the Technical Group on Theoretical Foundation of Computing  
   *		 of the Institute of Electronics, Information and Communication Engineers ({IEICE})},
   *	 year = {1998}
   *  }
   *
   *  KNN search with incremental distance computation based on:
   *  S. Arya and D. M. Mount, Algorithms for fast vector quantization,
   *  Proc. IEEE Data Compression Conference (DCC), eds. J. A.
   *  Storer and M. Cohn, IEEE Press, 381-390, 1993.
   *
   *  KNN search with Bounds-Overlap-Threshold (BOT) test based on:
   *  Ryusuke Sagawa, Tomohito Masuda, Katsushi Ikeuchi, "Effective Nearest
   *  Neighbor Search for Aligning and Merging Range Images," 
   *  Fourth International Conference on 3-D Digital Imaging and Modeling
   *  (3DIM '03), 79-87, 2003
   */
  template <std::size_t G_dimension, class G_subvalue, class G_value>
  class kdtree {
  public:
    enum { dimension = G_dimension };

    typedef kdtree<G_dimension, G_subvalue, G_value> this_t;
    typedef G_value  value_t;

    typedef G_subvalue subvalue_t;
    typedef typename sl::numeric_traits<subvalue_t>::T_floattype float_t; // For distances

    typedef kdcoords<G_dimension, G_subvalue> kdcoords_t;
    typedef kdnode<G_dimension, G_subvalue, G_value> node_t;
    typedef detail::kdt_L2<G_dimension, G_subvalue, G_value> distance_t;
    
    typedef value_t        value_type;
    typedef value_t&       reference;
    typedef const value_t& const_reference;
    typedef value_t*       pointer;
    typedef const value_t* const_pointer;

    typedef kditerator<node_t, value_t> iterator;
    typedef kditerator<const node_t, const value_t> const_iterator;

  protected:
 
    // Allocator
    fsb_allocator<node_t> node_allocator_;

    // Root of the tree
    node_t* root_;

    // Number of points
    std::size_t count_;
 
    // Maximum distance constant
    float_t maxd2_;
    float_t sqrt_maxd2_;

    // Search mode
    bool    is_priority_search_enabled_;
    std::size_t priority_search_max_visited_count_;

    mutable std::size_t stat_knn_query_count_;
    mutable std::size_t stat_knn_visited_count_;
    
  protected: // Random number generation
     
    random::std_irng_t rng_;

    /// Reset the pseudo random number generator
    inline void random_reset() {
      rng_.set_seed(); // default
    }

    /// An unsigned random number between l and h 
    inline unsigned int random_unsigned(unsigned int l, unsigned int h) {
      assert(h >= l);
      return rng_.value_within(l,h);
    }
      
    /// A random discriminant
    inline std::size_t random_discriminant() {
      return (std::size_t)(random_unsigned(0,dimension-1));
    }
      
    /// A random test returning true with probability a/b
    inline bool random_test(unsigned int a, unsigned int b) {
      assert(b>0);
      assert(b>=a);
      return random_unsigned(1,b) <= a;
    } 

  public:

    /// Create an empty kdtree
    inline kdtree() : 
      root_(0), count_(0), 
      maxd2_(sl::scalar_math<float_t>::upper_bound()),
      sqrt_maxd2_(std::sqrt(sl::scalar_math<float_t>::finite_upper_bound())),
      is_priority_search_enabled_(true),
      priority_search_max_visited_count_(std::size_t(-1)) {
      random_reset();
    }
      
    /// Create a kdtree as a deep copy of other
    kdtree(const this_t& other): 
      root_(0), count_(0), 
      maxd2_(sl::scalar_math<float_t>::upper_bound()),
      sqrt_maxd2_(std::sqrt(sl::scalar_math<float_t>::finite_upper_bound())),
      is_priority_search_enabled_(other.is_priority_search_enabled()),
      priority_search_max_visited_count_(other.is_priority_search_max_visited_count()) {
      random_reset();
      insert(other.begin(),other.end());
    }
    
    /// Create kdtree and insert all given elements
    template <class G_input_iterator>
    kdtree(G_input_iterator in_begin, G_input_iterator in_end) : 
      root_(0), count_(0), 
      maxd2_(sl::scalar_math<float_t>::upper_bound()),
      sqrt_maxd2_(std::sqrt(sl::scalar_math<float_t>::finite_upper_bound())),
      is_priority_search_enabled_(true),
      priority_search_max_visited_count_(std::size_t(-1)) {
      random_reset();
      insert(in_begin, in_end);
    }
      
    /// Destroy a kdtree
    inline ~kdtree() {
      clear();
    }

    /// Assign other to this kdtree (deep copy)
    this_t& operator=(const this_t& other) {
      if (this != &other) {
	clear();
	insert(other.begin(),other.end());
      }
      return *this;
    }

  public: // Simple queries

    /// The number of value points
    inline std::size_t count() const {
      return count_;
    }

    /// Is the tree empty
    inline bool is_empty() const {
      return count() == 0;
    }

  public: // KNN stats

    inline void stat_reset() {
      stat_knn_visited_count_ = 0;
      stat_knn_query_count_ = 0;
    }

    inline std::size_t stat_knn_query_count() const {
      return stat_knn_query_count_;
    }

    inline std::size_t stat_knn_visited_count() const {
      return stat_knn_visited_count_;
    }
    
  public: // KNN

    // Use priority search (x=true) or depth first search (x=false) for KNN queries
    void set_is_priority_search_enabled(bool x) {
      is_priority_search_enabled_ = x;
    }
    
    // Use priority search rather than depth first search for KNN queries
    bool is_priority_search_enabled() const {
      return is_priority_search_enabled_;
    }
 
    // Set max number of nodes visited during priority search
    void set_priority_search_max_visited_count(std::size_t x) {
      priority_search_max_visited_count_ = x;
    }
    
    // Max number of nodes visited during priority search
    std::size_t priority_search_max_visited_count() const {
      return priority_search_max_visited_count_;
    }
   
    /**
     *  Find a maximum of knn_maxsize points meeting the condition 
     *  that are closest to the query point p and are in the sphere 
     *  centered at p with radius knn_maxradius. The query is
     *  answered with a precision eps. 
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint, typename G_predicate>
    void k_approximately_nearest_neighbors_in(std::size_t&    knn_size,
					      const_iterator *knn_it,
					      float_t        *knn_dist2,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize,
					      float_t         knn_maxradius,
					      const G_predicate& condition) const {
      knn_size = 0;
      if (root_ && knn_maxsize>0) {
	float_t knn_maxdist2 = (knn_maxradius>=sqrt_maxd2_)?(maxd2_):(knn_maxradius*knn_maxradius);
	float_t sqr_one_plus_eps = sl::sqr(float_t(1.0)+eps);
	float_t sqr_bot = (bot>=sqrt_maxd2_)?(maxd2_):(bot*bot);
	
	internal_k_approximately_nearest_neighbors_in(knn_size,
						      knn_it,
						      knn_dist2,
						      p, 
						      sqr_one_plus_eps,
						      sqr_bot,
						      knn_maxsize,
						      knn_maxdist2,
						      condition);
      }
    }
    
    /**
     *  Find a maximum of knn_maxsize points meeting the condition 
     *  that are closest to the query point p and are in the sphere 
     *  centered at p with radius knn_maxradius. 
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     */
    template <typename G_querypoint, typename G_predicate>
    void k_nearest_neighbors_in(std::size_t&    knn_size,
				const_iterator *knn_it,
				float_t        *knn_dist2,
				const G_querypoint& p,
				std::size_t     knn_maxsize,
				float_t         knn_maxradius,
				const G_predicate& condition) const {
      const float_t eps = float_t(0.0);
      const float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, condition);
    }

    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint>
    void k_approximately_nearest_neighbors_in(std::size_t&    knn_size,
					      const_iterator *knn_it,
					      float_t     *knn_dist2,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize,
					      float_t         knn_maxradius) const {
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, kdalways());
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     */
    template <typename G_querypoint>
    void k_nearest_neighbors_in(std::size_t&    knn_size,
				const_iterator *knn_it,
				float_t        *knn_dist2,
				const G_querypoint& p,
				std::size_t     knn_maxsize,
				float_t         knn_maxradius) const {
      const float_t eps = float_t(0.0);
      const float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, kdalways());
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint>
    void k_approximately_nearest_neighbors_in(std::size_t&    knn_size,
					      const_iterator *knn_it,
					      float_t     *knn_dist2,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize) const {
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, maxd2_, kdalways());
    }

    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     */
    template <typename G_querypoint>
    void k_nearest_neighbors_in(std::size_t&    knn_size,
				const_iterator *knn_it,
				float_t     *knn_dist2,
				const G_querypoint& p, 
				std::size_t     knn_maxsize) const {
      const float_t eps = float_t(0.0);
      const float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, maxd2_, kdalways());
    }

    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint, typename G_predicate>
    void k_approximately_nearest_neighbors_in(std::vector<value_t>& knn,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize,
					      float_t      knn_maxradius,
					      const G_predicate& condition) const {
      std::size_t knn_size;
      const_iterator* knn_it = new const_iterator[knn_maxsize];
      float_t*     knn_dist2 = new float_t[knn_maxsize];
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, condition);
      knn.resize(knn_size);
      for (std::size_t k=0; k<knn_size; ++k) {
	knn[k] = *(knn_it[k]);
      }
      delete[] knn_it;
      delete[] knn_dist2;
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     */
    template <typename G_querypoint, typename G_predicate>
    void k_nearest_neighbors_in(std::vector<value_t>& knn,
				const G_querypoint& p, 
				std::size_t     knn_maxsize,
				float_t      knn_maxradius,
				const G_predicate& condition) const {
      const float_t eps = float_t(0.0);
      const float_t bot = sqrt_maxd2_;
      std::size_t knn_size;
      const_iterator* knn_it = new const_iterator[knn_maxsize];
      float_t*     knn_dist2 = new float_t[knn_maxsize];
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, condition);
      knn.resize(knn_size);
      for (std::size_t k=0; k<knn_size; ++k) {
	knn[k] = *(knn_it[k]);
      }
      delete[] knn_it;
      delete[] knn_dist2;
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint>
    void k_approximately_nearest_neighbors_in(std::vector<value_t>& knn,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize,
					      float_t      knn_maxradius) const {
      std::size_t knn_size;
      const_iterator* knn_it = new const_iterator[knn_maxsize];
      float_t*     knn_dist2 = new float_t[knn_maxsize];
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, kdalways());
      knn.resize(knn_size);
      for (std::size_t k=0; k<knn_size; ++k) {
	knn[k] = *(knn_it[k]);
      }
      delete[] knn_it;
      delete[] knn_dist2;
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     */
    template <typename G_querypoint>
    void k_nearest_neighbors_in(std::vector<value_t>& knn,
				const G_querypoint& p, 
				std::size_t     knn_maxsize,
				float_t      knn_maxradius) const {
      float_t eps = float_t(0.0);
      float_t bot = sqrt_maxd2_;
      std::size_t knn_size;
      const_iterator* knn_it = new const_iterator[knn_maxsize];
      float_t*     knn_dist2 = new float_t[knn_maxsize];
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, kdalways());
      knn.resize(knn_size);
      for (std::size_t k=0; k<knn_size; ++k) {
	knn[k] = *(knn_it[k]);
      }
      delete[] knn_it;
      delete[] knn_dist2;
    }

    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly.
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint>
    void k_approximately_nearest_neighbors_in(std::vector<value_t>& knn,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize) const {
      k_approximately_nearest_neighbors_in(knn, p, eps, bot, knn_maxsize, maxd2_);
    }

    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p.
     */
    template <typename G_querypoint>
    void k_nearest_neighbors_in(std::vector<value_t>& knn,
				const G_querypoint& p, 
				std::size_t     knn_maxsize) const {
      const float_t eps = float_t(0.0);
      const float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn, p, eps, bot, knn_maxsize, maxd2_);
    }
 
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint, typename G_predicate>
    void k_approximately_nearest_neighbors_in(std::size_t&    knn_size,
					      iterator       *knn_it,
					      float_t     *knn_dist2,
					      const G_querypoint& p,
					      float_t         eps,
					      float_t         bot,
					      std::size_t     knn_maxsize,
					      float_t      knn_maxradius,
					      const G_predicate& condition) {
      knn_size = 0;
      if (root_ && knn_maxsize>0) {
	float_t knn_maxdist2 = (knn_maxradius>=sqrt_maxd2_)?(maxd2_):(knn_maxradius*knn_maxradius);
	float_t sqr_one_plus_eps = sl::sqr(float_t(1.0)+eps);
	float_t sqr_bot = (bot>=sqrt_maxd2_)?(maxd2_):(bot*bot);

	internal_k_approximately_nearest_neighbors_in(knn_size,
						      knn_it,
						      knn_dist2,
						      p, 
						      sqr_one_plus_eps,
						      sqr_bot,
						      knn_maxsize,
						      knn_maxdist2,
						      condition);
      }
    }
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     */
    template <typename G_querypoint, typename G_predicate>
    void k_nearest_neighbors_in(std::size_t&    knn_size,
				iterator       *knn_it,
				float_t     *knn_dist2,
				const G_querypoint& p, 
				std::size_t     knn_maxsize,
				float_t      knn_maxradius,
				const G_predicate& condition) {
      const float_t eps = float_t(0.0);
      const float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2,
					   p,
					   eps,
					   bot,
					   knn_maxsize, knn_maxradius,
					   condition);
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly.
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint>
    void k_approximately_nearest_neighbors_in(std::size_t&    knn_size,
					      iterator *knn_it,
					      float_t     *knn_dist2,
					      const G_querypoint& p,
					      float_t             eps,
					      float_t             bot,
					      std::size_t     knn_maxsize,
					      float_t      knn_maxradius) {
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, kdalways());
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     */
    template <typename G_querypoint>
    void k_nearest_neighbors_in(std::size_t&    knn_size,
				iterator *knn_it,
				float_t     *knn_dist2,
				const G_querypoint& p, 
				std::size_t     knn_maxsize,
				float_t      knn_maxradius) {
      float_t eps = float_t(0.0);
      float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, knn_maxradius, kdalways());
    }

    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     *  A real value eps >=0 may be supplied. If so, then the ith nearest
     *  neighbor is a (1+eps) approximation to the true ith nearest neighbor.
     *  That is, the true (not squared) distance to this point may exceed
     *  the true distance to the real ith nearest neighbor of q by a factor
     *  of (1 + eps). If eps=0 then nearest neighbors are computed exactly
     *  if bot=oo.
     *  If the nearest neighbor point is far from a query, many of the
     *  leaf nodes must be examined during the search, which actually will not
     *  finish in logarithmic time. In order to speed-up search in those cases,
     *  during search, branches are pruned if the possible nearest neighbor
     *  point lies beyond a certain threshold bot. This may introduce errors
     *  for far queries.
     */
    template <typename G_querypoint>
    void k_approximately_nearest_neighbors_in(std::size_t&    knn_size,
					      iterator *knn_it,
					      float_t     *knn_dist2,
					      const G_querypoint& p,
					      float_t             eps,
					      float_t         bot,
					      std::size_t     knn_maxsize) {
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, maxd2_, kdalways());
    }
    
    /**
     *  Find a maximum of knn_maxsize points that are closest to the 
     *  query point p and are in the sphere centered at p with radius knn_maxradius.
     *  On exit, knn_size contains the number of neighbors found, 
     *  knn_it[] the iterator to the neighbors, and knn_dist the 
     *  square of the Euclidean distance to the
     *  query point. Neighbors are sorted by increasing distance.
     *  The knn_it, knn_dist2 arrays are required to contain at least
     *  knn_maxsize elements.
     */
    template <typename G_querypoint>
    void k_nearest_neighbors_in(std::size_t&    knn_size,
				iterator *knn_it,
				float_t     *knn_dist2,
				const G_querypoint& p, 
				std::size_t     knn_maxsize) {
      float_t eps = float_t(0.0);
      float_t bot = sqrt_maxd2_;
      k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2, p, eps, bot, knn_maxsize, maxd2_, kdalways());
    }

  protected:

    template <typename G_querypoint, typename G_predicate, typename G_iter>
    void internal_k_approximately_nearest_neighbors_in(std::size_t&    knn_size,    // OUT
						       G_iter         *knn_it,      // OUT
						       float_t        *knn_dist2,   // OUT
						       const G_querypoint& p,       // IN
						       float_t sqr_one_plus_eps,    // IN
						       float_t sqr_bot,             // IN
						       std::size_t knn_maxsize,     // IN
						       float_t& knn_maxdist2,       // IN/OUT
						       const G_predicate& condition // IN
						       ) const {

      ++stat_knn_query_count_;
      
      if (is_priority_search_enabled_) {
	internal_prioritized_k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2,
								  p, sqr_one_plus_eps, sqr_bot,
								  knn_maxsize,
								  knn_maxdist2,
								  condition);
      } else {
	float_t box_dist2 = float_t(0.0);
	kdcoords_t box_offset; // init to 0
	internal_depth_first_k_approximately_nearest_neighbors_in(knn_size, knn_it, knn_dist2,
								  root_,
								  p, sqr_one_plus_eps, sqr_bot,
								  knn_maxsize,
								  knn_maxdist2,
								  box_dist2,
								  box_offset,
								  condition);
      }
    }
    
    template <typename G_querypoint, typename G_predicate, typename G_iter>
    void internal_depth_first_k_approximately_nearest_neighbors_in(std::size_t&    knn_size,    // OUT
								   G_iter         *knn_it,      // OUT
								   float_t        *knn_dist2,   // OUT
								   node_t* node,                // IN
								   const G_querypoint& p,       // IN
								   float_t sqr_one_plus_eps,    // IN
								   float_t sqr_bot,             // IN
								   std::size_t knn_maxsize,     // IN
								   float_t& knn_maxdist2,       // IN/OUT
								   float_t box_dist2,           // IN
								   kdcoords_t &box_offset,      // IN/OUT
								   const G_predicate& condition // IN
								   ) const {
      assert(node);
      assert(knn_it);
      assert(knn_dist2);
      assert(knn_maxsize>0);
      assert(knn_maxdist2>=0);
      assert(knn_size<=knn_maxsize);
      
      // Sort children
      const value_t&    node_value   = node->value();
      const std::size_t discriminant = node->discriminant();
      const float_t cut_diff = (float_t(p[discriminant]) - 
				float_t(node_value[discriminant]));
      const float_t cut_diff2 = cut_diff*cut_diff;
      const bool    lo_first = cut_diff < float_t(0);
      
      node_t* near_node = (lo_first ? node->lo_child() : node->hi_child());
      node_t* far_node  = (lo_first ? node->hi_child() : node->lo_child());

      // Always visit closer child first
      if (near_node) {
	internal_depth_first_k_approximately_nearest_neighbors_in(knn_size,
								  knn_it,
								  knn_dist2,
								  near_node, 
								  p, 
								  sqr_one_plus_eps,
								  sqr_bot,
								  knn_maxsize,
								  knn_maxdist2,
								  box_dist2, 
								  box_offset,
								  condition);
      }

      // If cutting discriminant is closer than (updated)
      // cut distance, test value point and possibly far nodes
      if (cut_diff2 < knn_maxdist2) {
	// Compute distance to query point
	float_t dist2 = distance_t::value(node_value, p, knn_maxdist2);

	// As in ANN, we count the number of distance computations
	++stat_knn_visited_count_;

	// Add to knn if closer than cut distance and matching criterion
	if (dist2 < knn_maxdist2 && condition(node_value)) {
	  // Insert
	  if (knn_size<knn_maxsize) {
	    // Insert at end and keep cut distance
	    knn_it[knn_size] = G_iter(node);
	    knn_dist2[knn_size] = dist2;
	    ++knn_size;
	  } else {
	    // Replace worst
	    assert(knn_size==knn_maxsize);
	    assert(dist2<knn_maxdist2);
	    knn_it[knn_size-1] = G_iter(node);
	    knn_dist2[knn_size-1] = dist2;
	  }
	  // Keep knn array sorted by increasing distance from query point
	  std::size_t idx = knn_size-1; 
	  while (idx>0 && knn_dist2[idx]<knn_dist2[idx-1]) {
	    std::swap(knn_it[idx], knn_it[idx-1]);
	    std::swap(knn_dist2[idx], knn_dist2[idx-1]);
	    --idx;
	  }
	  // Update max distance
	  if (knn_size==knn_maxsize) {
	    knn_maxdist2 = knn_dist2[knn_size-1];
	  }
	}

	// Check far node and visit further child if close enough
	if (far_node) {
	  // Incrementally update distance to bounding box
	  const float_t old_off = box_offset[discriminant];
	  box_dist2 += cut_diff2 - old_off*old_off;

	  // Visit only if better than minimum of current best and bounds
	  // overlap test threshold
	  if ((box_dist2*sqr_one_plus_eps < knn_maxdist2) &&
	      (box_dist2*sqr_one_plus_eps < sqr_bot)) {
	    box_offset[discriminant] = cut_diff;
	    internal_depth_first_k_approximately_nearest_neighbors_in(knn_size,
								      knn_it,
								      knn_dist2,
								      far_node, 
								      p, 
								      sqr_one_plus_eps,
								      sqr_bot,
								      knn_maxsize,
								      knn_maxdist2,
								      box_dist2, 
								      box_offset,
								      condition);
	    box_offset[discriminant] = old_off;
	  }
	}
      }
    }	
    
    template <typename G_querypoint, typename G_predicate, typename G_iter>
    void internal_prioritized_k_approximately_nearest_neighbors_in(std::size_t&    knn_size,    // OUT
								   G_iter         *knn_it,      // OUT
								   float_t        *knn_dist2,   // OUT
								   const G_querypoint& p,       // IN
								   float_t sqr_one_plus_eps,    // IN
								   float_t sqr_bot,             // IN
								   std::size_t knn_maxsize,     // IN
								   float_t& knn_maxdist2,       // IN/OUT
								   const G_predicate& condition // IN
								   ) const {
      assert(root_);
      assert(knn_it);
      assert(knn_dist2);
      assert(knn_maxsize>0);
      assert(knn_maxdist2>=0);
      assert(knn_size<=knn_maxsize);
      
      typedef detail::kdt_priority_search_queue_element<G_dimension, G_subvalue, G_value> kdt_priority_element_t;

#if 0
      sl::reservable_priority_queue<kdt_priority_element_t> queue;
      // Prealloc queue to avoid most memory allocation in inner loop.
      // This speeds-up parallel search by reducing lock contention.
      // It seems, however, to slightly reduce sequential performance because
      // of time needed to to a large allocation upfront even for searches
      // that require a small queue
      queue.reserve(1024);
#else
      // Using temporary stack allocation to avoid heap operations -- this is a bit
      // more convoluted (and borderline wrt C++ standard), but provides increased
      // preformance in most cases.

      typedef sl::tmp_arena<1024*sizeof(kdt_priority_element_t)> tmp_arena_t;
      typedef sl::arena_allocator<tmp_arena_t, kdt_priority_element_t> tmp_alloc_t;
      typedef std::vector<kdt_priority_element_t, tmp_alloc_t> tmp_vector_t;
      typedef std::priority_queue<kdt_priority_element_t, tmp_vector_t> tmp_priority_queue_t;
    
      tmp_arena_t tmp_buffer; // stack-allocated memory for first elements, heap-allocated ones for later ones
      tmp_priority_queue_t queue = tmp_priority_queue_t(std::less<kdt_priority_element_t>(),
							tmp_vector_t(tmp_alloc_t(tmp_buffer)));
#endif
      
      queue.push(kdt_priority_element_t(root_, kdcoords_t(), float_t(0.0)));
      bool done = false;
      std::size_t visited_node_count = 0;
      while (!done && !queue.empty()) {
	const kdt_priority_element_t& top = queue.top();
	node_t*       node       = top.node_;
	kdcoords_t    box_offset = top.box_offset_;
	float_t       box_dist2  = top.box_dist2_;
	queue.pop();

	done = 
	  (visited_node_count>priority_search_max_visited_count_) || 
	  (((box_dist2*sqr_one_plus_eps > knn_maxdist2) ||
	    (box_dist2*sqr_one_plus_eps > sqr_bot)));

	// Visit near branch and push far
	while (node && !done) {
	  // Access node and sort children
	  const value_t& node_value   = node->value();
	  const std::size_t discriminant = node->discriminant();
	  const float_t cut_diff = (float_t(p[discriminant]) - 
				    float_t(node_value[discriminant]));
	  const float_t cut_diff2 = cut_diff*cut_diff;
	  bool    lo_first = (cut_diff < float_t(0));

	  // If cutting discriminant is closer than (updated)
	  // cut distance, test value point (which is on cut plane) and possibly far nodes
	  if (cut_diff2 < knn_maxdist2) {
	    // Compute distance to query point
	    float_t dist2 = distance_t::value(node_value, p, knn_maxdist2);

	    // As in ANN, we count the number of distance computations
	    ++visited_node_count;
	    ++stat_knn_visited_count_;
	    
	    // Add to knn if closer than cut distance and matching criterion
	    if ((dist2 < knn_maxdist2) && condition(node_value)) {
	      // Insert
	      if (knn_size<knn_maxsize) {
                // Insert at end and keep cut distance
		knn_it[knn_size] = G_iter(node);
		knn_dist2[knn_size] = dist2;
		++knn_size;
	      } else {
		// Replace worst
                assert(knn_size==knn_maxsize);
		assert(dist2<knn_maxdist2);
		assert(dist2<knn_dist2[knn_size-1]);
		knn_it[knn_size-1] = G_iter(node);
		knn_dist2[knn_size-1] = dist2;
	      }
	      // Keep knn array sorted by increasing distance from query point
	      std::size_t idx = knn_size-1; 
	      while (idx>0 && knn_dist2[idx]<knn_dist2[idx-1]) {
		std::swap(knn_it[idx], knn_it[idx-1]);
		std::swap(knn_dist2[idx], knn_dist2[idx-1]);
		--idx;
	      }
	      // Update max distance
	      if (knn_size==knn_maxsize) {
		knn_maxdist2 = knn_dist2[knn_size-1];
              }
	    }

	    // Check far node and push it in queue if close enough
	    node_t* far_node  = (lo_first ? node->hi_child() : node->lo_child());
	    if (far_node && (cut_diff2 < knn_maxdist2)) {

              // FIXME -- Hack to avoid enqueuing single-point subtrees
              if (node->tree_size() == 2) {
		// We force descending on far tree that contains a single point
                lo_first = !lo_first;
              } else { 
	        // Incrementally update distance to bounding box
	        const float_t old_off = box_offset[discriminant];
	        const float_t far_dist2 = box_dist2 + cut_diff2 - old_off*old_off;
	      
	         // Visit only if better than minimum of current best and bounds
	         // overlap test threshold
	         if ((far_dist2*sqr_one_plus_eps < knn_maxdist2) &&
		     (far_dist2*sqr_one_plus_eps < sqr_bot)) {
		   box_offset[discriminant] = cut_diff;
		   queue.push(kdt_priority_element_t(far_node, box_offset, far_dist2));
		   box_offset[discriminant] = old_off;
                 }
              }
	    }
	  } // cut_diff < knn_maxdist2

	  // Continue following  near branch
	  node_t* near_node = (lo_first ? node->lo_child() : node->hi_child());
	  node = near_node;
 	} // while node and not done
      } // while queue not empty and not done
    } 
    
  public: // Single nearest neighbor

    template <typename G_querypoint, typename G_predicate>
    std::pair<const_iterator,float_t> find_nearest(const G_querypoint& p,
						   float_t knn_maxradius,
						   const G_predicate& condition) const {
      std::size_t    knn_size;
      const_iterator knn_it[1];
      float_t     knn_dist2;
      k_nearest_neighbors_in(knn_size,
			     knn_it,
			     &knn_dist2,
			     p,
			     std::size_t(1),
			     knn_maxradius,
			     condition);
      if (knn_size) {
	assert(knn_size==1);
	return std::make_pair(knn_it[0],knn_dist2);
      } else {
	return std::make_pair(end(),float_t(-1)); // FIXME dist2 invalid
      }
    }

    template <typename G_querypoint>
    std::pair<const_iterator,float_t> find_nearest(const G_querypoint& p,
						   float_t knn_maxradius) const {
      return find_nearest(p, knn_maxradius, kdalways());
    }

    template <typename G_querypoint>
    std::pair<const_iterator,float_t> find_nearest(const G_querypoint& p) const {
      return find_nearest(p, sqrt_maxd2_, kdalways());
    }

    template <typename G_querypoint, typename G_predicate>
    std::pair<iterator,float_t> find_nearest(const G_querypoint& p,
					     float_t knn_maxradius,
					     const G_predicate& condition) {
      std::size_t    knn_size;
      iterator       knn_it[1];
      float_t     knn_dist2;
      k_nearest_neighbors_in(knn_size,
			     knn_it,
			     &knn_dist2,
			     p,
			     std::size_t(1),
			     knn_maxradius,
			     condition);
      if (knn_size) {
	assert(knn_size==1);
	return std::make_pair(knn_it[0],knn_dist2);
      } else {
	return std::make_pair(end(),float_t(-1)); // FIXME dist2 invalid
      }
    }

    template <typename G_querypoint>
    std::pair<iterator,float_t> find_nearest(const G_querypoint& p,
					     float_t knn_maxradius) {
      return find_nearest(p, knn_maxradius, kdalways());
    }

    template <typename G_querypoint>
    std::pair<iterator,float_t> find_nearest(const G_querypoint& p) {
      return find_nearest(p, sqrt_maxd2_, kdalways());
    }
    
  public: // Single approximate nearest neighbor

    template <typename G_querypoint, typename G_predicate>
    std::pair<const_iterator,float_t> find_approximately_nearest(const G_querypoint& p,
								 float_t eps,
								 float_t bot,
								 float_t knn_maxradius,
								 const G_predicate& condition) const {
      std::size_t    knn_size;
      const_iterator knn_it[1];
      float_t     knn_dist2;
      k_approximately_nearest_neighbors_in(knn_size,
					   knn_it,
					   &knn_dist2,
					   p,
					   eps,
					   bot,
					   std::size_t(1),
					   knn_maxradius,
					   condition);
      if (knn_size) {
	assert(knn_size==1);
	return std::make_pair(knn_it[0],knn_dist2);
      } else {
	return std::make_pair(end(),float_t(-1)); // FIXME dist2 invalid
      }
    }

    template <typename G_querypoint>
    std::pair<const_iterator,float_t> find_approximately_nearest(const G_querypoint& p,
								 float_t eps,
								 float_t bot,
								 float_t knn_maxradius) const {
      return find_approximately_nearest(p, eps, bot, knn_maxradius, kdalways());
    }

    template <typename G_querypoint>
    std::pair<const_iterator,float_t> find_approximately_nearest(const G_querypoint& p,
								 float_t eps,
								 float_t bot) const {
      return find_approximately_nearest(p, eps, bot, sqrt_maxd2_, kdalways());
    }

    template <typename G_querypoint, typename G_predicate>
    std::pair<iterator,float_t> find_approximately_nearest(const G_querypoint& p,
							   float_t eps,
							   float_t bot,
							   float_t knn_maxradius,
							   const G_predicate& condition) {
      std::size_t    knn_size;
      iterator       knn_it[1];
      float_t     knn_dist2;
      k_approximately_nearest_neighbors_in(knn_size,
					   knn_it,
					   &knn_dist2,
					   p,
					   eps,
					   bot,
					   std::size_t(1),
					   knn_maxradius,
					   condition);
      if (knn_size) {
	assert(knn_size==1);
	return std::make_pair(knn_it[0],knn_dist2);
      } else {
	return std::make_pair(end(),float_t(-1)); // FIXME dist2 invalid
      }
    }

    template <typename G_querypoint>
    std::pair<iterator,float_t> find_approximately_nearest(const G_querypoint& p,
							   float_t eps,
							   float_t bot,
							   float_t knn_maxradius) {
      return find_approximately_nearest(p, eps, bot, knn_maxradius, kdalways());
    }

    template <typename G_querypoint>
    std::pair<iterator,float_t> find_approximately_nearest(const G_querypoint& p,
							   float_t eps,
							   float_t bot) {
      return find_approximately_nearest(p, eps, bot, sqrt_maxd2_, kdalways());
    }

  public: // random access

    /// iterator to the ith element of the tree. Logarithmic time 
    iterator ith(std::size_t i) {
      // Preorder iteration
      std::size_t N = count();
      if (i>=N) {
	return end();
      } else {
	node_t*      n = root_;
	std::size_t  k= 0;
	while (k!=i) {
	  assert(n);
	  assert(k<i);
	  node_t* n_lo_child = n->lo_child();
	  node_t* n_hi_child = n->hi_child();
	  std::size_t sz_lo_child = (n_lo_child ? n_lo_child->tree_size() : std::size_t(0));
	  // std::size_t sz_hi_child = (n_hi_child ? n_hi_child->tree_size() : std::size_t(0));
	  // Skip current node
	  k+=1;
	  if (i<k+sz_lo_child) {
	    // Choose lo child
	    n=n_lo_child;
	  } else {
	    // Skip lo child and move to hi child
	    k+=sz_lo_child;
	    n=n_hi_child;
	  }
	}
	return iterator(n);
      }
    }

    /// iterator to the ith element of the tree. Logarithmic time 
    const_iterator ith(std::size_t i) const {
      return const_cast<this_t*>(this)->ith(i);
    }
    
  public: // find

    /**
     * Finds first value having the same exact location.
     */
    iterator find(const_reference p) {
      return find_nearest(p, float_t(1), kdsamecoords<dimension,value_t>(p)).first;
    }

    /**
     * Finds first value having the same exact location.
     */
    const_iterator find(const_reference p) const {
      return const_cast<this_t*>(this)->find(p);
    }
    
    /**
     * Finds first value equal to p
     */
    iterator find_exact(const_reference p) {
      return find_nearest(p, float_t(1), kdsame<value_t>(p)).first;
    }

    /**
     * Finds first value equal to p.
     */
    const_iterator find_exact(const_reference p) const {
      return const_cast<this_t*>(this)->find_exact(p);
    }

  public: // Insert/delete

    /// Clear the tree
    void clear();

    /// Insert a new value in the tree. The value is copied
    void insert(const_reference p);

    /// Insert all elements from in_begin to in_end
    template <class G_input_iterator>
    void insert(G_input_iterator in_begin, G_input_iterator in_end) {
      for (G_input_iterator it = in_begin; it != in_end; ++it) {
	insert(*it);
      }
    }

    /**
     * Erase first value equal to p from tree if present.
     * Returns true iff erased.
     */
    bool erase(const_reference p) {
      iterator it = find(p);
      bool found = (it != end());
      if (found) {
	erase(it);
      }
      return found;
    }

    /**
     * Erase first value equal to v from tree if present.
     * Returns true iff erased.
     */
    bool erase_exact(const_reference p) {
      iterator it = find_exact(p);
      bool found = (it != end());
      if (found) {
	erase(it);
      }
      return found;
    }

    /// Erase node pointed by iterator from tree
    void erase(iterator const& it) {
      node_t* n = it.node_ptr();
      if (n) {
	node_t* n_parent = n->parent();
	// Join two children of deleted node
	node_t* n_lo_child = n->lo_child();
	node_t* n_hi_child = n->hi_child();
	std::size_t n_discriminant = n->discriminant();

	n->set_lo_child(0);
	n->set_hi_child(0);

	node_t* n_joined = join(n_lo_child, n_hi_child, n_discriminant); 
	if (!n_parent) {
	  root_ = n_joined;
	} else {
	  if (n_parent->lo_child() == n) {
	    n_parent->set_lo_child(n_joined);
	  } else {
	    assert(n_parent->hi_child() == n);
	    n_parent->set_hi_child(n_joined);
	  }
	}

	// Update subtree sizes up to root
	while (n_parent) {
	  assert(n_parent->tree_size()>0);
	  n_parent->set_tree_size(n_parent->tree_size()-1);
	  n_parent = n_parent->parent();
	}

	// Free node memory and update node count
	assert(n->is_root());
	assert(n->is_leaf());
	delete_nodes(n);
      }
    }
 	
  public: // Traversal

    const_iterator begin() const {
      // Preorder iteration
      return const_iterator(root_);
    }

    const_iterator end() const {
      return const_iterator();
    }
    
    iterator begin() {
      // Preorder iteration
      return iterator(root_);
    }

    iterator end() {
      return iterator();
    }

    /// Apply action a to each datum in arbitrary order
    template <class action_t>
    inline void apply(action_t& a) const {
      if (root_) root_->apply(a);
    }

  public: // Statistics

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
    bool invariant() const;

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

    void tree_depth_update_in(const node_t* node,
                              std::size_t node_depth,
                              std::size_t* min_depth,
                              std::size_t* max_depth) const;

    bool node_invariant(const node_t* node) const;

    bool is_lo(const node_t* node,
               std::size_t i,
               subvalue_t x_i) const;

    bool is_hi(const node_t* node,
               std::size_t i,
               subvalue_t x_i) const;

    std::size_t external_path_length(const node_t* node,
                                     std::size_t depth) const;

    std::size_t internal_path_length(const node_t* node,
                                     std::size_t depth) const;

    
  protected:

    void split(node_t* t,
               std::size_t i,
               subvalue_t x_i,
               node_t** t_lo,
               node_t** t_hi);

    node_t *join(node_t* a,
                 node_t* b,
                 std::size_t i);

    node_t *insert(node_t* t, node_t* x);

    node_t *insert_at_root(node_t* t, node_t* x);

  }; // class kdtree
    			 
} // namespace sl

// -- Inline implementation

namespace sl {

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  void kdtree<G_dimension, G_subvalue, G_value>::delete_nodes(node_t* n) {
    if (n) {
      node_t* old_lo_child = n->lo_child();
      n->set_lo_child(0);
      delete_nodes(old_lo_child);
      
      node_t* old_hi_child = n->hi_child();
      n->set_hi_child(0);
      delete_nodes(old_hi_child);
#if 0
      delete n; 
#else
      n->~node_t();
      node_allocator_.deallocate(n,1);
#endif
      n = 0; 
      --count_;
    }
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  void kdtree<G_dimension, G_subvalue, G_value>::clear() {
    delete_nodes(root_);
    root_ = 0;
    random_reset();
    assert(is_empty());
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  void kdtree<G_dimension, G_subvalue, G_value>::insert(const_reference p) {
#if 0
    node_t* x = new node_t(p); // FIXME Allocator
#else
    node_t* x = new (node_allocator_.allocate(1)) node_t(p);
#endif
    x->set_discriminant(random_discriminant());
    if (root_ == 0) {
      // First node
      root_ = x; 
    } else {
      root_ = insert(root_, x);
    }
    ++count_;
  }
    
  template <std::size_t G_dimension, class G_subvalue,class G_value>
  typename kdtree<G_dimension, G_subvalue, G_value>::node_t* kdtree<G_dimension, G_subvalue, G_value>::insert(node_t* t, node_t* x) {
    assert(x);
    assert(x->is_leaf());

    if (t == 0) {
      return x;
    } else {
      std::size_t n =  t->tree_size();

      const bool select_x_as_root = random_test(1, n+1);

      if (select_x_as_root) {
	return insert_at_root(t,x);
      } else {
	const std::size_t discriminant = t->discriminant();
	if (x->value()[discriminant] <= t->value()[discriminant]) {
	  t->set_lo_child(insert(t->lo_child(), x));
	} else {
	  t->set_hi_child(insert(t->hi_child(), x));
	}
	t->set_tree_size(n+1);
	return t;
      }
    }
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  typename kdtree<G_dimension, G_subvalue, G_value>::node_t* kdtree<G_dimension, G_subvalue, G_value>::insert_at_root(node_t* t, node_t* x) {
    assert(t);
    assert(x);
    assert(x->is_leaf());

    node_t *t_lo;
    node_t *t_hi;
    split(t,x->discriminant(),x->discriminant_plane_constant(),&t_lo,&t_hi);
    if (t_lo) {
      x->set_lo_child(t_lo);
      x->set_tree_size(x->tree_size() + t_lo->tree_size());
    }
    if (t_hi) {
      x->set_hi_child(t_hi);
      x->set_tree_size(x->tree_size() + t_hi->tree_size());
    }
    return x;
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  void kdtree<G_dimension, G_subvalue, G_value>::split(node_t* t,
						   std::size_t i,
						   subvalue_t x_i,
						   node_t** t_lo,
						   node_t** t_hi) {
    assert(i<dimension);
    assert(t_lo);
    assert(t_hi);
    
    if (t == 0) {
      *t_lo = 0;
      *t_hi = 0;
    } else {
      std::size_t j = t->discriminant();
      subvalue_t y_i = t->value()[i];
      if (i==j) {
	if (y_i <= x_i) {
	  node_t *t_hi_lo;
	  node_t *t_hi_hi;
	  split(t->hi_child(), i, x_i, &t_hi_lo, &t_hi_hi);
	  t->set_hi_child(t_hi_lo);
	  *t_lo = t;
	  *t_hi = t_hi_hi;
	} else {
	  assert(y_i > x_i);
	  node_t *t_lo_lo;
	  node_t *t_lo_hi;
	  split(t->lo_child(), i, x_i, &t_lo_lo, &t_lo_hi);
	  t->set_lo_child(t_lo_hi);
	  *t_lo = t_lo_lo;
	  *t_hi = t;
	}
      } else {
	assert(i!=j);
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

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  typename kdtree<G_dimension, G_subvalue, G_value>::node_t *kdtree<G_dimension, G_subvalue, G_value>::join(node_t* a,
												    node_t* b,
												    std::size_t i) {
    assert(i<dimension);

    if (a == 0) {
      return b; 
    } else if (b == 0) {
      return a;
    } else {
      std::size_t n = a->tree_size();
      std::size_t m = b->tree_size();

      const bool select_a_as_root = random_test(n, n+m); // RANDOMIZED n>m
      // const bool select_a_as_root = n>m; // NON RANDOMIZED

      if (select_a_as_root) {
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

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  void kdtree<G_dimension, G_subvalue, G_value>::tree_depth_update_in(const node_t* node,
                                                                   std::size_t node_depth,
                                                                   std::size_t* min_depth,
                                                                   std::size_t* max_depth) const {
    assert(node);
    assert(min_depth);
    assert(max_depth);
    if (node->is_leaf()) {
      if (node_depth < *min_depth) *min_depth = node_depth;
      if (node_depth > *max_depth) *max_depth = node_depth;
    } else {
      if (node->lo_child()) tree_depth_update_in(node->lo_child(), node_depth+1, min_depth, max_depth);
      if (node->hi_child()) tree_depth_update_in(node->hi_child(), node_depth+1, min_depth, max_depth);
    }
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  void kdtree<G_dimension, G_subvalue, G_value>::tree_depth_in(std::size_t* min_depth,
							       std::size_t* max_depth) const {
    assert(min_depth);
    assert(max_depth);
    if (root_) {
      *min_depth = count()+1;
      *max_depth = 0;
      tree_depth_update_in(root_, 1, min_depth, max_depth);
    } else {
      *min_depth = 0;
      *max_depth = 0;
    }      
  }  
  
  template <std::size_t G_dimension, class G_subvalue,class G_value>
  bool kdtree<G_dimension, G_subvalue, G_value>::is_lo(const node_t* node,
                                                    std::size_t i,
                                                    subvalue_t x_i) const {
    bool result = true;
    if (node) {
      result = result && node->value()[i] <= x_i;
      result = result && is_lo(node->lo_child(), i, x_i);
      result = result && is_lo(node->hi_child(), i, x_i);
    }
    return result;
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  bool kdtree<G_dimension, G_subvalue, G_value>::is_hi(const node_t* node,
                                                    std::size_t i,
                                                    subvalue_t x_i) const {
    bool result = true;
    if (node) {
      result = result && node->value()[i] > x_i;
      result = result && is_hi(node->lo_child(), i, x_i);
      result = result && is_hi(node->hi_child(), i, x_i);
    }
    return result;
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  bool kdtree<G_dimension, G_subvalue, G_value>::invariant() const {
    bool result = true;

    const std::size_t tree_sz = (root_ == 0) ? 0 : root_->tree_size();

    if (result) {
      result = (count() == tree_sz);

      if (!result) SL_TRACE_OUT(0) << "kdtree::invariant() : count() [" << count() << "] != tree_sz [" << tree_sz << "]!" << std::endl;
    }

    if (result) {
      result = node_invariant(root_);
    }

    return result;
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  bool kdtree<G_dimension, G_subvalue, G_value>::node_invariant(const node_t* node) const {
    bool result = true;
    if (node) {
      const node_t* l = node->lo_child();
      const node_t* h = node->hi_child();
      
      if (result) {
        const std::size_t l_sz = (l==0) ? 0 : l->tree_size();
        const std::size_t h_sz = (h==0) ? 0 : h->tree_size();
      
        result = (node->tree_size() == (1 + l_sz + h_sz));

        if (!result) SL_TRACE_OUT(0) << "kdtree::invariant(node) : inconsistent tree size !" << std::endl;
      }

      if (result) {
	result = ((node == root_) == (node->parent() == 0));
        if (!result) SL_TRACE_OUT(0) << "kdtree::invariant(node) : root" << std::endl;
      }

      if (result) {
	if (node->lo_child()) {
	  result = (node->lo_child()->parent() == node);
	}
        if (!result) SL_TRACE_OUT(0) << "kdtree::invariant(node) : inconsistent lo child parent" << std::endl;
      }

      if (result) {
	if (node->hi_child()) {
	  result = (node->hi_child()->parent() == node);
	}
        if (!result) SL_TRACE_OUT(0) << "kdtree::invariant(node) : inconsistent hi child parent" << std::endl;
      }

      if (result) {
        result = is_lo(node->lo_child(), node->discriminant(), node->discriminant_plane_constant());

        if (!result) SL_TRACE_OUT(0) << "kdtree::invariant(node) : lo split discriminant test failed !" << std::endl;
      }

      if (result) {
        result = is_hi(node->hi_child(), node->discriminant(), node->discriminant_plane_constant());

        if (!result) SL_TRACE_OUT(0) << "kdtree::invariant(node) : hi split discriminant test failed !" << std::endl;
      }
	
      if (result) {
        // Check children invariant
        result = result && node_invariant(l);
        result = result && node_invariant(h);
      }
    }
    return result;
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  std::size_t kdtree<G_dimension, G_subvalue, G_value>::external_path_length() const {
    return external_path_length(root_, 1);
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  std::size_t kdtree<G_dimension, G_subvalue, G_value>::internal_path_length() const {
    return internal_path_length(root_, 1);
  }

  template <std::size_t G_dimension, class G_subvalue,class G_value>
  std::size_t kdtree<G_dimension, G_subvalue, G_value>::external_path_length(const node_t* node,
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


  template <std::size_t G_dimension, class G_subvalue,class G_value>
  std::size_t kdtree<G_dimension, G_subvalue, G_value>::internal_path_length(const node_t* node,
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

}

#endif
