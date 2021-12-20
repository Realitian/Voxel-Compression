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
#ifndef SL_RANDOMIZED_KD_TREE_HPP
#define SL_RANDOMIZED_KD_TREE_HPP

#include <sl/kd_tree.hpp>
#include <sl/random.hpp>

/////////////////////////////// DEPRECATED CLASS ////////////////////////
#warning \
Class sl::randomized_kd_tree from sl/randomized_kd_tree.hpp is obsolete \
and will be removed from future sl releases. Please use class sl::kdtree \
defined in sl/kdtree.hpp.
/////////////////////////////// DEPRECATED CLASS ////////////////////////

namespace sl {

  /**
   *  Randomized kd-trees of dimension G_dimension containing
   *  entities of type G_data with coordinates of
   *  type G_value. The class G_data
   *  is constrained to export G_value operator[](std::size_t i)
   *  to provide read access to coordinates 0 to dimension.
   *  The insertion procedure uses a randomized approach in
   *  order to keed the tree (probabilistically) balanced
   *  independently from the data.
   *
   *  Based on:
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
   */
  template <std::size_t G_dimension, 
            class G_value,
	    class G_data>
  class randomized_kd_tree: public kd_tree<G_dimension, G_value, G_data> {
  public:
    enum { dimension = G_dimension };

    typedef kd_tree<G_dimension, G_value, G_data> super_t;
    typedef randomized_kd_tree<G_dimension, G_value, G_data> this_t;
    typedef G_value value_t;
    typedef G_data  data_t;

  protected: // Random number generation

    typedef detail::kd_coords<G_dimension, G_value> kd_coords_t;
    typedef detail::kd_node<G_dimension, G_value, G_data> node_t;

    random::std_irng_t rng_;

    /// Reset the pseudo random number generator
    inline void random_reset() {
      rng_.set_seed(); // default
    }

    /// An unsigned random number between l and h 
    inline unsigned int random_unsigned(unsigned int l, unsigned int h) {
      SL_REQUIRE("Good parameter", h >= l);
      return l + rng_.value() % (h-l+1);
    }
      
    /// A random discriminant
    inline std::size_t random_discriminant() {
      return (std::size_t)(random_unsigned(0,dimension-1));
    }
      
    /// A random test returning true with probability a/b
    inline bool random_test(unsigned int a, unsigned int b) {
      SL_REQUIRE("Good parameter", b>0);
      SL_REQUIRE("Good parameter", b>=a);
      return random_unsigned(1,b) <= a;
    } 

  public:

    inline randomized_kd_tree(const std::size_t memory_pool_chunk_first_count = 0,
                              const std::size_t memory_pool_chunk_grow_factor = 1)
    : super_t( memory_pool_chunk_first_count,memory_pool_chunk_grow_factor) {
      random_reset();
    }

    /// Clear the tree
    virtual void clear() {
      super_t::clear();
      random_reset();
    }

    /// Insert a new data point in the tree. The point is copied.
    virtual void insert(const data_t& p) {
      node_t* x = pool_new_node(p);
      x->set_discriminant(random_discriminant());
      if ((this->root_) == 0) {
	// First node
	(this->root_) = x;
	// Init bbox
	for (std::size_t i = 0; i<dimension; ++i) {
	  (this->bbox_lo_)[i] = p[i];
	  (this->bbox_hi_)[i] = p[i];
	}
      } else {
	// Update bbox
	for (std::size_t i = 0; i<dimension; ++i) {
	  if (p[i]<(this->bbox_lo_)[i]) (this->bbox_lo_)[i] = p[i];
	  if (p[i]>(this->bbox_hi_)[i]) (this->bbox_hi_)[i] = p[i];
	}
	
	(this->root_) = randomized_insert((this->root_), x);
      }
      ++(this->count_);
    }

  protected:
    
    node_t *randomized_insert(node_t* t, 
			      node_t* x) {
      SL_REQUIRE("Good parameter", x);
      SL_REQUIRE("Good parameter", x->is_leaf());

      if (t == 0) {
	return x;
      } else {
	std::size_t n =  t->tree_size();
#if 1
	bool insert_x_at_root = random_test(1, n+1);
#else
	std::size_t n_lo =  (t->lo_child() == 0) ? 0 : t->lo_child()->tree_size();
	std::size_t n_hi = n-1-n_lo;

	bool insert_x_at_root = n>10 && (n_lo/(double(n-1)) < 0.3) || (n_hi/(double(n-1)) < 0.3);
#endif
	if (insert_x_at_root) {
	  return randomized_insert_at_root(t,x);
	} else {
	  const std::size_t discriminant = t->discriminant();
	  if (x->data()[discriminant] <= t->data()[discriminant]) {
	    t->set_lo_child(randomized_insert(t->lo_child(), x));
	  } else {
	    t->set_hi_child(randomized_insert(t->hi_child(), x));
	  }
	  t->set_tree_size(n+1);
	  return t;
	}
      }
    }

    node_t* randomized_insert_at_root(node_t* t, node_t* x) {
      SL_REQUIRE("Good parameter", t);
      SL_REQUIRE("Good parameter", x);
      SL_REQUIRE("Good parameter", x->is_leaf());

      node_t *t_lo;
      node_t *t_hi;
      randomized_split(t,x->discriminant(),x->discriminant_plane_constant(),&t_lo,&t_hi);
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

    void randomized_split(node_t* t,
			  std::size_t i,
			  value_t x_i,
			  node_t** t_lo,
			  node_t** t_hi) {
      SL_REQUIRE("Good parameter", i<dimension);
      SL_REQUIRE("Good parameter", t_lo);
      SL_REQUIRE("Good parameter", t_hi);

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
	    randomized_split(t->hi_child(), i, x_i, &t_hi_lo, &t_hi_hi);
	    t->set_hi_child(t_hi_lo);
	    *t_lo = t;
	    *t_hi = t_hi_hi;
	  } else {
	    SL_CHECK("Ordered", y_i > x_i);
	    node_t *t_lo_lo;
	    node_t *t_lo_hi;
	    randomized_split(t->lo_child(), i, x_i, &t_lo_lo, &t_lo_hi);
	    t->set_lo_child(t_lo_hi);
	    *t_lo = t_lo_lo;
	    *t_hi = t;
	  }
	} else {
	  SL_CHECK("Different dimension", i!=j);
	  node_t *t_lo_lo;
	  node_t *t_lo_hi;
	  randomized_split(t->lo_child(), i, x_i, &t_lo_lo, &t_lo_hi);
	  node_t *t_hi_lo;
	  node_t *t_hi_hi;
	  randomized_split(t->hi_child(), i, x_i, &t_hi_lo, &t_hi_hi);
	  if (y_i <= x_i) {
	    t->set_lo_child(t_lo_lo);
	    t->set_hi_child(t_hi_lo);
	    *t_lo = t;
	    *t_hi = randomized_join(t_lo_hi,t_hi_hi,j);
	  } else {
	    *t_lo = randomized_join(t_lo_lo,t_hi_lo,j);
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

    node_t *randomized_join(node_t* a,
			    node_t* b,
			    std::size_t i) {
      SL_REQUIRE("Good parameter", i<dimension);

      if (a == 0) {
	return b; 
      } else if (b == 0) {
	return a;
      } else {
	std::size_t n = a->tree_size();
	std::size_t m = b->tree_size();
#if 1
	if (random_test(n,n+m)) {
#else
	if (n>m) {
#endif
	  // select a as root of the join
	  std::size_t j_a = a->discriminant();
	  if (i==j_a) {
	    a->set_hi_child(randomized_join(a->hi_child(), b, i));
	  } else {
	    node_t* b_lo;
	    node_t* b_hi;
	    randomized_split(b, j_a, a->discriminant_plane_constant(), &b_lo, &b_hi);
	    a->set_lo_child(randomized_join(a->lo_child(),b_lo, i));
	    a->set_hi_child(randomized_join(a->hi_child(),b_hi, i));
	  }
	  a->set_tree_size(n+m);
	  return a;
	} else {
	  // select b as root of the join
	  std::size_t j_b = b->discriminant();
	  if (i==j_b) {
	    b->set_lo_child(randomized_join(a, b->lo_child(), i));
	  } else {
	    node_t* a_lo;
	    node_t* a_hi;
	    randomized_split(a, j_b, b->discriminant_plane_constant(), &a_lo, &a_hi);
	    b->set_lo_child(randomized_join(a_lo, b->lo_child(), i));
	    b->set_hi_child(randomized_join(a_hi, b->hi_child(), i));
	  }
	  b->set_tree_size(n+m);
	  return b;
	}  
      }
    }

  }; // class randomized_kd_tree

} // namespace sl

#endif
