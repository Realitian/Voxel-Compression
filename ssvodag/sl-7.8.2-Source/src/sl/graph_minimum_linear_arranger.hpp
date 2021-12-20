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
#ifndef SL_GRAPH_MINIMUM_LINEAR_ARRANGER_HPP
#define SL_GRAPH_MINIMUM_LINEAR_ARRANGER_HPP

#include <sl/random.hpp>
#include <cassert>
#include <vector>
#include <iostream>

namespace sl {

  class static_weighted_graph {
  protected:
    std::vector<std::size_t> node_degree_;
    std::vector<std::size_t> node_adjacency_;
    std::vector<std::size_t> edge_target_;
    std::vector<float>       edge_weight_;

  protected:

    void rebuild_adjacency();
    
  public:

    inline static_weighted_graph() {}
    inline ~static_weighted_graph() {}

    void clear();

    void init(const std::vector<std::size_t>& graph_node_degree,
	      const std::vector<std::size_t>& graph_edge_target,
	      const std::vector<float>&       graph_edge_weight);
    
    void init(const std::vector<std::size_t>& graph_node_degree,
	      const std::vector<std::size_t>& graph_edge_target);

    void load_gra(std::istream& is);
    void load_gra(const std::string& filename);

    void save_gra(std::ostream& os) const;
    void save_gra(const std::string& filename) const;

  public:
 
    inline std::size_t node_count() const { return node_degree_.size(); }

    inline std::size_t edge_count() const { return edge_target_.size(); }

    inline std::size_t node_degree(std::size_t u) const { assert(u<node_count()); return node_degree_[u];}

    inline std::size_t node_adjacency(std::size_t u) const { assert(u<node_count()); return node_adjacency_[u]; }

    inline std::size_t edge_target(std::size_t e) const { assert(e<edge_count()); return edge_target_[e]; }

    inline float       edge_weight(std::size_t e) const { assert(e<edge_count()); return edge_weight_[e]; }

  public:

    std::size_t max_degree() const;
    
    void depth_first_visit_in(std::vector<bool>& visited, std::size_t u) const;
				
    bool are_connected(std::size_t u, std::size_t v) const;

    float edge_weight(std::size_t u, std::size_t v) const;
    
    std::size_t component_count() const;

  public:
       
    void make_undirected();

    void reorder(const std::vector<std::size_t>& ordering);
  };
  
  /**
   *  Approximate solution of We concentrate here of the weighted minimum
   *  linear arrangement problem (MinLA), which consists of placing the N
   *  vertices of a graph at positions 1 . . . n on a line, so as to
   *  minimize the sum of the (weighted) edge lengths.
   *
   *  Implementation loosely based on:
   *
   *  Koren, Yehuda, and David Harel. "A multi-scale algorithm for the linear
   *  arrangement problem." Graph-Theoretic Concepts in Computer Science.
   *  Springer Berlin Heidelberg, 2002.
   *
   */
  class graph_minimum_linear_arranger {
  public:
    typedef graph_minimum_linear_arranger this_t;

    typedef enum { CANONICAL, DFS, BFS, GREEDY } init_kind_t;
  protected:

    init_kind_t init_kind_;
    std::size_t restart_k_;
    std::size_t relax_k1_;
    std::size_t relax_k2_;
    float       relax_beta_;
    std::size_t localopt_n_;
    std::size_t localopt_k_;
    mutable random::std_irng_t irnd_;
    
  public:

    graph_minimum_linear_arranger();
    ~graph_minimum_linear_arranger();

    void arrangement_in(std::vector<std::size_t>&        ordering,
			const sl::static_weighted_graph& graph);

    float cost(const std::vector<std::size_t>&  ordering,
	       const sl::static_weighted_graph& graph) const;
 
    void reciprocal_in(std::vector<std::size_t>&  out_ordering,
		       const std::vector<std::size_t>&  in_ordering) const;
    
  protected: // Internals

    void initialize_in(std::vector<std::size_t>&        ordering,
		       const sl::static_weighted_graph& graph);

    void canonical_initialize_in(std::vector<std::size_t>&        ordering,
				 const sl::static_weighted_graph& graph);

    void dfs_initialize_in(std::vector<std::size_t>&        ordering,
			   const sl::static_weighted_graph& graph);

    void bfs_initialize_in(std::vector<std::size_t>&        ordering,
			   const sl::static_weighted_graph& graph);

    void greedy_initialize_in(std::vector<std::size_t>&        ordering,
			      const sl::static_weighted_graph& graph);
    
    void greedy_improve(std::vector<std::size_t>&        ordering,
			const sl::static_weighted_graph& graph);
   
    void coarsen_in(std::vector<std::size_t>&        coarse_ordering,
		    sl::static_weighted_graph&       coarse_graph,
		    std::vector<std::size_t>&        fine_to_coarse_vertex_index,
		    std::vector<std::pair<std::size_t,std::size_t> >& coarse_to_fine_vertex_index,
		    const std::vector<std::size_t>&  fine_ordering,
		    const sl::static_weighted_graph& fine_graph);

    void uncoarsen_in(std::vector<std::size_t>&        fine_ordering,
		      const sl::static_weighted_graph& fine_graph,
		      const std::vector<std::size_t>&  coarse_ordering,
		      const std::vector<std::size_t>&             fine_to_coarse_vertex_index,
		      const std::vector<std::pair<std::size_t,std::size_t> >& coarse_to_fine_vertex_index);

    float relaxed_l1_irwls(std::size_t u,
			   const std::vector< std::pair<float,std::size_t> >& placement,
			   const sl::static_weighted_graph&  graph) const;

    float relaxed_l1_median(std::size_t u,
			 const std::vector< std::pair<float,std::size_t> >& placement,
			 const sl::static_weighted_graph&  graph) const;

    float relaxed_l2(std::size_t u,
			 const std::vector< std::pair<float,std::size_t> >& placement,
			 const sl::static_weighted_graph&  graph) const;
    
    void relax(std::vector<std::size_t>&        graph_node_ordering,
	       const sl::static_weighted_graph& graph,
	       float beta, 
	       std::size_t sweep_count);

    void optimize_small_local_sequence(std::vector<std::size_t>&        ordering,
				       std::vector<std::size_t>&        sorted_nodes, // reciprocal of ordering
				       const sl::static_weighted_graph& graph,
				       std::size_t seq_start, std::size_t seq_sz);
    
    void v_cycle(std::vector<std::size_t>&        ordering,
		 const sl::static_weighted_graph& graph);

    void local_optimization_sweep(std::vector<std::size_t>&        ordering,
				  const sl::static_weighted_graph& graph,
				  std::size_t sweep_count);

  };
    
} // namespace sl

#endif
