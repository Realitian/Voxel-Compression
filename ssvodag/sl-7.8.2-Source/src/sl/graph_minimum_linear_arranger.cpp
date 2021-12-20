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
// ============ FIXME
// #undef NDEBUG
// ============ FIXME

#include <sl/graph_minimum_linear_arranger.hpp>
#include <sl/permutations.hpp>
#include <sl/math.hpp>
#include <fstream>
#include <vector>
#include <cassert>
#include <map>
#include <set>
#include <algorithm>

// ----------------------------------------------------------------------------------
// Static weighted graph
// ----------------------------------------------------------------------------------

namespace sl {

#ifndef NDEBUG

  static void DBG_DUMP(const std::string& msg, const std::vector<std::size_t>& o,
		       std::size_t seq_start, std::size_t seq_sz) {
    std::cerr << msg << "[";
    for (std::size_t i=seq_start; i<seq_start+seq_sz; ++i) {
      std::cerr << " "<< o[i];
    }
    std::cerr << " ]";
  }
  
  static void DBG_DUMP(const std::string& msg, const std::vector<std::size_t>& o) {
    DBG_DUMP(msg, o, 0, o.size());
  }

  static bool DBG_CHECK_UNDIRECTED(const std::string& msg, const sl::static_weighted_graph& graph) {
    bool result=true;
    
    const std::size_t N = graph.node_count();
    for (std::size_t u=0; u<N; ++u) {
      const std::size_t e_bgn = graph.node_adjacency(u);
      const std::size_t e_end = e_bgn + graph.node_degree(u);
      for (std::size_t e=e_bgn; e<e_end; ++e) {
	const std::size_t v=graph.edge_target(e);
	assert(graph.are_connected(u,v));
	if (!graph.are_connected(v,u)) {
	  std::cerr << msg << "!!! not connected(" << u << v << ")!" << std::endl;
	  result=false;
	} else if (graph.edge_weight(u,v) != graph.edge_weight(v,u)) {
	  std::cerr << msg << "!!! weight(" << u << v << ")=" << graph.edge_weight(u,v) << " !=" << graph.edge_weight(v,u) << std::endl;
	  result=false;
	}
      }
    }
    return result;
  }
#endif
  

			
  void static_weighted_graph::init(const std::vector<std::size_t>& graph_node_degree,
				   const std::vector<std::size_t>& graph_edge_target,
				   const std::vector<float>&       graph_edge_weight) {
    assert(graph_edge_target.size() == graph_edge_weight.size());
    
    node_degree_ = graph_node_degree;
    edge_target_ = graph_edge_target;
    edge_weight_ = graph_edge_weight;
    rebuild_adjacency();
  }

  void static_weighted_graph::clear() {
    node_degree_.clear();
    node_adjacency_.clear();
    edge_target_.clear();
    edge_weight_.clear();
  }
  
  void static_weighted_graph::rebuild_adjacency() {
    const std::size_t N = node_degree_.size();
    node_adjacency_.resize(N);
    std::size_t edge_offset=0;
    for (std::size_t u=0; u<N; ++u) {
      assert(edge_offset<edge_count());
      node_adjacency_[u] = edge_offset;
      edge_offset+= node_degree(u);
    }
  }

  void static_weighted_graph::init(const std::vector<std::size_t>& graph_node_degree,
				   const std::vector<std::size_t>& graph_edge_target) {
    const std::size_t N = graph_node_degree.size();
    std::vector<float> graph_edge_weight(N, 1.0f);
    init(graph_node_degree, graph_edge_target, graph_edge_weight);
  }

  void static_weighted_graph::load_gra(std::istream& is) {
    std::size_t N;
    std::size_t M;
    is >> N;
    is >> M;
    node_degree_.resize(N);
    M=0;
    for (std::size_t i=0; i<N; ++i) {
      is >> node_degree_[i];
      M+= node_degree_[i];
    }
    edge_target_.resize(M);
    for (std::size_t i=0; i<M; ++i) {
      is >> edge_target_[i];
    }
    edge_weight_.resize(M);
    for (std::size_t i=0; i<M; ++i) {
      edge_weight_[i] = 1.0f;
    }
    rebuild_adjacency();
  }
  
  void static_weighted_graph::load_gra(const std::string& filename) {
    clear();
    std::ifstream ifs(filename.c_str());
    if (ifs.good()) {
      load_gra(ifs);
    }
  }

  void static_weighted_graph::save_gra(std::ostream& os) const {
    const std::size_t N = node_count();
    const std::size_t M = edge_count();
    os << N << std::endl;
    os << M << std::endl;
    for (std::size_t i=0; i<N; ++i) {
      os << " " << node_degree_[i];
    }
    os << std::endl;
    for (std::size_t i=0; i<M; ++i) {
      os << " " << edge_target_[i];
    }
    os << " -1" << std::endl;
  }
  
  void static_weighted_graph::save_gra(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs.good()) {
      save_gra(ofs);
    }
  }

  std::size_t static_weighted_graph::max_degree() const {
    std::size_t result=0;

    const std::size_t N = node_count();
    for (std::size_t u=0; u<N; ++u) {
      result = std::max(result, node_degree(u));
    }
    return result;
  }

  bool static_weighted_graph::are_connected(std::size_t u, std::size_t v) const {
    bool result=false;

    const std::size_t e_bgn = node_adjacency(u);
    const std::size_t e_end = e_bgn + node_degree(u);
    for (std::size_t e=e_bgn; e<e_end && !result; ++e) {
      result = (edge_target(e) == v);
    }

    return result;
  }
  
  float static_weighted_graph::edge_weight(std::size_t u, std::size_t v) const {
    const std::size_t e_bgn = node_adjacency(u);
    const std::size_t e_end = e_bgn + node_degree(u);
    for (std::size_t e=e_bgn; e<e_end; ++e) {
      if (edge_target(e) == v) return edge_weight(e);
    }
    return 0.0f;
  }
 
  void static_weighted_graph::depth_first_visit_in(std::vector<bool>& visited, std::size_t u) const {
    assert(visited.size() == node_count());
    assert(u<node_count());

    visited[u]=true;
    const std::size_t e_bgn = node_adjacency(u);
    const std::size_t e_end = e_bgn + node_degree(u);
    for (std::size_t e=e_bgn; e<e_end; ++e) {
      const std::size_t v=edge_target(e);
      if (!visited[v]) {
	depth_first_visit_in(visited, v);
      }
    }
  }

  std::size_t static_weighted_graph::component_count() const {
    const std::size_t N = node_count();

    std::vector<bool> visited(N, false);

    std::size_t result=0;
    for (std::size_t u=0; u<N; ++u) {
      if (!visited[u]) {
	++result;
	depth_first_visit_in(visited, u);
      }
    }
    return result;
  }

  void static_weighted_graph::make_undirected() {
    // Construct undirected edge map
    const std::size_t N = node_count();
    std::vector< std::map<std::size_t, float> > edge_weight_map(N);

    std::size_t M = 0;
    for (std::size_t u=0; u<N; ++u) {
      const std::size_t e_bgn = node_adjacency(u);
      const std::size_t e_end = e_bgn + node_degree(u);
      for (std::size_t e=e_bgn; e<e_end; ++e) {
	assert(e<edge_count());
	const std::size_t v=edge_target(e);
	assert(v<node_count());
	const float w_uv = edge_weight(e);

	if (edge_weight_map[u].find(v) == edge_weight_map[u].end()) {
	  edge_weight_map[u][v] = w_uv;
	  ++M;
	} else {
	  edge_weight_map[u][v] += w_uv;
	}	  
	if (edge_weight_map[v].find(u) == edge_weight_map[v].end()) {
	  edge_weight_map[v][u] = w_uv;
	  ++M;
	} else {
	  edge_weight_map[v][u] += w_uv;
	}	  
      }
    }

    // Reconstruct graph from undirected edge map
    node_degree_.resize(N);
    node_adjacency_.resize(N);
    edge_target_.resize(M);
    edge_weight_.resize(M);
    assert(node_count() == N);
    assert(edge_count() == M);
    
    std::size_t edge_offset=0;
    for (std::size_t u=0; u<N; ++u) {
      assert(u<node_count());
      
      const std::size_t M_u = edge_weight_map[u].size();
      node_degree_[u] = M_u;
      node_adjacency_[u] = edge_offset;
      assert(edge_offset<edge_count());
      assert(edge_offset+M_u<=edge_count());
      std::size_t k=0;
      for (std::map<std::size_t, float>::iterator it = edge_weight_map[u].begin();
	   it != edge_weight_map[u].end();
	   ++it) {
	std::size_t e = edge_offset + k;
	assert(e<edge_count());
	edge_target_[e] = it->first;
 	edge_weight_[e] = 0.5f*it->second;
	++k;
      }
      edge_offset += M_u;
    }

  }
  
  void static_weighted_graph::reorder(const std::vector<std::size_t>& ordering) {
    const std::size_t N = node_count();
    const std::size_t M = edge_count();
    
    std::vector<std::size_t> new_node_degree(N);
    std::vector<std::size_t> new_node_adjacency(N);
    std::vector<std::size_t> new_edge_target(M);
    std::vector<float>       new_edge_weight(M);

    // Rebuild node data
    for (std::size_t u=0; u<N; ++u) {
      const std::size_t new_u = ordering[u];
      assert(new_u<node_count());
      new_node_degree[new_u] = node_degree(u);
    }
    std::size_t new_edge_offset=0;
    for (std::size_t new_u=0; new_u<N; ++new_u) {
      new_node_adjacency[new_u] = new_edge_offset;
      assert(new_edge_offset<edge_count());
      assert(new_edge_offset+new_node_degree[new_u]<=edge_count());
      new_edge_offset += new_node_degree[new_u];
    }

    // Rebuild edge data
    for (std::size_t u=0; u<N; ++u) {
      const std::size_t new_u = ordering[u];
      assert(new_u<node_count());
      
      const std::size_t e_bgn = node_adjacency(u);
      const std::size_t e_end = e_bgn + node_degree(u);

      const std::size_t new_e_bgn = new_node_adjacency[new_u];
      
      for (std::size_t e=e_bgn; e<e_end; ++e) {
	assert(e<edge_count());
	
	const std::size_t new_e = new_e_bgn + (e-e_bgn);
	assert(new_e<edge_count());
	       
	const std::size_t v=edge_target(e);
	assert(v<node_count());

	const std::size_t new_v = ordering[v];
	assert(new_v<node_count());

	const float w_uv = edge_weight(e);
	
	new_edge_target[new_e] = new_v;
	new_edge_weight[new_e] = w_uv;
      }
    }

    std::swap(node_degree_, new_node_degree);
    std::swap(node_adjacency_, new_node_adjacency);
    std::swap(edge_target_, new_edge_target);
    std::swap(edge_weight_, new_edge_weight);
  }

}

// ----------------------------------------------------------------------------------
// MINLA
// ----------------------------------------------------------------------------------

namespace sl {

  graph_minimum_linear_arranger::graph_minimum_linear_arranger() {
    init_kind_  = GREEDY;
    restart_k_  = 8;
    relax_k1_   = 40;
    relax_k2_   = 10;
    relax_beta_ = 1.0f;
    localopt_n_ = 5;
    localopt_k_ = 6;
  }
  
  graph_minimum_linear_arranger::~graph_minimum_linear_arranger() {
  }
    
  float graph_minimum_linear_arranger::cost(const std::vector<std::size_t>& ordering,
					    const sl::static_weighted_graph& graph) const {
    const std::size_t N = graph.node_count();

    float result = 0.0f;
    for (std::size_t u=0; u<N; ++u) {
      const std::size_t e_bgn = graph.node_adjacency(u);
      const std::size_t e_end = e_bgn + graph.node_degree(u);
      for (std::size_t e=e_bgn; e<e_end; ++e) {
	const std::size_t v=graph.edge_target(e);
	const float w_uv = graph.edge_weight(e);
    	result += w_uv*sl::abs(float(ordering[u])-float(ordering[v]));
      }
    }
    return 0.5f*result; // We sum over directed edges
  }

  void graph_minimum_linear_arranger::reciprocal_in(std::vector<std::size_t>&  out_ordering,
						    const std::vector<std::size_t>&  in_ordering) const {
    const std::size_t N = in_ordering.size();
    out_ordering.resize(N);
    for (std::size_t u=0; u<N; ++u) {
      out_ordering[in_ordering[u]] = u;
    }
  }
  
  void graph_minimum_linear_arranger::arrangement_in(std::vector<std::size_t>& best_ordering,
						     const sl::static_weighted_graph& graph) {
    sl::static_weighted_graph undirected_graph = graph;
    undirected_graph.make_undirected();
    std::cerr << "*** MINLA(" << undirected_graph.node_count() << " nodes/" <<  undirected_graph.node_count() << " edges)" << std::endl;

    float best_cost = 1e30f;
    for (std::size_t rk=0; rk<restart_k_; ++rk) {
      std::vector<std::size_t> ordering;

      std::cerr << std::endl;

      initialize_in(ordering, undirected_graph);
      float c= cost(ordering, graph); if ((rk==0) || (c<best_cost)) { best_cost = c; best_ordering=ordering; }
      std::cerr << rk+1 << "/" << restart_k_ << "|    [INIT] COST=" << c << " BEST=" << best_cost << std::endl;

      greedy_improve(ordering, undirected_graph);
      c= cost(ordering, graph); if (c<best_cost) { best_cost = c; best_ordering=ordering; }
      std::cerr << rk+1 << "/" << restart_k_ << "|    [GREEDY] COST=" << c << " BEST=" << best_cost << std::endl;
      
      relax(ordering, undirected_graph, relax_beta_, relax_k1_);
      c= cost(ordering, graph); if (c<best_cost) { best_cost = c; best_ordering=ordering; }
      std::cerr << rk+1 << "/" << restart_k_ << "|    [RELAX(" << relax_beta_ << ", " << relax_k1_ << ")] COST=" <<  c << " BEST=" << best_cost << std::endl;
    
      for (std::size_t k=0; k<relax_k2_; ++k) {
	relax(ordering, undirected_graph, relax_beta_, relax_k2_-k);
	c= cost(ordering, graph); if (c<best_cost) { best_cost = c; best_ordering=ordering; }
	std::cerr << rk+1 << "/" << restart_k_ << "|    " << k+1 << "/" << relax_k2_ << ": [RELAX(" << relax_beta_ << ", " << relax_k2_-k << ")] COST=" << c << " BEST=" << best_cost << std::endl;
	v_cycle(ordering, undirected_graph);
	c= cost(ordering, graph); if (c<best_cost) { best_cost = c; best_ordering=ordering; }
	std::cerr << rk+1 << "/" << restart_k_ << "|    " << k+1 << "/" << relax_k2_ << ": [V_CYCLE    ] COST=" << c << " BEST=" << best_cost << std::endl;
      }
    }

    // Fine-tune best ordering
    std::cerr << std::endl;

    std::vector<std::size_t> ordering = best_ordering;
    bool improved = true;
    while (improved) {
      improved = false;
      local_optimization_sweep(ordering, undirected_graph, localopt_k_);
      float c= cost(ordering, graph); if (c<best_cost) { improved=true; best_cost = c; best_ordering=ordering; }
      std::cerr << "FINAL FINE_TUNE COST=" << c << " BEST=" << best_cost << std::endl;
    }
  }
 
  void graph_minimum_linear_arranger::initialize_in(std::vector<std::size_t>&        ordering,
						    const sl::static_weighted_graph& graph) {
    switch (init_kind_) {
    case CANONICAL: canonical_initialize_in(ordering, graph); break;
    case DFS      : dfs_initialize_in(ordering, graph); break;
    case BFS      : bfs_initialize_in(ordering, graph); break;
    case GREEDY   : greedy_initialize_in(ordering, graph); break;
    default:
      std::cerr << "(EEE) Unknown init kind: " << (int)init_kind_ << std::endl;
      bfs_initialize_in(ordering, graph);
    }
  }

  void graph_minimum_linear_arranger::canonical_initialize_in(std::vector<std::size_t>&        ordering,
							      const sl::static_weighted_graph& graph) {
    const  std::size_t N = graph.node_count();
    
    // Trivial lexial ordering
    ordering.resize(N);
    for (std::size_t u=0; u<N; ++u) {
      ordering[u] = u;
    }
  }

  
  void graph_minimum_linear_arranger::dfs_initialize_in(std::vector<std::size_t>&        ordering,
							      const sl::static_weighted_graph& graph) {
    const  std::size_t N = graph.node_count();
    
    // Trivial lexial ordering
    ordering.resize(N);

    std::vector<std::size_t> dfs_stack;
    std::vector<bool> dfs_visited(N, false);

    std::vector<std::size_t> shuffled(N);
    for (std::size_t u=0; u<N; ++u) {
      shuffled[u] = u;
    }
    std::random_shuffle(shuffled.begin(), shuffled.end());
    
    std::size_t o=0;
    for (std::size_t root_n=0; root_n<N; ++root_n) {
      std::size_t root_u = shuffled[root_n];
      if (!dfs_visited[root_n]) {
	dfs_stack.push_back(root_u);
	while (!dfs_stack.empty()) {
	  std::size_t u = dfs_stack.back();
	  dfs_stack.pop_back();
	  if (!dfs_visited[u]) {
	    dfs_visited[u]=true;
	  
	    ordering[u]=o; ++o;

	    const std::size_t e_bgn = graph.node_adjacency(u);
	    const std::size_t e_end = e_bgn + graph.node_degree(u);
	    for (std::size_t e=e_bgn; e<e_end; ++e) {
	      const std::size_t v=graph.edge_target(e);
	      if (!dfs_visited[v]) {
		dfs_stack.push_back(v);
	      }
	    }
	  }
	}
      }
    }
  }
  
  void graph_minimum_linear_arranger::bfs_initialize_in(std::vector<std::size_t>&        ordering,
							const sl::static_weighted_graph& graph) {
    const  std::size_t N = graph.node_count();
    
    // Trivial lexial ordering
    ordering.resize(N);

    std::vector<bool> bfs_visited(N, false);

    std::vector<std::size_t> shuffled(N);
    for (std::size_t u=0; u<N; ++u) {
      shuffled[u] = u;
    }
    std::random_shuffle(shuffled.begin(), shuffled.end());
    
    std::size_t o=0;
    for (std::size_t root_n=0; root_n<N; ++root_n) {
      std::size_t u = shuffled[root_n];
      if (!bfs_visited[u]) {
	bfs_visited[u]=true;
	ordering[u]=o; ++o;
	
	const std::size_t e_bgn = graph.node_adjacency(u);
	const std::size_t e_end = e_bgn + graph.node_degree(u);
	for (std::size_t e=e_bgn; e<e_end; ++e) {
	  const std::size_t v=graph.edge_target(e);
	  if (!bfs_visited[v]) {
	    bfs_visited[v]=true;
	    ordering[v]=o; ++o;
	  }
	}
      }
    }
  }

  void graph_minimum_linear_arranger::greedy_initialize_in(std::vector<std::size_t>&        ordering,
							   const sl::static_weighted_graph& graph) {
    const  std::size_t N = graph.node_count();
   
    ordering.resize(N);

    std::vector<std::size_t> shuffled(N);
    std::vector<bool> labeled(N, false);
    std::set<std::size_t> unlabeled;
    for (std::size_t u=0; u<N; ++u) {
      shuffled[u] = u;
      unlabeled.insert(u);
    }
    std::random_shuffle(shuffled.begin(), shuffled.end());
      
    std::vector<int> layout(N);
    int l=-1; int r=0;
    for (std::size_t root_n=0; root_n<N; ++root_n) {
      std::size_t root_u = shuffled[root_n];
      if (!labeled[root_u]) {
	// New root -- as it has no connection with currently
	// labeled nodes, arbitrarily put it to the right
	labeled[root_u]=true;
	unlabeled.erase(root_u);
	layout[root_u]=r; ++r; 
	
	while (!unlabeled.empty()) {
	  std::size_t best_u = std::size_t(-1);
	  float best_sf = 0.0f; float best_cl = 0.0f; float best_cr = 0.0f;
	  for (std::set<std::size_t>::iterator unlabeled_it = unlabeled.begin();
	       unlabeled_it != unlabeled.end();
	       ++unlabeled_it) {
	    // FIXME we should actually loop only on unlabeled neighbors of labeled...
	    std::size_t u = *unlabeled_it;
	    float sf_u = 0.0f;
	    float cl_u = 0.0f;
	    float cr_u = 0.0f;
	    const std::size_t e_bgn = graph.node_adjacency(u);
	    const std::size_t e_end = e_bgn + graph.node_degree(u);
	    for (std::size_t e=e_bgn; e<e_end; ++e) {
	      const std::size_t v=graph.edge_target(e);
	      const float       w_uv = graph.edge_weight(e);
	      if (labeled[v]) {
		sf_u -= w_uv;

		cl_u+= w_uv * sl::abs(float(layout[v])-float(l));
		cr_u+= w_uv * sl::abs(float(layout[v])-float(r));
	      } else {
		sf_u += w_uv;
	      }
	    }
	    if ((best_u == std::size_t(-1)) || (sf_u<best_sf)) {
	      // First or improvement, accept it
	      best_u = u;
	      best_sf = sf_u; best_cl = cl_u; best_cr = cr_u;
	    } else if (sf_u==best_sf) {
	      // Same -- accept based on cost
	      if (std::min(cl_u,cr_u)<std::min(best_cr, best_cl)) {
		best_u = u;
		best_sf = sf_u; best_cl = cl_u; best_cr = cr_u;
	      }
	    }
	  }
	  labeled[best_u]=true;
	  unlabeled.erase(best_u);
	  if (best_cl<best_cr) {
	    layout[best_u]=l; --l;
	  } else {
	    layout[best_u]=r; ++r;
	  }
	} // for each unlabeled
      } // if root not already labeled
    } // for each root

    // Normalization into ordering
    for (std::size_t u=0; u<N; ++u) {
      ordering[u]=std::size_t(layout[u]-(l+1));
    }    
  }

  void graph_minimum_linear_arranger::greedy_improve(std::vector<std::size_t>&        ordering,
						     const sl::static_weighted_graph& graph) {
    const  std::size_t N = graph.node_count();

    std::vector<int> layout(N);
    int l=-1; int r=1;
    layout[ordering[0]]=0;
    for (std::size_t i=1; i<N; ++i) {
      std::size_t u=ordering[i];
      // Left 
      layout[u]=l; float cl=0;
      for (std::size_t j=0; j<i; ++j) {
	std::size_t v=ordering[j];
	if (graph.are_connected(u,v)) {	
	  cl+= graph.edge_weight(u,v) * sl::abs(float(layout[u])-float(layout[v]));
	}
      }
      // Right */
      layout[u]=r; float cr=0;
      for (std::size_t j=0; j<i; ++j) {
	std::size_t v= ordering[j];
	if (graph.are_connected(u,v)) {
	  cr+= graph.edge_weight(u,v) * sl::abs(float(layout[u])-float(layout[v]));
	}
      }
      // Decision 
      if (cl<cr) {
	layout[u]=l; --l;
      } else {
	++r;
      }
    }
    // Normalization into ordering
    for (std::size_t u=0; u<N; ++u) {
      ordering[u]=std::size_t(layout[u]-(l+1));
    }    
  }
  
  void graph_minimum_linear_arranger::coarsen_in(std::vector<std::size_t>&        coarse_ordering,
						 sl::static_weighted_graph&       coarse_graph,
						 std::vector<std::size_t>&        fine_to_coarse_vertex_index,
						 std::vector<std::pair<std::size_t,std::size_t> >& coarse_to_fine_vertex_index,
						 const std::vector<std::size_t>&  fine_ordering,
						 const sl::static_weighted_graph& fine_graph) {
    const std::size_t N_fine = fine_graph.node_count();

    std::size_t unmatched_idx = std::size_t(-1); // no unmatched
    if (N_fine%2) {
      // odd vertex count, one node should remain unmatched
      unmatched_idx = std::size_t(irnd_.value_within(0, N_fine-1));
      if (unmatched_idx%2) {
	// We want the unmatched node to be even
	unmatched_idx = (unmatched_idx+1)%N_fine;
      }
    }
    //std::cerr << "COARSEN: unmatched="<< unmatched_idx << std::endl;
    
    // Compute mapping
    fine_to_coarse_vertex_index.resize(N_fine);
    coarse_to_fine_vertex_index.clear();

    std::vector<std::size_t> sorted_fine_graph_nodes;
    reciprocal_in(sorted_fine_graph_nodes, fine_ordering);
    for (std::size_t i=0; i<N_fine; ++i) {
      std::size_t cidx = coarse_to_fine_vertex_index.size();
      std::size_t u0 = sorted_fine_graph_nodes[i];
      fine_to_coarse_vertex_index[u0] = cidx;
      if ((i!=unmatched_idx) && ((i+1)<N_fine)) {
	++i;
	std::size_t u1 = sorted_fine_graph_nodes[i];
	fine_to_coarse_vertex_index[u1] = cidx;
	coarse_to_fine_vertex_index.push_back(std::make_pair(u0,u1));
	//std::cerr << "(" << cidx << ": " << u0 << " " << u1 << ")";
      } else {
	coarse_to_fine_vertex_index.push_back(std::make_pair(u0,u0));
	//std::cerr << "(" << cidx << ": " << u0 << ")";
      }
    }
    std::size_t N_coarse = coarse_to_fine_vertex_index.size();
    //std::cerr << std::endl;
    
    assert(DBG_CHECK_UNDIRECTED("fine_graph", fine_graph));
        
    // Builde edge map using desired mapping
    std::size_t M_coarse = 0;
    std::vector< std::map<std::size_t, float> > coarse_edge_weight_map(N_coarse);
    for (std::size_t u_fine=0; u_fine<N_fine; ++u_fine) {
      const std::size_t u_coarse = fine_to_coarse_vertex_index[u_fine];
      
      const std::size_t e_fine_bgn = fine_graph.node_adjacency(u_fine);
      const std::size_t e_fine_end = e_fine_bgn + fine_graph.node_degree(u_fine);
      for (std::size_t e_fine=e_fine_bgn; e_fine<e_fine_end; ++e_fine) {
	const std::size_t v_fine=fine_graph.edge_target(e_fine);
	const std::size_t v_coarse = fine_to_coarse_vertex_index[v_fine];
	const float w_uv = fine_graph.edge_weight(e_fine);
	
	if (coarse_edge_weight_map[u_coarse].find(v_coarse) == coarse_edge_weight_map[u_coarse].end()) {
	  coarse_edge_weight_map[u_coarse][v_coarse] = w_uv;
	  ++M_coarse;
	} else {
	  coarse_edge_weight_map[u_coarse][v_coarse] += w_uv;
	}
      }
    }

    // Construct coarse graph
    
    coarse_ordering.resize(N_coarse);

    std::vector<std::size_t>  coarse_graph_node_degree(N_coarse);
    std::vector<std::size_t>  coarse_graph_edge_target(M_coarse);
    std::vector<float>        coarse_graph_edge_weight(M_coarse);

    std::size_t  coarse_edge_offset=0;
    for (std::size_t i=0; i<N_coarse; ++i) {
      coarse_ordering[i] = i;
      const std::size_t M_i = coarse_edge_weight_map[i].size();
      coarse_graph_node_degree[i] = M_i;
      std::size_t k=0;
      for (std::map<std::size_t, float>::iterator it = coarse_edge_weight_map[i].begin();
	   it != coarse_edge_weight_map[i].end();
	   ++it) {
	std::size_t e = coarse_edge_offset + k;
	coarse_graph_edge_target[e] = it->first;
 	coarse_graph_edge_weight[e] = it->second;
	++k;
      }
      coarse_edge_offset += M_i;
    }

    coarse_graph.init(coarse_graph_node_degree,
		      coarse_graph_edge_target,
		      coarse_graph_edge_weight);

    assert(DBG_CHECK_UNDIRECTED("coarse_graph", coarse_graph));
  }

  void graph_minimum_linear_arranger::uncoarsen_in(std::vector<std::size_t>&        fine_ordering,
						   const sl::static_weighted_graph& fine_graph,
						   const std::vector<std::size_t>&  coarse_ordering,
						   const std::vector<std::size_t>&             fine_to_coarse_vertex_index,
						   const std::vector<std::pair<std::size_t,std::size_t> >& coarse_to_fine_vertex_index) {
    assert(coarse_ordering.size()==coarse_to_fine_vertex_index.size());
    assert(fine_ordering.size()>=coarse_ordering.size());

    const std::size_t N_coarse = coarse_ordering.size();

    
    //std::cerr << "UNCOARSEN" << std::endl;

    std::size_t fidx = 0;
    for (std::size_t cidx=0; cidx<N_coarse; ++cidx) {
      assert(fidx<N_coarse);
      const std::size_t u_coarse = coarse_ordering[cidx];
      const std::size_t u0_fine = coarse_to_fine_vertex_index[u_coarse].first;
      const std::size_t u1_fine = coarse_to_fine_vertex_index[u_coarse].second;
      
      //std::cerr << "(" << fidx << ":" << u0_fine << ")";
      fine_ordering[u0_fine] = fidx; ++fidx;
      if (u1_fine!=u0_fine)  {
	assert(fidx<fine_ordering.size());
	//std::cerr << "(" << fidx << ":" << u1_fine << ")";
	fine_ordering[u1_fine] = fidx; ++fidx;
	
      }
    }
    assert(fidx == fine_ordering.size());
    //std::cerr  << std::endl;

#if 0
    // As there is an ambiguity here for the placement of the fine nodes contained in the same coarse node,
    // we do a quick optimization by evaluating the cost of swapping them.
    // Look at Safro weighted MinLA for a better approach
    const std::size_t N_fine   = fine_ordering.size();

    std::vector<std::size_t> sorted_fine_graph_nodes;
    reciprocal_in(sorted_fine_graph_nodes, fine_ordering);
    for (std::size_t i_fine_bgn=0; i_fine_bgn<N_fine; ) {
      std::size_t i_fine_end=i_fine_bgn;
      while ((i_fine_end<N_fine) && (fine_to_coarse_vertex_index[i_fine_end] == fine_to_coarse_vertex_index[i_fine_bgn])) {
	++i_fine_end;
      }
      //std::cerr << "Optimizing " << i_fine_bgn << ".." << i_fine_end << std::endl;

      optimize_small_local_sequence(fine_ordering,
				    sorted_fine_graph_nodes,
				    fine_graph, 
				    i_fine_bgn, i_fine_end-i_fine_bgn);
      i_fine_bgn = i_fine_end;
    }
    //std::cerr << "v_cycle: COST AFTER FINE TWEAKING: " << N_fine << " nodes -- cost = " << cost(fine_ordering, fine_graph) << std::endl;
#endif
  }

  float graph_minimum_linear_arranger::relaxed_l2(std::size_t u,
						       const std::vector< std::pair<float,std::size_t> >& placement,
						       const sl::static_weighted_graph&  graph) const {
    const std::size_t M_u = graph.node_degree(u);
    if (M_u==0) {
      return placement[u].first;
    } else {
      // PLAIN AVERAGE
      float p_sum = 0.0f;
      float w_sum = 0.0f;
      const std::size_t e_bgn = graph.node_adjacency(u);
      const std::size_t e_end = e_bgn + graph.node_degree(u);
      for (std::size_t e=e_bgn; e<e_end; ++e) {
	const std::size_t v=graph.edge_target(e);
	const float w_uv = graph.edge_weight(e);
	p_sum += w_uv * placement[v].first;
	w_sum += w_uv;
      }
      return (p_sum / w_sum); // weighted average of neighbors
    }
  }
  
  float graph_minimum_linear_arranger::relaxed_l1_irwls(std::size_t u,
							const std::vector< std::pair<float,std::size_t> >& placement,
							const sl::static_weighted_graph&  graph) const {
    const std::size_t M_u = graph.node_degree(u);
    float placement_u = placement[u].first;

    if (M_u==0) {
      return placement_u;
    } else {
      // L1 norm solution approximated with irwls
      const std::size_t irwls_k   = 4;
      const float       irwls_eps = 1e-6f;
      const std::size_t e_bgn = graph.node_adjacency(u);
      const std::size_t e_cnt = graph.node_degree(u);
      const std::size_t e_end = e_bgn + e_cnt;
      std::vector<float> irwls_w(e_cnt, 1.0);
      for (std::size_t k=0; k<irwls_k; ++k) {
	float p_sum = 0.0f;
	float w_sum = 0.0f;
	for (std::size_t e=e_bgn; e<e_end; ++e) {
	  const std::size_t v=graph.edge_target(e);
	  const float w_uv = graph.edge_weight(e) * irwls_w[e-e_bgn];
	  p_sum += w_uv * placement[v].first;
	  w_sum += w_uv;
	}
	placement_u = (p_sum / w_sum);
	// Update weights from residual
	for (std::size_t e=e_bgn; e<e_end; ++e) {
	  const std::size_t v=graph.edge_target(e);
	  const float w_uv = graph.edge_weight(e);
	  irwls_w[e-e_bgn] = 1.0f/std::max(irwls_eps, w_uv * sl::abs(placement_u-placement[v].first));
	}
      }
      
      return placement_u;
    }
  }
						    
  float graph_minimum_linear_arranger::relaxed_l1_median(std::size_t u,
						      const std::vector< std::pair<float,std::size_t> >& placement,
						      const sl::static_weighted_graph&  graph) const {
    // FIXME not working... 
    const std::size_t M_u = graph.node_degree(u);
    if (M_u<3) {
      return relaxed_l2(u, placement, graph);
    } else {
      std::vector< std::pair<float,float> > p_w; p_w.reserve(M_u);
      double w_sum = 0.0f;
      const std::size_t e_bgn = graph.node_adjacency(u);
      const std::size_t e_end = e_bgn + graph.node_degree(u);
      for (std::size_t e=e_bgn; e<e_end; ++e) {
	const std::size_t v=graph.edge_target(e);
	const float w_uv = graph.edge_weight(e);
	p_w.push_back(std::make_pair(placement[v].first, w_uv));
	w_sum += w_uv;
      }
      std::sort(p_w.begin(), p_w.end()); // Sort by placement
      
      // Extract weighted median
      const double w_mid = 0.5*w_sum;
      std::size_t idx_geq = 0;
      double w_geq = p_w[idx_geq].second;
      while ((w_geq < w_mid) && (idx_geq+1<M_u))  {
	++idx_geq; w_geq += p_w[idx_geq].second;
      }
      if (idx_geq==0) {
	// Empty set
	return (p_w[0].first);
      } else if (w_geq == w_mid && (idx_geq+1<M_u)) {
	std::size_t idx_gt = idx_geq+1;
	return 0.5f * (p_w[idx_geq].first + p_w[idx_gt].first);
      } else {
	return (p_w[idx_geq].first);
      }
    } // if median applicable
  }
  
  void graph_minimum_linear_arranger::relax(std::vector<std::size_t>&         ordering,
					    const sl::static_weighted_graph&  graph,
					    float beta,
					    std::size_t sweep_count) {
    const std::size_t N = graph.node_count();

    // Iteratively place each node at the barycenter of the neighbors
    std::vector< std::pair<float,std::size_t> > placement(N);
    for (std::size_t i=0; i<N; ++i) {
      placement[i] = std::make_pair(float(ordering[i]), i);
    }
    for (std::size_t s=0; s<sweep_count; ++s) {
      std::vector<std::size_t> shuffled(N);
      for (std::size_t u=0; u<N; ++u) {
	shuffled[u] = u;
      }
      std::random_shuffle(shuffled.begin(), shuffled.end());
      for (std::size_t uu=0; uu<N; ++uu) {
	const std::size_t u = shuffled[uu];
#if 0
	const float relaxed_placement_u = relaxed_l2(u, placement, graph);
	
#elif 0
	const float relaxed_placement_u = relaxed_l1_median(u, placement, graph);
#else
	const float relaxed_placement_u = relaxed_l1_irwls(u, placement, graph);
#endif

	placement[u].first = (1.0f-beta) * placement[u].first + beta * relaxed_placement_u;
      } // foreach node u
    } // foreach sweep s

    // Convert placement to linear ordering
    std::sort(placement.begin(), placement.end());
    for (std::size_t i=0; i<N; ++i) {
      ordering[placement[i].second] = i;
    }
  }
  
  void graph_minimum_linear_arranger::v_cycle(std::vector<std::size_t>&        ordering,
					      const sl::static_weighted_graph& graph) {
    const std::size_t N = graph.node_count();

    if (N<=localopt_n_) {
      std::vector<std::size_t> sorted_nodes;
      reciprocal_in(sorted_nodes, ordering);
      optimize_small_local_sequence(ordering, sorted_nodes, graph, 
				    0, N);
    } else {
      std::vector<std::size_t>             coarse_ordering;
      std::vector<std::size_t>             fine_to_coarse_vertex_index;
      std::vector<std::pair<std::size_t,std::size_t> > coarse_to_fine_vertex_index;
      sl::static_weighted_graph            coarse_graph;

      coarsen_in(coarse_ordering, coarse_graph,
		 fine_to_coarse_vertex_index, coarse_to_fine_vertex_index,
		 ordering, graph);
  
      v_cycle(coarse_ordering, coarse_graph);

      uncoarsen_in(ordering,
		   graph,
		   coarse_ordering,
		   fine_to_coarse_vertex_index, coarse_to_fine_vertex_index);
      
      local_optimization_sweep(ordering, graph, localopt_k_);
    }
  }

  void graph_minimum_linear_arranger::local_optimization_sweep(std::vector<std::size_t>&        ordering,
							       const sl::static_weighted_graph& graph,
							       std::size_t sweep_count) {
    const std::size_t N = graph.node_count();
    std::vector<std::size_t> sorted_nodes;
    reciprocal_in(sorted_nodes, ordering);
    for (std::size_t k=0; k<sweep_count; ++k) {
      std::vector<std::size_t> shuffled(N);
      for (std::size_t j=0; j<N; ++j) {
	shuffled[j] = j;
      }
      std::random_shuffle(shuffled.begin(), shuffled.end());
      for (std::size_t j=0; j<N; ++j) {
	optimize_small_local_sequence(ordering, sorted_nodes, graph, 
				      shuffled[j], std::min(localopt_n_,N-shuffled[j]));
      }
    }
  }
  
								    
  void graph_minimum_linear_arranger::optimize_small_local_sequence(std::vector<std::size_t>&        ordering,
								    std::vector<std::size_t>&        sorted_nodes, // reciprocal of ordering
								    const sl::static_weighted_graph& graph,
								    std::size_t seq_start, std::size_t seq_sz) {
    if (seq_sz>1) {
#ifndef NDEBUG
      const float pre_cost = cost(ordering, graph);
#endif 

      // Try best combination of subsequence
	
      float best_cost = 0.0f;
      std::vector<std::size_t> best_subsequence(seq_sz);
      std::copy(sorted_nodes.begin()+seq_start, sorted_nodes.begin()+seq_start+seq_sz, best_subsequence.begin());

#if 0
      DBG_DUMP("ORIG", sorted_nodes, seq_start, seq_sz); std::cerr << " => C=" << pre_cost << std::endl;
      DBG_DUMP("ORIG/FULL", sorted_nodes, 0, sorted_nodes.size()); std::cerr << std::endl;
#endif
      
      sl::star_swap_permutator perm(seq_sz);
      while (!perm.off()) {
	// Find next permutation of subsequence and update ordering
	std::swap(sorted_nodes[seq_start], sorted_nodes[seq_start+perm.value()]);
	ordering[sorted_nodes[seq_start]]              = seq_start;
	ordering[sorted_nodes[seq_start+perm.value()]] = seq_start+perm.value();

	// =================================
	// Find cost of this permutation
	// NOTE: We could incrementally compute costs by just recomputing the cost of the swap.
	// We don't do so for simplicity (and we assume that seq_sz is small...)
	float local_cost = 0.0f;
	for (std::size_t s=seq_start; s<seq_start+seq_sz; ++s) {
	  const std::size_t u = sorted_nodes[s];
	  const std::size_t e_bgn = graph.node_adjacency(u);
	  const std::size_t e_end = e_bgn + graph.node_degree(u);
	  for (std::size_t e=e_bgn; e<e_end; ++e) {
	    const std::size_t v=graph.edge_target(e);
	    float w_uv = graph.edge_weight(e);
	    if (ordering[v]<seq_start || ordering[v]>=seq_start+seq_sz) w_uv *= 2.0f; // Count backward link
	    local_cost += w_uv*sl::abs(float(ordering[u])-float(ordering[v]));
	  }
	}
#if 0
	DBG_DUMP("PERM", sorted_nodes, seq_start, seq_sz); std::cerr << " => C=" << cost(ordering, graph) << " LC=" << local_cost << std::endl;
#endif

	// Update best permutation
	if ((perm.value() == 0) || (local_cost<best_cost)) {
	  std::copy(sorted_nodes.begin()+seq_start, sorted_nodes.begin()+seq_start+seq_sz, best_subsequence.begin());
	  best_cost = local_cost;
	}

	// Compute next swap
	++perm;
      }

      // Update best subsequence and ordering
      std::copy(best_subsequence.begin(), best_subsequence.end(), sorted_nodes.begin()+seq_start);
      for (std::size_t u=seq_start; u<seq_start+seq_sz; ++u) {
	ordering[sorted_nodes[u]] = u;
      }
#if 0
      DBG_DUMP("BEST", sorted_nodes, seq_start, seq_sz); std::cerr << " => C=" << cost(ordering, graph) << " LC=" << best_cost << std::endl;
      DBG_DUMP("BEST/FULL", sorted_nodes, 0, sorted_nodes.size()); std::cerr << std::endl;
#endif
#ifndef NDEBUG
      float post_cost = cost(ordering, graph);
      if (pre_cost<post_cost) std::cerr << "!!! BAD LOCAL OPT: " << graph.node_count() << " nodes -- cost = " << pre_cost << "->" << post_cost << std::endl;
#endif
    }
  } 

} // namespace sl

