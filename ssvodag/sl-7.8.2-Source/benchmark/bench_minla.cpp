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
#include <sl/tester.hpp>
#include <sl/graph_minimum_linear_arranger.hpp>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cmath>

int main(int argc, const char* argv[]) {
  if (argc!=2) {
    std::cerr << "Usage: " << argv[0] << " <file.gra>" << std::endl;
    return 1;
  }

  std::string filename = argv[1];

  sl::static_weighted_graph graph;

  graph.load_gra(filename);
  if (graph.node_count()==0) {
    std::cerr<< "Graph is empty: aborting." << std::endl;
  } else {

    // Trivial ordering
    std::vector<std::size_t> ordering(graph.node_count());
    for (std::size_t i=0; i<graph.node_count(); ++i) ordering[i] = i;

    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "MINLA OF GRAPH WITH " << graph.node_count() << " NODES AND " << graph.edge_count() << " EDGES"<< std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "Optimizing..." << std::endl;
    
    sl::graph_minimum_linear_arranger minla;

    float c0 = minla.cost(ordering, graph);
    
    std::vector<std::size_t> minla_ordering;
    minla.arrangement_in(minla_ordering, graph);
  
    float c1 = minla.cost(minla_ordering, graph);
    
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "Initial cost: " << c0 << std::endl;
    std::cout << "Final cost: " << c1 << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
  }

  return 0;
}
