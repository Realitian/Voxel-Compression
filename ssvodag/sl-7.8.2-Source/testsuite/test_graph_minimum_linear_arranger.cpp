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

static std::size_t failed_test_count = 0;


void test_simple() {
  sl::tester tester("graph_minimum_linear_arranger");
  
  // Simple graph from Jordi Petit testsuite
  std::string simple_gra =
    "5\n"
    "8\n"
    "3 3 3 3 4\n"
    "1 3 4 0 2 4 1 3 4 0 2 4 0 1 2 3 -1\n";
  
  sl::static_weighted_graph graph;
  std::istringstream iss(simple_gra);
  graph.load_gra(iss);

  tester.test("Node count", graph.node_count(), std::size_t(5));
  tester.test("Edge count", graph.edge_count(), std::size_t(16));
  tester.test("Component count", graph.component_count(), std::size_t(1));
  
  // Trivial ordering
  std::vector<std::size_t> ordering(graph.node_count());
  for (std::size_t i=0; i<graph.node_count(); ++i) ordering[i] = i;

  // MINLA
  sl::graph_minimum_linear_arranger minla;

  float c0 = minla.cost(ordering, graph);
  tester.test("Initial cost", c0, sl::intervalf("16.00"));

  std::vector<std::size_t> minla_ordering;
  minla.arrangement_in(minla_ordering, graph);
  
  float c1 = minla.cost(minla_ordering, graph);
  tester.test("final cost", c1, sl::intervalf("14.00"));
  
  failed_test_count += tester.failed_test_count();
}

int main() {
  test_simple();
  
  return failed_test_count;
}


