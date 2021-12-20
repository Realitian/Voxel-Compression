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
#include <sl/triangle_mesh_stripifier.hpp>
#include <sl/triangle_mesh_vertex_clusterer.hpp>
#include <sl/triangle_mesh_edge_collapser.hpp>
#include <sl/triangle_mesh_refiner.hpp>
#include <sl/triangle_mesh_interior_culler.hpp>
#include <sl/triangle_mesh_distance_sampler.hpp>
#include <sl/triangle_mesh_uniform_remesher.hpp>
#include <sl/clock.hpp>
#include <algorithm>
#include <cassert>

static std::size_t failed_test_count = 0;

class vdata_t {
public:
  typedef float       value_t;
  typedef sl::point3f point_t;
protected:
  point_t position_;

public:
  vdata_t() {}
  
  public: // Triangle Mesh interface
    inline point_t&  position()       { return position_; }
    inline const point_t&  position() const { return position_; }

    inline vdata_t lerp(vdata_t& other, float t) const {
      vdata_t result = vdata_t();
      result.position() = position().lerp(other.position(), t);
      return result;
    }
};

struct edata_t {};
struct tdata_t {};
    
typedef sl::triangle_mesh<vdata_t, edata_t, tdata_t> tmesh_t;

tmesh_t* new_test_mesh(std::size_t tricount,
                       bool regular) {
  typedef tmesh_t::vertex_data_t::point_t point_t;
  typedef point_t::value_t value_t;
  tmesh_t* result = new tmesh_t;
  std::vector<tmesh_t::edge_t> edges;

  std::size_t M = sl::isqrt(tricount);
  std::size_t N = tricount/M;
  if (regular) {
    N/=2;
  }
  M = M+1;
  N = N+1;
  
  for (std::size_t i=0; i<M-1; ++i) {
    for (std::size_t j=0; j<N-1; ++j) {
      std::size_t i00 = (i+0)+(j+0)*M;
      std::size_t i01 = (i+0)+(j+1)*M;
      std::size_t i10 = (i+1)+(j+0)*M;
      std::size_t i11 = (i+1)+(j+1)*M;
      result->insert_triangle(tmesh_t::triangle_t(i00,i10,i11));
      result->insert_triangle(tmesh_t::triangle_t(i00,i11,i01));
      result->vertex_data(i00)->position() = point_t(value_t(i+0),value_t(j+0),value_t(0));
      result->vertex_data(i01)->position() = point_t(value_t(i+0),value_t(j+1),value_t(0));
      result->vertex_data(i10)->position() = point_t(value_t(i+1),value_t(j+1),value_t(0));
      result->vertex_data(i11)->position() = point_t(value_t(i+1),value_t(j+1),value_t(1));
      edges.push_back(tmesh_t::edge_t(i00,i10));
      edges.push_back(tmesh_t::edge_t(i10,i11));
      edges.push_back(tmesh_t::edge_t(i11,i00));
    }
  }

  if (!regular) {
    std::random_shuffle(edges.begin(),edges.end());
    while ((result->triangle_count()>tricount) &&
           !edges.empty()) {
      // Randomly collapse edge
      tmesh_t::edge_t e = edges[edges.size()-1];
      edges.pop_back();
      if (result->has_edge(e)) {
        result->collapse_edge(e[0],e[1]);
      }
    }
  }  
  return result;
}

std::size_t strip_cache_miss_count(std::size_t cache_size,
                                   const std::vector<std::size_t>& stitched_strip) {
  sl::fifo_cache_simulator<std::size_t> cache;
  cache.resize(std::max(std::size_t(2), cache_size));
  for (std::size_t i=0; i<stitched_strip.size(); ++i) {
    cache.insert(stitched_strip[i]);
  }
  return cache.miss_count();
}


std::size_t strip_extra_count(const tmesh_t* m,
                              const std::vector<std::size_t>& stitched_strip) {
  std::size_t extra_count = 0;
  std::set<tmesh_t::triangle_t> strip_tris;
  for (std::size_t i=2; i<stitched_strip.size(); ++i) {
    std::size_t v0= stitched_strip[i-2];
    std::size_t v1= stitched_strip[i-1];
    std::size_t v2= stitched_strip[i-0];
    tmesh_t::triangle_t t(v0, (i%2==0)?v1:v2, (i%2==0)?v2:v1);
    if (t.is_valid()) {
      strip_tris.insert(t);
      if (!m->has_triangle(t)) {
        ++extra_count;
        std::cerr << "Strip contains extra tri: " << t[0] << " " << t[1] << " " << t[2] << std::endl;
      }
    }
  }
  return extra_count;
}

std::size_t strip_missed_count(const tmesh_t* m,
                              const std::vector<std::size_t>& stitched_strip) {
  std::set<tmesh_t::triangle_t> strip_tris;
  for (std::size_t i=2; i<stitched_strip.size(); ++i) {
    std::size_t v0= stitched_strip[i-2];
    std::size_t v1= stitched_strip[i-1];
    std::size_t v2= stitched_strip[i-0];
    tmesh_t::triangle_t t(v0, (i%2==0)?v1:v2, (i%2==0)?v2:v1);
    if (t.is_valid()) {
      strip_tris.insert(t);
    }
  }
  std::size_t missed_count = 0;
  for (tmesh_t::triangle_map_t::const_iterator it = m->triangle_map().begin();
       it != m->triangle_map().end();
       ++it) {
    const tmesh_t::triangle_t& t = it->first;
    if (strip_tris.find(t) == strip_tris.end()) {
      ++missed_count;
      std::cerr << "Strip missed mesh tri: " << t[0] << " " << t[1] << " " << t[2] << std::endl;
    }
  }
  return missed_count;
}

std::size_t strip_double_count(const tmesh_t* m,
                               const std::vector<std::size_t>& stitched_strip) {
  std::set<tmesh_t::triangle_t> strip_tris;
  for (std::size_t i=2; i<stitched_strip.size(); ++i) {
    std::size_t v0= stitched_strip[i-2];
    std::size_t v1= stitched_strip[i-1];
    std::size_t v2= stitched_strip[i-0];
    tmesh_t::triangle_t t(v0, (i%2==0)?v1:v2, (i%2==0)?v2:v1);
    if (t.is_valid()) {
      strip_tris.insert(t);
    }
  }
  std::size_t double_count = strip_tris.size() > m->triangle_count() ? (strip_tris.size() - m->triangle_count()) : 0;
  return double_count;
}

void test_stripification(const std::string& test_name,
			 sl::triangle_mesh_stripifier* stripifier,
			 const tmesh_t* m) {
  sl::tester tester(test_name);
  
  sl::real_time_clock ck;
  ck.restart();

  stripifier->begin_input();
  for (tmesh_t::triangle_map_t::const_iterator it = m->triangle_map().begin();
       it!= m->triangle_map().end();
       ++it) {
    const tmesh_t::triangle_t& t = it->first;
    //std::cerr << "M: " << t[0] << " " << t[1] << " " << t[2] << std::endl;
    stripifier->insert_input_triangle(t[0], t[1], t[2]);
  }
  stripifier->end_input();
  
  std::vector<std::size_t> stitched_strip;
  std::size_t strip_count = 0;
  stripifier->begin_output();
  while (stripifier->has_output_strip()) {
    stripifier->begin_output_strip();
    ++strip_count;
    std::size_t v0 = stripifier->get_output_strip_vertex();
    if (stitched_strip.size() > 0) {
      stitched_strip.push_back(stitched_strip[stitched_strip.size()-1]);
      stitched_strip.push_back(v0);
      if (stitched_strip.size()%2==1) stitched_strip.push_back(v0);
    }
    assert(stitched_strip.size() % 2 == 0);
    stitched_strip.push_back(v0);
    while (stripifier->has_output_vertex()) {
      stitched_strip.push_back(stripifier->get_output_strip_vertex());
    }
    stripifier->end_output_strip();
  }
  stripifier->end_output();

  sl::time_duration strip_time = ck.elapsed();
  SL_USEVAR(strip_time);
  
  tester.test("Missed triangle count", strip_missed_count(m, stitched_strip), std::size_t(0));
  tester.test("Double triangle count", strip_double_count(m, stitched_strip), std::size_t(0));
  tester.test("Extra triangle count",  strip_extra_count(m, stitched_strip), std::size_t(0));
  tester.test("Strip efficiency", double(stitched_strip.size())/double(m->triangle_count()),
              sl::intervald(0.0, 3.0));
  if (stripifier->cache_size()>2) {
    tester.test("Cache efficiency", double(strip_cache_miss_count(stripifier->cache_size(), stitched_strip))/double(m->triangle_count()),
                sl::intervald(0.0, 1.0));
  }
	
  std::cerr << "Efficiency      : " << double(stitched_strip.size())/double(m->triangle_count()) << std::endl;
  std::cerr << "Cache efficiency: " <<double(strip_cache_miss_count(std::max(stripifier->cache_size(),std::size_t(2)), stitched_strip))/double(m->triangle_count()) << std::endl;
	    
  failed_test_count += tester.failed_test_count();
}
  
			 
			 
////////////////////////////////////////////////////////////////////////

void test_stripification(const std::string& test_name,
                         const tmesh_t* m,
                         std::size_t cache_size) {
  sl::triangle_mesh_greedy_stripifier gstripifier; // FIXME
  gstripifier.set_cache_size(cache_size);

  test_stripification("triangle_mesh_greedy_stripifier: " + test_name + " -- cache=" + sl::to_string(cache_size),
		      &gstripifier,
		      m);

  sl::triangle_mesh_backtracking_stripifier bstripifier; // FIXME
  bstripifier.set_cache_size(cache_size);

  test_stripification("triangle_mesh_backtracking_stripifier: " + test_name + " -- cache=" + sl::to_string(cache_size),
		      &bstripifier,
		      m);
}

void test_vertex_clusterer() {
  tmesh_t* m = new_test_mesh(2000, true);
  sl::tester tester("triangle_mesh_vertex_clusterer");
  std::size_t old_vertex_count = m->vertex_count();
  std::size_t old_boundary_vertex_count = m->boundary_vertex_count();
  std::size_t old_non_boundary_vertex_count = old_vertex_count-old_boundary_vertex_count;
  
  sl::triangle_mesh_cluster_non_boundary_vertices(*m, 2.0f);

  std::size_t new_vertex_count = m->vertex_count();
  std::size_t new_boundary_vertex_count = m->boundary_vertex_count();
  std::size_t new_non_boundary_vertex_count = new_vertex_count-new_boundary_vertex_count;
  
  tester.test("Vertex count", old_vertex_count>new_vertex_count, true);
  tester.test("Boundary vertex count", old_boundary_vertex_count, new_boundary_vertex_count);
  tester.test("Interior vertex count", old_non_boundary_vertex_count>new_non_boundary_vertex_count, true);
  
  delete m;
}

void test_edge_collapser() {
  tmesh_t* m = new_test_mesh(2000, true);
  sl::tester tester("triangle_mesh_edge_collapser");
  std::size_t old_vertex_count = m->vertex_count();
  std::size_t old_boundary_vertex_count = m->boundary_vertex_count();
  std::size_t old_non_boundary_vertex_count = old_vertex_count-old_boundary_vertex_count;
  
  sl::triangle_mesh_collapse_non_boundary_edges_to_target_edge_length(*m, 2.0f);

  std::size_t new_vertex_count = m->vertex_count();
  std::size_t new_boundary_vertex_count = m->boundary_vertex_count();
  std::size_t new_non_boundary_vertex_count = new_vertex_count-new_boundary_vertex_count;
  
  tester.test("Vertex count", old_vertex_count>new_vertex_count, true);
  tester.test("Boundary vertex count", old_boundary_vertex_count, new_boundary_vertex_count);
  tester.test("Interior vertex count", old_non_boundary_vertex_count>new_non_boundary_vertex_count, true);
  
  delete m;
}

void test_refiner() {
  tmesh_t* m = new_test_mesh(200, true);
  sl::tester tester("triangle_mesh_refiner");
  std::size_t old_vertex_count = m->vertex_count();
  std::size_t old_boundary_vertex_count = m->boundary_vertex_count();
  std::size_t old_non_boundary_vertex_count = old_vertex_count-old_boundary_vertex_count;
  
  sl::triangle_mesh_refine_non_boundary_to_target_tricount(*m, std::size_t(500));

  std::size_t new_vertex_count = m->vertex_count();
  std::size_t new_boundary_vertex_count = m->boundary_vertex_count();
  std::size_t new_non_boundary_vertex_count = new_vertex_count-new_boundary_vertex_count;
  
  tester.test("Vertex count", old_vertex_count<new_vertex_count, true);
  tester.test("Boundary vertex count", old_boundary_vertex_count, new_boundary_vertex_count);
  tester.test("Interior vertex count", old_non_boundary_vertex_count<new_non_boundary_vertex_count, true);
  
  delete m;
}

void test_uniform_remesher() {
  tmesh_t* m = new_test_mesh(200, true);
  sl::tester tester("triangle_mesh_uniform_remesher");
  
  sl::triangle_mesh_remesh_to_target_edge_length(*m, 0.5f);
  tester.test("Average edge length", sl::triangle_mesh_average_edge_length(*m), sl::intervalf(0.25f,0.75f));
  // FIXME tester.test("Minimum edge length", sl::triangle_mesh_minimum_edge_length(*m), sl::intervalf("0.5"));
  // FIXME tester.test("Maximum edge length", sl::triangle_mesh_maximum_edge_length(*m), sl::intervalf("0.5"));
  
  delete m;
}

void test_interior_culler() {
  tmesh_t* m = new_test_mesh(2000, true);
  sl::tester tester("triangle_mesh_interior_culler");
  std::size_t old_vertex_count = m->vertex_count();
  
  // FIXME: For now, only a stupid test on a fully visible mesh
  sl::triangle_mesh_cull_interior(*m);
  std::size_t new_vertex_count = m->vertex_count();
  tester.test("Vertex count", old_vertex_count, new_vertex_count);
  
  delete m;
}

void test_distance_sampler() {
  tmesh_t* m1 = new_test_mesh(2000, true);
  tmesh_t* m2 = new_test_mesh(2000, true); /*FIXME*/ m2->vertex_data(0)->position() += sl::vector3f(0.0f, 0.0f, 1.0f);
  sl::tester tester("triangle_mesh_distance_sampler");
  tester.test("Self Hausdorff distance", sl::triangle_mesh_symmetric_hausdorff_distance(*m1, *m1), sl::intervalf("0.000"));  
  tester.test("Self Hausdorff distance", sl::triangle_mesh_symmetric_hausdorff_distance(*m2, *m2), sl::intervalf("0.000"));  
  tester.test("Hausdorff distance", sl::triangle_mesh_symmetric_hausdorff_distance(*m1, *m2), sl::intervalf("1.0"));  
  delete m1;
  delete m2;
}

int main() { 
  tmesh_t* m_regular   = new_test_mesh(2000, true);
  tmesh_t* m_arbitrary = new_test_mesh(2000, false);
  
  test_stripification("regular grid", m_regular, 0);
  test_stripification("regular grid", m_regular, 32);

  test_stripification("arbitrary mesh", m_arbitrary, 0);
  test_stripification("arbitrary mesh", m_arbitrary, 32);

  test_vertex_clusterer();
  test_edge_collapser();
  test_refiner();
  test_uniform_remesher();
  test_interior_culler();
  test_distance_sampler();
  
  return failed_test_count;
}




