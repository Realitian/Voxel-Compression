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
/////// ALWAYS TEST IN DEBUG MODE
#if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////
#include <sl/fixed_size_point.hpp>
#include <sl/fixed_size_vector.hpp>
#include <vector>
#include <sl/kdtree.hpp>
#include <sl/tester.hpp>
#include <stdlib.h> // for rand()

static std::size_t failed_test_count = 0;

static inline void random_seed(unsigned int seed) {
  srand(seed);
}

static inline sl::point3f random_point() {
  float x = (rand() / (float)(RAND_MAX));
  float y = (rand() / (float)(RAND_MAX));
  float z = (rand() / (float)(RAND_MAX));
  return sl::point3f(x,y,z);
}

/**
 *  A simple particle
 */
class particle {
protected:
  std::size_t id_;
  sl::point3f position_;
  sl::vector3f velocity_;
  float mass_;
public:

  bool operator==(const particle& other) const {
    return id_ == other.id_;
  }

  inline particle() {
    id_ = 0;
    mass_=1.0f;   
  }
  
  inline particle(std::size_t id,
		  const sl::point3f& pos, 
		  const sl::vector3f& vel, 
		  float m) 
    :
    id_(id), position_(pos), velocity_(vel), mass_(m) {
  }

  inline std::size_t id() const {
    return id_;
  }

  inline float operator[](std::size_t i) const {
    return position_[i];
  }

  inline const sl::point3f& position() const {
    return position_;
  }

  inline const sl::vector3f& velocity() const {
    return velocity_;
  }

  inline void set_position(const sl::point3f& pos){
    position_ = pos;  
  }

  inline void set_velocity(const sl::vector3f& vel){
    velocity_ = vel;
  }

  inline float mass() const {
    return mass_;
  }

  inline void set_mass(float m) {
    mass_ = m;
  }

  // ...
};


static void show_stats(sl::kdtree<3,float,particle>& kdt) {
  std::size_t kd_deepest;
  std::size_t kd_shortest;
  kdt.tree_depth_in(&kd_shortest, &kd_deepest);
  std::size_t kd_count = kdt.count();
  
  std::size_t kd_internal = kdt.internal_path_length();
  std::size_t kd_external = kdt.external_path_length();

  double kd_cost         = (kd_internal + kd_external) / (kd_count == 0 ? 1.0 : (double)kd_count);
  double kd_optimum_cost = kd_count == 0 ? 0.0 : std::log((double)kd_count)/std::log(2.0);
								    
  std::cerr << "    ----------------------------" << std::endl;
  std::cerr << "    KD Tree Balance Analysis: " << std::endl <<
    "        " << kd_count << " nodes, " << kd_shortest << " min depth, " << kd_deepest << " max depth." << std::endl <<
    "        " << kd_internal << " internal path length" << std::endl <<
    "        " << kd_external << " external path length" << std::endl <<
    "        " << "avg search cost = " << kd_cost << " (expected) vs. " << kd_optimum_cost << " (optimum)" <<
    std::endl;
  std::cerr << "    ----------------------------" << std::endl;
}

static void test_knn_particles(const char * id,
		     sl::kdtree<3,float,particle>& kdt) {
  typedef sl::kdtree<3,float,particle> kdt_t;

  sl::tester tester("kd-tree");

  std::cerr << std::endl;
  std::cerr << "=====================================================" << std::endl;
  std::cerr << "Test: " << id << std::endl;
  std::cerr << "=====================================================" << std::endl;
  std::cerr << std::endl;

  std::vector<particle> particles;

  // Create a random particle system
  random_seed(1);
  for (std::size_t i=0; i<20; ++i) {
    sl::vector3f v;
    float m = 5.0f;
    particles.push_back(particle(i+1,random_point(),v,m));
  }

  // Build a space partition - this *copies* the particles
  for (std::size_t i=0; i<particles.size(); ++i) {
    tester.test("find before insert", (kdt.find(particles[i]) == kdt.end()));
    kdt.insert(particles[i]);
    tester.test("kdt invariant after insert", kdt.invariant());
    tester.test("find after insert", (kdt.find(particles[i]) != kdt.end()));
  }
  
  // Test iterator
  std::size_t i=0;
  for (kdt_t::iterator it=kdt.begin(); it!=kdt.end(); ++it,++i) {
    const particle& p_it = *it;
    const particle& p_ith = *kdt.ith(i);
    tester.test("iterator vs. ith", p_it.id(), p_ith.id());
  }
  show_stats(kdt);

  // Loop on each particle and compute nearest neighbors
  std::vector<particle> knn;
  for (std::size_t i=0; i<particles.size(); ++i) {
    std::size_t pid = particles[i].id();
    sl::point3f p = particles[i].position();

    kdt.k_nearest_neighbors_in(knn,
			       p,
			       5, // Maximum number of neighbors
			       14.0f); // Maximum distance
    
    // Loop on neighbors of particle i (e.g. to compute forces)
    SL_TRACE_OUT(1) << 
      "### " << i << ": " << "p = " << p[0] << " " << p[1] << " " << p[2] << " " << std::endl;
    for (std::vector<particle>::iterator it = knn.begin();
	 it != knn.end();
	 ++it) {
      // If this is not the same particle, handle it (here, print)
      if (it->id() != pid) {	
	sl::point3f q = it->position();
	
	SL_TRACE_OUT(1) <<
	  "    " << "d2 = " << q.distance_squared_to(p) << std::endl <<
	  "    " << "     q = " << q[0] << " " << q[1] << " " << q[2] << " " << std::endl;
      }
    }
  }

  // Test erase 

  for (std::size_t i=0; i<particles.size(); i+=3) {
    tester.test("find before erase", (kdt.find(particles[i]) != kdt.end()));
    kdt.erase_exact(particles[i]);
    tester.test("kdt invariant after erase", kdt.invariant());
    tester.test("find after erase", (kdt.find(particles[i]) == kdt.end()));
  }

  kdt.clear(); 
  tester.test("kdt invariant after clear", kdt.invariant());

  failed_test_count += tester.failed_test_count();
}

std::string to_string(const sl::point2f& p) {
  return std::string() + "[" + " " + sl::to_string(p[0]) + " " + sl::to_string(p[1]) + " ]";
}

void dump_subtree(sl::kdnode<2,float,sl::point2f>* n, std::size_t level) {
  for (std::size_t l=0; l<level; ++l) std::cerr << "  ";
  std::cerr << "d=" << n->discriminant() << " p= [" << n->value()[0]<< " " << n->value()[1] << "]" << std::endl;
  if (n->lo_child()) dump_subtree(n->lo_child(), level+1);
  if (n->hi_child()) dump_subtree(n->hi_child(), level+1);
}

void test_knn_2d() {
  typedef sl::kdtree<2,float,sl::point2f> kdt_t;
  kdt_t kdt;
  
  const std::size_t N=8;
  std::vector<sl::point2f> pts;
  pts.resize(N*N);
  
  for (std::size_t i=0; i<N; ++i) {
    for (std::size_t j=0; j<N; ++j) {
      std::size_t idx = i*N+j;
      pts[idx] = sl::point2f(float(i),float(j));
      kdt.insert(pts[idx]);
    }
  }
	       
  for (std::size_t pass=0; pass<2; ++pass) {
    sl::tester tester(pass == 0 ? "kdtree priority search" : "kdtree depth first search");
    
    kdt.set_is_priority_search_enabled(pass == 0 ? true : false);

    for (std::size_t i=0; i<N-1; ++i) {
      for (std::size_t j=0; j<N-1; ++j) {
	std::size_t idx00 = (i+0)*N+(j+0); sl::point2f p00 = pts[idx00];
	std::size_t idx01 = (i+0)*N+(j+1); sl::point2f p01 = pts[idx01];
	std::size_t idx10 = (i+1)*N+(j+0); sl::point2f p10 = pts[idx10];
	std::size_t idx11 = (i+1)*N+(j+1); sl::point2f p11 = pts[idx11];

	sl::point2f pquery = (p00.lerp(p10,0.4f)).lerp((p01.lerp(p11,0.4f)),0.3f);

	std::vector<sl::point2f> knn;
	kdt.k_nearest_neighbors_in(knn,pquery,4,1e30f);

	tester.test("knn size", knn.size(), std::size_t(4)); 
	if (knn.size()==4) {
	  tester.test("knn 00", to_string(knn[0]), to_string(p00));
	  tester.test("knn 10", to_string(knn[1]), to_string(p10));
	  tester.test("knn 01", to_string(knn[2]), to_string(p01));
	  tester.test("knn 11", to_string(knn[3]), to_string(p11));
	  // i=j=N; // FIXME 
	}
      }
    }
    
    failed_test_count += tester.failed_test_count();
    if (tester.failed_test_count()) {
      std::cerr << "=== DUMP OF TREE ON ERROR" << std::endl;
      dump_subtree(kdt.begin().node_ptr(), 0);

    }
  }	  
}


int main() {
  sl::kdtree<3,float,particle> rkdt;

  test_knn_particles("RANDOMIZED KD TREE", rkdt);
  test_knn_2d();
  
  return (int)failed_test_count;
}
