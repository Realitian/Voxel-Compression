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
#ifndef SL_TRIANGLE_MESH_DISTANCE_SAMPLER_HPP
#define SL_TRIANGLE_MESH_DISTANCE_SAMPLER_HPP

#include <sl/uniform_grid_container.hpp>
#include <sl/triangle_mesh.hpp>
#include <sl/triangle_mesh_utilities.hpp>
#include <sl/random.hpp>
#include <sl/ball.hpp>

namespace sl {

  template <class T_Mesh>
  class triangle_mesh_distance_sampler {
  public:
    typedef triangle_mesh_distance_sampler<T_Mesh> this_t;
    typedef T_Mesh   mesh_t;
    
    typedef typename mesh_t::vertex_data_t::point_t point_t;
    typedef typename mesh_t::vertex_data_t::value_t value_t;
    enum {dimension = point_t::dimension };
    
    typedef typename mesh_t::const_triangle_iterator_t const_triangle_iterator_t;
    typedef typename mesh_t::const_vertex_iterator_t const_vertex_iterator_t;
    typedef typename mesh_t::const_edge_iterator_t const_edge_iterator_t;
    
    typedef fixed_size_vector<sl::column_orientation,dimension,value_t> vector_t;
    typedef axis_aligned_box<dimension,value_t>                         aabox_t;

    typedef uniform_grid_container<dimension, value_t, const_triangle_iterator_t> ugrid_t;
    
    typedef enum face_sampling_mode {
      MONTECARLO_FACE_SAMPLING_MODE,
      SUBDIVISION_FACE_SAMPLING_MODE,
      SIMILAR_TRIANGLES_FACE_SAMPLING_MODE
    } face_sampling_mode_t;

  protected:
    // data structures
    const mesh_t&   S1_; 
    const mesh_t&   S2_;
    ugrid_t         gS2_;

    // parameters
    value_t         n_samples_per_area_unit_;
    int             n_samples_target_;

    bool            vertex_sampling_enabled_;
    bool            edge_sampling_enabled_;
    bool            face_sampling_enabled_;
    face_sampling_mode_t face_sampling_mode_;
    
    // results
    unsigned long   n_total_samples_;
    unsigned long   n_total_area_samples_;
    unsigned long   n_total_edge_samples_;
    unsigned long   n_total_vertex_samples_;
    value_t         max_dist_;
    value_t         mean_dist_; 
    value_t         RMS_dist_;
    value_t         volume_;
    value_t         area_S1_;
    
    // globals
    int             n_samples_;

    // Random
    mutable sl::random::uniform_closed<value_t>  rng_;
    
  public:

    triangle_mesh_distance_sampler(const mesh_t& m1, const mesh_t& m2)
        : S1_(m1), S2_(m2) {
      vertex_sampling_enabled_ = true;   // FIXME
      edge_sampling_enabled_   = false;
      face_sampling_enabled_   = true;
      face_sampling_mode_      = MONTECARLO_FACE_SAMPLING_MODE;
      set_sample_count_on_target(std::max(S1_.triangle_count(),
                                          S2_.triangle_count()));
      area_S1_ = triangle_mesh_area(S1_);
      // FIXME Init others
      if ( area_S1_ < value_t(1e-6) ) {
	SL_TRACE_OUT( -1 ) << "computed mesh area " << area_S1_ << " changed to " << 1e-6 << std::endl;
	area_S1_ = 1e-6;	
      }
    }

  protected:
    
    value_t computed_distance_point_surface(const point_t& p,
                                            const mesh_t&  m,
                                            const typename ugrid_t::object_list_t& tlist,
                                            const value_t& dist_upper_bound) const {
#if 0
      std::size_t fixme_same = 0;
      std::size_t fixme_count = 0;
#endif
      value_t result = dist_upper_bound;
      for (typename ugrid_t::object_list_t::const_iterator tit = tlist.begin();
           tit != tlist.end();
           ++tit) {
        const point_t& p0 = m.vertex_data((*tit)->first[0])->position();
        const point_t& p1 = m.vertex_data((*tit)->first[1])->position();
        const point_t& p2 = m.vertex_data((*tit)->first[2])->position();
        
        value_t d = p.distance_to(p0,p1,p2);
        if (d<result) result = d;
#if 0
        // FIXME CHECK
        ++fixme_count;
        if (p==p0 || p==p1 || p==p2) {
          ++fixme_same;
          if (d > 1e-7) {
            SL_TRACE_OUT(-1) << "Distance error!: " <<
              "p = " << p[0] << " " << p[1] << " " << p[2] << std::endl <<
              "p0 = " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl <<
              "p1 = " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl <<
              "p2 = " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl <<
              "=====> d = " << d << " vs. " << 0.0f << std::endl;
          }
        }
#endif
      }
      return result;
    }
    
    value_t computed_distance_point_surface(const point_t& p,
                                            const mesh_t&  m,
                                            const ugrid_t& ugrid_m,
                                            const value_t& dist_upper_bound) const {
      typename ugrid_t::subscript_t idx0 = ugrid_m.subscript(p);
      value_t result = computed_distance_point_surface(p,m,ugrid_m.item(idx0),dist_upper_bound);
      value_t voxdim = ugrid_m.voxel_dimensions()[ugrid_m.voxel_dimensions().iamin()];
      int step = 1;
      if ( voxdim < value_t( 1e-6 ) ) {
	SL_TRACE_OUT(-1) << "voxdim = " << voxdim << " changed to " << value_t( 1e-6 ) << std::endl;
	voxdim = value_t( 1e-6 );
      }
      while ((step-1)*voxdim <= result) {
        typename ugrid_t::subscript_t idx_l;
        typename ugrid_t::subscript_t idx_h;
        for (std::size_t i=0; i<dimension; ++i) {
	  if (int(idx0[i]) >= step) {
	    idx_l[i] = idx0[i]-step;
	  } else {
	    idx_l[i] = 0;
	  }
          idx_h[i] = idx0[i]+step;
	  if (idx_h[i] >= ugrid_m.extent()[i]) {
	    idx_h[i] = ugrid_m.extent()[i] - 1;
	  }
        }
        for (typename ugrid_t::subscript_t idx = idx_l; idx.all_le(idx_h); idx.increment(idx_l,idx_h)) {
          if (ugrid_m.good_subscript(idx)) {
            bool outer_layer = false;
            for (std::size_t i=0; i<dimension && !outer_layer; ++i) {
              outer_layer = (idx[i] == idx_l[i]) || (idx[i] == idx_h[i]);
            }
            if (outer_layer) {
              if (ugrid_m.squared_point_cell_distance(p, idx) < result*result) {
                result = computed_distance_point_surface(p,m,ugrid_m.item(idx),result);
              }
            }
          }
        }
        ++step;
	if ( step > 10000000 ) {
	  SL_TRACE_OUT(-1) << "high step value: " << step << ", result: " << result << ", voxdimstd: " << voxdim << std::endl;
	  break;
	}
      }
      return result;
    }
    
  protected:
    
    value_t add_sample(const point_t &p) {
      value_t dist_upper_bound = gS2_.bounding_box().diagonal().two_norm();
      value_t dist = computed_distance_point_surface(p, S2_, gS2_, dist_upper_bound);
      if (dist == dist_upper_bound) return -sl::scalar_math<value_t>::one();     
      if(dist > max_dist_)
        max_dist_ = dist;        // L_inf
      mean_dist_ += dist;	 // L_1
      RMS_dist_  += dist*dist;   // L_2
      n_total_samples_++;

      return dist;    
    }
    
    void add_random_sample(const mesh_t&  m,
                           const_triangle_iterator_t& t_it) {
      const point_t& p0 = m.vertex_data(t_it->first[0])->position();
      const point_t& p1 = m.vertex_data(t_it->first[1])->position();
      const point_t& p2 = m.vertex_data(t_it->first[2])->position();

      const vector_t e1 = p1-p0;
      const vector_t e2 = p2-p0;

      value_t rnd_1 = rng_.value();
      value_t rnd_2 = rng_.value();
      if (rnd_1 + rnd_2 > sl::scalar_math<value_t>::one()) {
        rnd_1 = sl::scalar_math<value_t>::one() - rnd_1;
        rnd_2 = sl::scalar_math<value_t>::one() - rnd_2;
      }

      add_sample(p0+rnd_1*e1+rnd_2*e2);
      
      n_total_area_samples_++;
    }

    void add_edge_samples(const point_t & v0, const point_t & v1, int n_samples_per_edge) {
      // uniform sampling of the segment v0v1.
      vector_t    e = (v1-v0)/(value_t)(n_samples_per_edge+1);
      for(int i=1; i <= n_samples_per_edge; i++) {
        add_sample(v0 + value_t(i)*e);
        n_total_edge_samples_++;
      }
    }

  protected:
    
    void do_vertex_sampling() {
      for(const_vertex_iterator_t v_it=S1_.vertex_map().begin();
          v_it!=S1_.vertex_map().end();
          ++v_it) {
        add_sample(v_it->second.data()->position());
        n_total_vertex_samples_++;
      }
    }
    
    void do_edge_sampling() {
      value_t n_samples_per_length_unit =
        face_sampling_enabled() ? std::sqrt(n_samples_per_area_unit_) : n_samples_per_area_unit_;
      value_t n_samples_decimal = sl::scalar_math<value_t>::zero();
      for(const_edge_iterator_t e_it=S1_.edge_map().begin();
          e_it!=S1_.edge_map().end();
          ++e_it) {
        const point_t& p0 = S1_.vertex_data(e_it->first[0])->position();
        const point_t& p1 = S1_.vertex_data(e_it->first[1])->position();

        n_samples_decimal += p0.distance_to(p1) * n_samples_per_length_unit;
        n_samples_      = (int) n_samples_decimal;
        add_edge_samples(p0, p1, (int) n_samples_);
        n_samples_decimal -= (value_t) n_samples_;
      }
    }
    
    void do_face_subdivision_sampling(const point_t &v0,
                                      const point_t &v1,
                                      const point_t &v2,
                                      int maxdepth) {
      if (maxdepth == 0) {
        add_sample(v0 + (v1-v0 + v2-v0)/3.0f);
        n_total_area_samples_++;
        n_samples_++;
      } else {
        // Split at longest edge
        value_t  maxd01 = v0.distance_squared_to(v1);
        value_t  maxd12 = v1.distance_squared_to(v2);
        value_t  maxd20 = v2.distance_squared_to(v0);
        int     ilongest;
        if(maxd01 > maxd12) 
          if(maxd01 > maxd20)     ilongest = 0;
          else                    ilongest = 2;
        else
          if(maxd12 > maxd20)     ilongest = 1;
          else                    ilongest = 2;
        switch(ilongest) {
        case 0 : {
          point_t pp = v0+0.5f*(v1-v0); 
          do_face_subdivision_sampling(v0,pp,v2,maxdepth-1);
          do_face_subdivision_sampling(pp,v1,v2,maxdepth-1);
        } break;
        case 1 : {
          point_t pp = v1+0.5f*(v2-v1);
          do_face_subdivision_sampling(v0,v1,pp,maxdepth-1);
          do_face_subdivision_sampling(v0,pp,v2,maxdepth-1);
        } break;
        case 2 : {
          point_t pp = v2+0.5f*(v0-v2); 
          do_face_subdivision_sampling(v0,v1,pp,maxdepth-1);
          do_face_subdivision_sampling(pp,v1,v2,maxdepth-1);
        } break;
        }
      }
    }

    void do_similar_triangles_sampling(const point_t &v0,
                                       const point_t &v1,
                                       const point_t &v2,
                                       int n_samples_per_edge) {
      vector_t e1 = (v1-v0)/(value_t)(std::max(n_samples_per_edge-1,1));
      vector_t e2 = (v2-v0)/(value_t)(std::max(n_samples_per_edge-1,1));

      for(int i=1; i < n_samples_per_edge-1; i++) {
        for(int j=1; j < n_samples_per_edge-1-i; j++) {
          add_sample(v0 + value_t(i)*e1 + value_t(j)*e2);
          n_total_area_samples_++;
          n_samples_++;
        }
      }
    }
    
    void do_montecarlo_face_sampling() {
      value_t  n_samples_decimal = 0.0;
      for (const_triangle_iterator_t t_it = S1_.triangle_map().begin();
           t_it != S1_.triangle_map().end();
           ++t_it) {
        const point_t& p0 = S1_.vertex_data(t_it->first[0])->position();
        const point_t& p1 = S1_.vertex_data(t_it->first[1])->position();
        const point_t& p2 = S1_.vertex_data(t_it->first[2])->position();

        const value_t A = std::sqrt(sl::triangle_area_squared(p0,p1,p2));

        n_samples_decimal += A * n_samples_per_area_unit_;
        n_samples_         = (int) n_samples_decimal;
        for(int i=0; i < n_samples_; i++) {
          add_random_sample(S1_, t_it);
        }

        n_samples_decimal -= (value_t) n_samples_;
      }
    }
    
    void do_face_subdivision_sampling() {
      value_t  n_samples_decimal = 0.0;
      for (const_triangle_iterator_t t_it = S1_.triangle_map().begin();
           t_it != S1_.triangle_map().end();
           ++t_it) {
        const point_t& p0 = S1_.vertex_data(t_it->first[0])->position();
        const point_t& p1 = S1_.vertex_data(t_it->first[1])->position();
        const point_t& p2 = S1_.vertex_data(t_it->first[2])->position();

        const value_t A = std::sqrt(sl::triangle_area_squared(p0,p1,p2));

        n_samples_decimal += A * n_samples_per_area_unit_;
        n_samples_         = (int) n_samples_decimal;
        if(n_samples_) {
          // face sampling.
          int maxdepth = (int)(std::log((value_t)n_samples_)/log(2.0));
          n_samples_ = 0;
          do_face_subdivision_sampling(p0,p1,p2, maxdepth);
        }
        n_samples_decimal -= n_samples_;
      }
    }
    
    void do_similar_face_sampling() {
      value_t  n_samples_decimal = 0.0;
      for (const_triangle_iterator_t t_it = S1_.triangle_map().begin();
           t_it != S1_.triangle_map().end();
           ++t_it) {
        const point_t& p0 = S1_.vertex_data(t_it->first[0])->position();
        const point_t& p1 = S1_.vertex_data(t_it->first[1])->position();
        const point_t& p2 = S1_.vertex_data(t_it->first[2])->position();

        const value_t A = std::sqrt(sl::triangle_area_squared(p0,p1,p2));

        n_samples_decimal += A * n_samples_per_area_unit_;
        n_samples_         = (int) n_samples_decimal;
        if (n_samples_) {
          int n_samples_per_edge = (int)((std::sqrt(1.0+8.0*(value_t)n_samples_)+5.0)/2.0);
          n_samples_ = 0;
          do_similar_triangles_sampling(p0, p1, p2, n_samples_per_edge);
        }
        n_samples_decimal -= (value_t) n_samples_;
      }
    }

    void do_face_sampling() {
      switch (face_sampling_mode()) {
      case MONTECARLO_FACE_SAMPLING_MODE: do_montecarlo_face_sampling(); break;
      case SUBDIVISION_FACE_SAMPLING_MODE: do_face_subdivision_sampling(); break;
      case SIMILAR_TRIANGLES_FACE_SAMPLING_MODE: do_similar_face_sampling(); break;
      default: SL_FAIL("Unknown face sampling mode"); break;
      }
    }
    
    void do_sampling() {
      // initialize sampling statistics.
      n_total_area_samples_   = 0;
      n_total_edge_samples_   = 0;
      n_total_vertex_samples_ = 0;
      n_total_samples_        = 0;
      n_samples_              = 0;
      max_dist_               = sl::scalar_math<value_t>::lower_bound();
      mean_dist_              = 0;
      RMS_dist_               = 0;

      SL_TRACE_OUT(1) << "Sampling start: n_samples_to_go: " << n_samples_target_ << std::endl;

      // Sample!!
      if (vertex_sampling_enabled()) {
        do_vertex_sampling();
        n_samples_target_ -= (int) n_total_samples_;
        SL_TRACE_OUT(1) << "Sampled vertices: n_samples_to_go: " << n_samples_target_ << std::endl;
      }
      
      if (n_samples_target_ > 0) {
        n_samples_per_area_unit_ = (value_t)(n_samples_target_) / area_S1_;
        if (edge_sampling_enabled()) {
          do_edge_sampling();
          n_samples_target_ -= (int) n_total_samples_;
          SL_TRACE_OUT(1) << "Sampled edges: n_samples_to_go: " << n_samples_target_ << std::endl;
        }
        if (n_samples_target_ > 0) {
          if (face_sampling_enabled()) {
            n_samples_per_area_unit_  = value_t(n_samples_target_) / area_S1_;
            do_face_sampling();
            n_samples_target_ -= (int) n_total_samples_;
            SL_TRACE_OUT(1) << "Sampled faces: n_samples_to_go: " << n_samples_target_ << std::endl;
          }
        }
      }

      // compute statistics
      n_samples_per_area_unit_ = value_t(n_total_samples_) / area_S1_;
      volume_                  = mean_dist_ / n_samples_per_area_unit_ / value_t(2.0f);
      mean_dist_              /= value_t(n_total_samples_);
      RMS_dist_                = std::sqrt(RMS_dist_ / n_total_samples_);
    }
    
  public:
    
    void compute() {
      // --  Find a good bbox
      aabox_t bbox;
      bbox = triangle_mesh_aabox(S1_);
      bbox.merge(triangle_mesh_aabox(S2_));
      vector_t hsl = 1.01f*bbox.half_side_lengths();
      value_t hsl_max = hsl[hsl.iamax()];
      value_t hsl_eps = std::max(value_t(1e-5), value_t(hsl_max * 0.0001));
      for (std::size_t i=0; i<dimension; ++i) {
        hsl[i] = std::max(hsl[i], hsl_eps);
      }
      bbox = aabox_t(bbox.center()-hsl, bbox.center()+hsl);

      // -- Find a good cell size and engrid S2      
      value_t cell_sz = std::sqrt(2.0f/sqrt(3.0f)*area_S1_/sl::max(S2_.triangle_count(), std::size_t(1)));
      SL_TRACE_OUT(1) << "Resize: " << cell_sz << std::endl;
      gS2_.clear();
      gS2_.resize(bbox,cell_sz);

      SL_TRACE_OUT(1) << "Grid: " << gS2_.extent() << std::endl;
      SL_TRACE_OUT(1) << "VOX : " << gS2_.voxel_dimensions()[0] << " " << gS2_.voxel_dimensions()[1] << " " << gS2_.voxel_dimensions()[2] << std::endl;
      
      for (const_triangle_iterator_t t_it = S2_.triangle_map().begin();
           t_it != S2_.triangle_map().end();
           ++t_it) {
        const point_t& p0 = S2_.vertex_data(t_it->first[0])->position();
        const point_t& p1 = S2_.vertex_data(t_it->first[1])->position();
        const point_t& p2 = S2_.vertex_data(t_it->first[2])->position();
        aabox_t t_box;
        t_box.to(p0);
        t_box.merge(p1);
        t_box.merge(p2);
        gS2_.insert(t_box, t_it);
      }
      SL_TRACE_OUT(1) << "Inserted : " << S2_.triangle_count() << " triangles." << std::endl;

      //-- Sample!
      do_sampling();

      SL_TRACE_OUT(1) << " Area     : " << area() << std::endl;
      SL_TRACE_OUT(1) << " Dist max : " << max_distance()  << std::endl;
      SL_TRACE_OUT(1) << " Dist mean: " << mean_distance() << std::endl;
      SL_TRACE_OUT(1) << " Dist rms : " << rms_distance() << std::endl;
      SL_TRACE_OUT(1) << " Dist vol : " << distance_volume() << std::endl;
      
      SL_TRACE_OUT(1) << " Total sample count   : " << total_sample_count() << std::endl;
      SL_TRACE_OUT(1) << " Area samples count   : " << area_samples_count() << std::endl;
      SL_TRACE_OUT(1) << " Edge samples count   : " << edge_samples_count() << std::endl;
      SL_TRACE_OUT(1) << " Vertex samples count : " << vertex_samples_count() << std::endl;
      SL_TRACE_OUT(1) << " Samples/area unit    : " << sample_count_per_area_unit() << std::endl;
      SL_TRACE_OUT(1) << " Samples/target       : " << sample_count_on_target() << std::endl;
      
      //-- Cleanup!
      gS2_.clear();
    }
    
    value_t          area() const  { return area_S1_; }
    value_t          max_distance() const { return max_dist_; }
    value_t          mean_distance() const { return mean_dist_; }
    value_t          rms_distance() const  { return RMS_dist_; }
    value_t          distance_volume() const  { return volume_; }

    std::size_t      total_sample_count() const   { return n_total_samples_; }
    std::size_t      area_samples_count() const   { return n_total_area_samples_; }
    std::size_t      edge_samples_count() const   { return n_total_edge_samples_; }
    std::size_t      vertex_samples_count() const { return n_total_vertex_samples_; }

    value_t          sample_count_per_area_unit() const { return n_samples_per_area_unit_; }
    int              sample_count_on_target() const { return n_samples_target_; }

  public: // Config

    bool vertex_sampling_enabled() const {
      return vertex_sampling_enabled_;
    }
    
    void set_vertex_sampling_enabled(bool b = true) {
      vertex_sampling_enabled_ = b;
    }
    
    bool edge_sampling_enabled() const {
      return edge_sampling_enabled_;
    }

    void set_edge_sampling_enabled(bool b = true) {
      edge_sampling_enabled_ = b;
    }
    
    bool face_sampling_enabled() const {
      return face_sampling_enabled_;
    }

    void set_face_sampling_enabled(bool b = true) {
      face_sampling_enabled_ = b;
    }

    face_sampling_mode_t face_sampling_mode() const {
      return face_sampling_mode_;
    }
    
    void  set_face_sampling_mode(const face_sampling_mode_t& x) {
      face_sampling_mode_ = x;
    }
    
    void set_sample_count_on_target(int n_samp) {
      n_samples_target_        = n_samp;
      n_samples_per_area_unit_ = (value_t) n_samples_target_ / area_S1_;
    }

    void set_sample_count_per_area_unit(value_t n_samp) {
      n_samples_per_area_unit_ = n_samp;
      n_samples_target_        = (int)((value_t) n_samples_per_area_unit_ * area_S1_);
    }
  };


  template <class mesh_t>
  typename triangle_mesh_distance_sampler<mesh_t>::value_t triangle_mesh_forward_hausdorff_distance(const mesh_t& M1,
                                                                                                           const mesh_t& M2) {
    typedef typename triangle_mesh_distance_sampler<mesh_t>::value_t value_t;
    if (M1.triangle_count() == 0) {
      if (M2.triangle_count() == 0) {
        // Null mesh to null mesh: assume zero
        return scalar_math<value_t>::zero();
      } else {
        // Mesh to null mesh: assume a "collapse" to bsphere center
        return triangle_mesh_minimum_enclosing_ball(M2).radius();
      }
    } else if (M2.triangle_count() == 0) {
        // Mesh to null mesh: assume a "collapse" to bsphere center
      return triangle_mesh_minimum_enclosing_ball(M1).radius();
    } else {
      triangle_mesh_distance_sampler<mesh_t> forward_sampler(M1,M2);
      forward_sampler.compute();
      return forward_sampler.max_distance();
    }
  }

  template <class mesh_t>
  typename triangle_mesh_distance_sampler<mesh_t>::value_t triangle_mesh_symmetric_hausdorff_distance(const mesh_t& M1,
                                                                                                             const mesh_t& M2) {
    typedef typename triangle_mesh_distance_sampler<mesh_t>::value_t value_t;
    if (M1.triangle_count() == 0) {
      if (M2.triangle_count() == 0) {
        // Null mesh to null mesh: assume zero
        return scalar_math<value_t>::zero();
      } else {
        // Mesh to null mesh: assume a "collapse" to bsphere center
        return triangle_mesh_minimum_enclosing_ball(M2).radius();
      }
    } else if (M2.triangle_count() == 0) {
      // Mesh to null mesh: assume a "collapse" to bsphere center
      return triangle_mesh_minimum_enclosing_ball(M1).radius();
    } else {
      return std::max(triangle_mesh_forward_hausdorff_distance(M1,M2),
                      triangle_mesh_forward_hausdorff_distance(M2,M1));
    }
  }

}

#endif

