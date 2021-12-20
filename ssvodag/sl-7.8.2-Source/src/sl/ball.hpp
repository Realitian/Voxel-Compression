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
#ifndef SL_BALL_HPP
#define SL_BALL_HPP

#include <sl/bounding_volume.hpp>
#include <sl/linear_map.hpp>
#include <sl/fixed_size_plane.hpp>
#include <sl/interval.hpp>
#include <sl/random.hpp>
#include <sl/projective_map.hpp>
#include <list>

namespace sl {

  namespace detail {
    template <size_t N_dimension, class T> 
    struct ball_data {
      typedef fixed_size_point<N_dimension,T> point_t;
      typedef typename point_t::value_t       value_t;
      
      point_t center_;
      value_t radius_;

      inline ball_data() : radius_(scalar_math<value_t>::zero()) {};
      inline ball_data(const point_t& p,
                       const value_t& r) : center_(p), radius_(r) {
      }
    public: // Serialization

      void store_to(output_serializer& s) const {
        s << center_ << radius_;
      }
    
      void retrieve_from(input_serializer& s) {
        s >> center_ >> radius_;
      }
    };
  }

  
  /**
   * A N-dimensional ball
   */
  template <size_t N_dimension, class T> 
  class ball: 
    public 
      bounding_volume_base<
        N_dimension,
        T,
        ball<N_dimension,T>,
        detail::ball_data<N_dimension,T> 
      > 
  {
  public:
    enum { dimension = N_dimension };

    typedef bounding_volume_base<
        N_dimension,
        T,
        ball<N_dimension,T>,
        detail::ball_data<N_dimension,T> 
    > super_t;

    typedef typename super_t::derived_t         this_t;
    typedef typename super_t::value_t           value_t;
    typedef typename super_t::data_t            data_t;
    typedef typename super_t::point_t           point_t;
    typedef typename super_t::line_t            line_t;
    typedef typename super_t::ray_t             ray_t;
    typedef typename super_t::line_segment_t    line_segment_t;

    typedef fixed_size_vector<column_orientation,N_dimension,T> vector_t;
    typedef fixed_size_plane<N_dimension,T>                     plane_t;


  protected: // Helpers
    
    static inline size_t factorial(size_t x) {
      size_t result = 1;
      for (size_t i=1; i<=x; ++i) {
        result *= i;
      }
      return result;
    }

    static inline size_t double_factorial(size_t x) {
      size_t result = 1;
      for (size_t i=2-(x%2); i<=x; i+=2) {
        result *= i;
      }
      return result;
    }

  public: // Creation & Destruction

    /// Default init (positioned at the origin, zero size)
    inline ball() : super_t(data_t()) { 
    }

    /// Init from params
    inline ball(const point_t& p,
                const value_t& r): super_t(data_t(p,r)) {
    }

    /// The ball containing the single point p
    inline explicit ball(const point_t& p): super_t(data_t(p,scalar_math<value_t>::zero())) {
    }

  public: // Access

    /// The ball center
    inline point_t& center() {
      return (this->data_).center_; 
    }
	
    /// The ball center
    inline const point_t& center() const { 
      return (this->data_).center_;
    }

    /// The ball radius
    inline value_t& radius() {
      return (this->data_).radius_; 
    }
	
    /// The ball radius
    inline const value_t& radius() const { 
      return (this->data_).radius_;
    }

  public: // Queries
		
    /// True iff the ball is empty (i.e. radius is negative)
    inline bool is_empty() const {
      return sl::is_negative(radius());
    }

    /// True iff the ball is a single point (i.e. radius is zero)
    inline bool is_point() const {
      return sl::is_zero(radius());
    }
    
    /// the volume of the ball
    value_t volume() const {
      value_t result = scalar_math<value_t>::zero();
      if (!is_empty()) {
        // N-dimensional formula...
        // see http://mathworld.wolfram.com/Hypersphere.html

        value_t S_n;
        if (dimension % 2 == 1) {
          // Odd dimension
          S_n =
            value_t(pow(2.0,(dimension+1)/2))* pow(scalar_math<value_t>::Pi(),(dimension-1)/2) /
            value_t(double_factorial(sl::max((int)dimension,2)-2));
        } else {
          // Even dimension
          S_n =
            scalar_math<value_t>::two() * pow(scalar_math<value_t>::Pi(),dimension/2) /
            value_t(factorial(dimension/2-1));
        }
        
        result = S_n * pow(radius(),dimension)/value_t(dimension);
      }

      return result;
    }

    /// The surface area of the ball
    value_t surface_area() const {
      value_t result = scalar_math<value_t>::zero(); 
      SL_FAIL("Not implemented");
      return result;
    }

    /// The projected area of the ball in direction d
    value_t projected_area(const vector_t& d) const {
      SL_USEVAR(d);
      value_t result = scalar_math<value_t>::zero();
      SL_FAIL("Not implemented");
      return result;
    }

    /// The square of the minimum Euclidean distance to the point p
    value_t square_distance_to(const point_t& p) const {
      return sqr(distance_to(p));
    }

    /// The minimum Euclidean distance to the point p
    inline value_t distance_to(const point_t& p) const {
      value_t result = scalar_math<value_t>::zero(); 
      value_t d2_pc = p.distance_squared_to(center());
      if (d2_pc >= radius()*radius()) {
        result = std::sqrt(d2_pc) - radius();
      }
      return result;
    }

    /// The range of signed distances to the plane hp
    inline interval<value_t> signed_distance_range(const plane_t& hp) const {
      value_t d_center = hp.signed_distance(center());
      return interval<value_t>(d_center - radius(), d_center + radius());
    }

    /// Is the ball fully on the side of the plane indicated by its normal? 
    inline bool is_fully_above(const plane_t& hp) const {
      return sl::is_positive(hp.signed_distance(center()) - radius());
    }

    /// Is the ball fully on the side of the plane indicated by its normal? 
    inline bool is_fully_below(const plane_t& hp) const {
      return sl::is_non_positive(hp.signed_distance(center()) + radius());
    }

  public: // Containment 

    /// True iff p is inside this
    bool contains(const point_t& p) const {
      return (!is_empty()) && (p.distance_squared_to(center()) <= sqr(radius()));
    }

    /// True if this scaled by 1+eps contains p 
    bool epsilon_contains(const point_t& p, const value_t& eps = value_t(0.001)) const {
      const value_t One = scalar_math<value_t>::one();
      return (!is_empty()) && (p.distance_squared_to(center()) <= sqr(radius()*(One+eps)));
    }

    /// True iff other is inside this
    bool contains(const this_t& other) const {
      if (other.radius() > radius()) {
        return false;
      } else {
        return (other.center() - center()).two_norm_squared() <= sqr(radius()-other.radius());
      }
    }

    /// True if this scaled by 1+eps contains other 
    bool epsilon_contains(const this_t& other, const value_t& eps = value_t(0.001)) const {
      if (other.radius() > radius()) {
        return false;
      } else {
        const value_t One = scalar_math<value_t>::one();
        return (other.center() - center()).two_norm_squared() <= sqr(radius()*(One+eps)-other.radius());
      }
    }
    
  public: // Ray-tracing

    /**
     * Clip a linear component with this axis aligned ball.
     * on entry, tmin and tmax determine the extent of
     * the linear component. On exit, off is true if
     * the linear component is outside of the ball, and
     * tmin, tmax contain the portion of the component
     * inside the ball.
     */
    void clip(const point_t& org,
              const vector_t& dir,
              bool&    off,
              value_t& tmin,
              value_t& tmax) const {
      value_t A = dir.two_norm_squared();

      value_t B = scalar_math<value_t>::zero();
      for (std::size_t d=0; d<dimension; ++d) {
        B += dir[d] * (org[d] - center()[d]); 
      }
      B+= B;
            
      value_t C = org.distance_squared_to(center()) - radius()*radius(); 

      value_t Delta = B*B - 4*A*C;
      off = is_negative(Delta);
      if (!off) {
        value_t Sqrt_Delta = std::sqrt(Delta);
        value_t t0 = (-B - Sqrt_Delta) / (2*A);
        value_t t1 = (-B + Sqrt_Delta) / (2*A);

        tmin = sl::max(tmin, t0);
        tmax = sl::min(tmax, t1);
        off = (tmin > tmax);
      }
    }

    /// True iff the ray intersects the volume
    inline bool intersection_exists(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }


    /// True iff the line intersects the volume
    inline bool intersection_exists(const line_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }

    /// True iff the line segment intersects the volume
    inline bool intersection_exists(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return !off;
    }

    /**
     * ray - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest ray parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      if (off) {
        return std::pair<bool,value_t>(false,tmin);
      } else if (is_zero(tmin)) {
        // No` near intersection farther than epsilon away, look at far intersection
        return std::pair<bool,value_t>(true,tmax);
      } else {
        return std::pair<bool,value_t>(true,tmin);
      }
    }

    /**
     * line - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest line parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const line_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

    /**
     * segment - hollow volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    inline std::pair<bool,value_t> intersection(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      if (off) {
        return std::pair<bool,value_t>(false,tmin);
      } else if (is_zero(tmin)) {
        return std::pair<bool,value_t>(!is_one(tmax),tmax);
      } else {
        return std::pair<bool,value_t>(true,tmin);
      }
    }

    /**
     * ray - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest ray parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const ray_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

    /**
     * line - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest line parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const line_t& r) const {
      return intersection(r);
    }

    /**
     * segment - volume intersection. The routine returns the 
     * pair <false,???> if there is no intersection, and the pair
     * <true, t> if an intersection occurs. t is the value of the
     * smallest segment parameter for an intersection point.
     */
    inline std::pair<bool,value_t> bounds_intersection(const line_segment_t& r) const {
      bool    off  = false;
      value_t tmin = r.tmin();
      value_t tmax = r.tmax();
      clip(r.origin(),r.direction(),off,tmin,tmax);
      return std::pair<bool,value_t>(!off,tmin);
    }

  public: // Commands

    /// Set this to the zero-sized ball containing p
    inline void to(const point_t& p) {
      center() = p;
      radius() = scalar_math<value_t>::zero();
    }

    /// Set this to the empty ball
    inline void to_empty() {
      center() =  scalar_math<value_t>::zero();
      radius() = -scalar_math<value_t>::one();
    }

    /// Set this to the zero-sized ball containing the origin
    inline void to_zero() {
      to(point_t());
    }
	
    /// Set this to the largest representable ball
    inline void to_huge() {
      center() = scalar_math<value_t>::zero();
      radius() = std::sqrt(scalar_math<value_t>::upper_bound());
    }

    /// Add point p to this, eventually adjusting the center and radius of this
    void merge(const point_t& p) {
      if (is_empty()) {
        to(p);
      } else {
        vector_t pc = p - center();
        value_t  pc_len2 = pc.two_norm_squared();
        if (pc_len2 > radius()*radius()) {
          // this point is outside of current sphere
          value_t pc_len = std::sqrt(pc_len2);
          (this->data_).radius_ = value_t(0.5) * ((this->data_).radius_ + pc_len); 
          (this->data_).center_ = (this->data_).center_ + (pc_len - (this->data_).radius_)/pc_len * pc;
        }
      }
    }
      
    /// Merge ball other with this, eventually adjusting the center and radius of this
    void merge(const this_t& other) {
      if (other.is_empty()) {
        // nothing to do
      } else if (is_empty()) {
        (*this) = other;
      } else {
#if 0
        vector_t pc = other.center() - center();
        value_t  pc_len = pc.two_norm();
        if (pc_len+other.radius() > radius()) {
          (this->data_).radius_ = value_t(0.5) * ((this->data_).radius_ + pc_len + other.radius()); 
          (this->data_).center_ = (this->data_).center_ + (pc_len + other.radius() - (this->data_).radius_)/pc_len * pc;
        }
#else
        vector_t c1c2 = other.center() - center();
        value_t c1c2_len2 = c1c2.two_norm_squared();
        value_t dr = other.radius() - radius();

        if ( dr*dr >= c1c2_len2 ) {
          if (sl::is_positive(dr)) {
            *this = other;
          }
        } else {
          value_t c1c2_len = std::sqrt(c1c2_len2);
          if (sl::is_positive(c1c2_len)) {
            (this->data_).center_ += (c1c2_len + dr)/(scalar_math<value_t>::two()*c1c2_len)*c1c2;
          } 
          (this->data_).radius_ = (c1c2_len + radius() + other.radius())/scalar_math<value_t>::two();
        }
#endif
      }
    }

    /// Add point p to this, eventually adjusting the radius of this
    void grow(const point_t& p) {
      if (is_empty()) {
        to(p);
      } else {
        value_t d2_pc = p.square_distance_to(center());
        if (d2_pc > radius()*radius()) {
          (this->data_).radius_ = std::sqrt(d2_pc);
        }
      }
    }
      
    /// Merge ball other with this, eventually adjusting the radius of this
    void grow(const this_t& other) {
      if (other.is_empty()) {
        // nothing to do
      } else if (is_empty()) {
        (*this) = other;
      } else {
        value_t r = center().distance_to(other.center())+other.radius();
        if (r > radius()) {
          (this->data_).radius_ = r;
        }
      }
    }

  public: // Transformations

    /// the ball containing this transformed by map 
    template <class T_tag>
    this_t transformed_by(const sl::linear_map_base<T_tag, dimension, value_t>& m) const {
      SL_USEVAR(m);
      this_t result;
      SL_FAIL("Not implemented");
      return result;
    }
    
  public: // Comparison

    /// -1 if this < t2, +1 if this > t2, 0 otherwise (sequential element comparison)
    inline int compare(const this_t& t2) const {
      int result = center().compare(t2.center());
      if (result == 0) {
        if (radius() < t2.radius()) {
          result = -1;
        } else if (radius() > t2.radius()) {
          result = 1;
        }
      }
      return result;
    }

    /// is this < t2 (sequential element comparison
    inline bool operator<(const this_t& t2) const {
      return compare(t2) < 0;
    }

    /// is this equal to t2?
    inline bool operator == (const this_t& t2) const {
      return compare(t2) == 0;
    }

    SL_OP_COMPARABLE1(this_t);
    SL_OP_EQUALITY_COMPARABLE1(this_t);

  }; // ball

  template <std::size_t DIM, typename OUT_ET>
  class conv_to< ball<DIM, OUT_ET> > {
  public:
    typedef ball<DIM, OUT_ET> result_t;

    // Explicit conversion from arrays of another type
    template <typename IN_ET> 
    inline static result_t from(const ball<DIM, IN_ET>& in) {
      return result_t(conv_to< fixed_size_point<DIM,OUT_ET> >(in.center()),
		      static_cast<OUT_ET>(in.radius()));
    }
  };
  
} // namespace sl

/// the axis_aligned ball b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::ball<N_dimension,T> transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
                                              const sl::ball<N_dimension,T>& b) {
  return b.transformed_by(map);
}

/// the axis_aligned ball b transformed by the inverse of map
template <class T_tag, size_t N_dimension, class T>
inline sl::ball<N_dimension,T> inverse_transformation(const sl::linear_map_base<T_tag, N_dimension, T>& map, 
                                                      const sl::ball<N_dimension,T>& b) {
  return b.transformed_by(map.inverse());
}

/// the ball b transformed by map
template <class T_tag, size_t N_dimension, class T>
inline sl::ball<N_dimension,T> operator * (const sl::linear_map_base<T_tag, N_dimension, T>& map, 
                                           const sl::ball<N_dimension,T>& b) {
  return b.transformed_by(map);
}

// I/O

template <size_t N_dimension, class T>
std::ostream& operator <<(std::ostream& s, const sl::ball<N_dimension,T>& b) {
  s << b.center() << std::endl;;
  s << b.radius() << std::endl;;
  return s;
}
    
template <size_t N_dimension, class T>
std::istream& operator >>(std::istream& s, sl::ball<N_dimension,T>& b) {
  s >> b.center();
  s >> b.radius();
  return s;
}

namespace sl {

  /// An object for building bounding balls from point/ball clouds
  template <size_t N_dimension, class T> 
  class ball_builder: 
    public bounding_volume_builder<ball<N_dimension,T> > {
  public:

    typedef bounding_volume_builder<ball<N_dimension,T> > super_t;
    typedef ball_builder<N_dimension, T> this_t;

    enum { dimension = super_t::dimension };

    typedef typename super_t::value_t           value_t;
    typedef typename super_t::point_t           point_t;
    typedef typename point_t::vector_t          vector_t;
    typedef typename super_t::bounding_volume_t bounding_volume_t;
    typedef bounding_volume_t                   ball_t;
    
    typedef std::list<point_t>                     point_list_t;
    typedef typename point_list_t::iterator        point_list_iterator_t;
    typedef typename point_list_t::const_iterator  point_list_const_iterator_t;
    
  public: // Constraints

    SL_COMPILE_TIME_CHECK("Non null dimension", dimension > 0);

  protected: // Minball support

    class mbb_basis {
    protected:
      const std::vector<ball_t>* balls_;

      ball_t                     solution_;
      std::vector<std::size_t>   support_;
    public:

      mbb_basis(const std::vector<ball_t>* balls)
        : balls_(balls) {
      }

      ~mbb_basis() {
      }

      std::size_t support_count() const {
        return support_.size();
      }

      std::size_t support(std::size_t i) const {
        SL_REQUIRE("good index", i<support_count());
        return support_[i];
      }

      const ball_t& solution() const {
        return solution_;
      }
      
      bool is_violated_by(std::size_t i) const {
        return
          (support_count() == 0) ||
          (!solution().epsilon_contains((*balls_)[i]));
      }

      void expand(std::size_t new_constraint) {
        SL_REQUIRE("violated", is_violated_by(new_constraint));
               
        if (support_count() == 0) {
          support_.push_back(new_constraint);
          solution_ = (*balls_)[new_constraint];
        } else {
          // Generate all possible subsets containing
          // existing support + new ball and look for
          // best solution
          
          std::vector<std::size_t> best_support;
          ball_t                   best_ball;
          best_ball.to_empty();
          
          std::size_t N = support_count();
          std::size_t subset_count = 0;
          for(std::size_t i=0;i<N;i++) subset_count|=1<<i;
          std::vector<std::size_t> test_support;
          std::vector<std::size_t> test_others;
          for(std::size_t k=0;k<=subset_count;k++) {
            test_support.clear();
            test_others.clear();
            for(std::size_t i=0;i<N;i++) {
              if(k&(1<<i)) {
                test_support.push_back(support(i));
              } else {
                test_others.push_back(support(i));
              }
            }
            test_support.push_back(new_constraint);
            if (test_support.size() <= dimension+1) {
              ball_t test_ball = restricted_mbb(test_support);
              if (test_ball.is_empty()) {
                // Linearly dependent basis detected
              } else if (best_ball.is_empty() || (test_ball.radius() < best_ball.radius())) {
                if (check_containment(test_ball, test_others)) {
                  // Update best solution
                  best_support = test_support;
                  best_ball = test_ball;
                }
              }
            }
          }

          // FIXME: Should check whether because of numerical errors no
          // solution was found
          support_  = best_support;
          solution_ = best_ball;
        }
      }

      bool check_containment(const ball_t& b,
                             const std::vector<std::size_t>& others) const {
        bool result = true;
        for (std::size_t i=0; result && (i<others.size()); ++i) {
          result = b.epsilon_contains((*balls_)[others[i]]);
        }
        return result;                              
      }
      
      ball_t restricted_mbb(const std::vector<std::size_t>& support) {
        SL_REQUIRE("Good dim", support.size() <= dimension+1);
        
        typedef fixed_size_square_matrix<dimension,value_t>             matrix_t;
        typedef fixed_size_vector<column_orientation,dimension,value_t> vector_t;
         
        ball_t result;
        result.to_empty();
        
        std::size_t N = support.size();
        switch (N) {
        case 0:
          break;
        case 1:
          result = (*balls_)[support[0]];
          break;
#if 0
        case 2:
          if ((*balls_)[support[0]].contains((*balls_)[support[1]]) ||
              (*balls_)[support[1]].contains((*balls_)[support[0]])) {
            // Invalid support
          } else {
            result = (*balls_)[support[0]];
            result.merge((*balls_)[support[1]]);
          }
          break;
#endif
        default:
          // General case
          value_t r0 = (*balls_)[support[0]].radius();
          const point_t& c0 = (*balls_)[support[0]].center();
          
          vector_t d[dimension];
          for (std::size_t i=0; i<N-1; ++i) {
            d[i] = (*balls_)[support[i+1]].center() - c0;
          }
          matrix_t Q;
          Q.to_identity();
          for (std::size_t i=0; i<N-1; ++i) {
            for (std::size_t j=0; j<N-1; ++j) {
              Q(i,j) = d[i].dot(d[j]);
            }
          }
          matrix_t Qinv; bool invertible;
          Q.invert_to(Qinv, &invertible);
          if (!invertible) {
            // Invalid support: linearly dependent centers
          } else {
            // Compute ball touching all elements
            const value_t One  = scalar_math<value_t>::one();
            const value_t Two  = scalar_math<value_t>::two();
            const value_t Four = Two*Two;
            const value_t One_Plus_Eps = One+value_t(0.001); // FIXME
            vector_t r;
            vector_t b;
            value_t rmax = r0;
            
            for (std::size_t i=0; i<N-1; ++i) {
              value_t ri = (*balls_)[support[i+1]].radius();
              rmax = sl::max(rmax, ri);
              r[i] = r0 - ri;
              b[i] = (d[i].dot(d[i]) + r0*r0 - ri*ri)/Two;
            }
            rmax /= One_Plus_Eps;
            
            value_t A = r.dot(Qinv*r)-One;
            value_t B = -Two*(r.dot(Qinv*b)-r0);
            value_t C = b.dot(Qinv*b)-r0*r0;

            // Find radius
            value_t rho_star = -One;
            if (sl::is_zero(A)) {
              if (sl::is_zero(B)) {
                // Infeasible
              } else {
                // Valid if positive
                value_t rho = -C/B;
                if (rho >= rmax) rho_star = rho; 
              }
            } else {
              value_t D = B*B-Four*A*C;
              if (sl::is_negative(D)) {
                // Infeasible
              } else {
                value_t sqrtD = std::sqrt(D);
                value_t rho1 = (-B + sqrtD)/(Two*A);
                value_t rho2 = (-B - sqrtD)/(Two*A);
                if (rho1 < rmax) {
                  if (rho2 < rmax) {
                    // Invalid
                  } else {
                    rho_star = rho2;
                  }
                } else if (rho2 < rmax) {
                  SL_REQUIRE("Good rho1", rho1 >= rmax);
                  rho_star = rho1;
                } else {
                  SL_REQUIRE("Good rho1", rho1 >= rmax);
                  SL_REQUIRE("Good rho2", rho2 >= rmax);
              
                  rho_star = sl::min(rho1, rho2);
                }
              }
            }

            // Find center
            if (sl::is_non_negative(rho_star)) {            
              vector_t l_star = Qinv*(b - rho_star*r);
              value_t l_star_0 = One;
              point_t c_star;
              for (std::size_t i=0; i<N-1; ++i) {
                const point_t& ci = (*balls_)[support[i+1]].center();
                l_star_0 -= l_star[i]; 
                c_star   += l_star[i]*ci.as_vector();
              }
              c_star += l_star_0*c0.as_vector();

              // Assign result
              result = ball_t(c_star, rho_star);
            }
          }
          break;
        }

        return result;
      }
    };

    
    class mbp_basis {
    protected:
      // data members
      int                  m, s;   // size and number of support points
      point_t              q0;
      
      value_t              z[dimension+1];
      value_t              f[dimension+1];
      vector_t             v[dimension+1];
      value_t              a[dimension+1][dimension];
      
      point_t              c[dimension+1];
      value_t              sqr_r[dimension+1];
      
      point_t*             current_c;      // points to some c[j]
      value_t              current_sqr_r;
      
    public:
      
      mbp_basis() {
        reset();
      }
      
      // access
      const point_t& center() const {
        return *current_c;
      }
      
      value_t squared_radius() const {
        return current_sqr_r;
      }
      
      int size() const {
        return m;
      }
      
      int support_size() const {
        return s;
      }
      
      value_t excess(const point_t& p) const {
        return -current_sqr_r + p.distance_squared_to(center());
      }
      
    public: // modification
      
      // generates empty sphere with m=s=0
      void reset() {
        m = s = 0;
        // we misuse c[0] for the center of the empty sphere
        c[0] = point_t();
        current_c = &(c[0]);
        current_sqr_r = -1;
      }
      
      bool push (const point_t& p) {
        int i, j;
        value_t eps = scalar_math<value_t>::epsilon();
        if (m==0) {
          q0   = p;
          c[0] = q0;
          sqr_r[0] = 0;
        } else {
          // set v_m to Q_m
          for (i=0; i<dimension; ++i)
            v[m][i] = p[i]-q0[i];
          
          // compute the a_{m,i}, i< m
          for (i=1; i<m; ++i) {
            a[m][i] = 0;
            for (j=0; j<dimension; ++j)
              a[m][i] += v[i][j] * v[m][j];
            a[m][i]*=(scalar_math<value_t>::two()/z[i]);
          }
          
          // update v_m to Q_m-\bar{Q}_m
          for (i=1; i<m; ++i) {
            for (j=0; j<dimension; ++j)
              v[m][j] -= a[m][i]*v[i][j];
          }
          
          // compute z_m
          z[m]=0;
          for (j=0; j<dimension; ++j)
            z[m] += sqr(v[m][j]);
          z[m]*=scalar_math<value_t>::two();
          
          // reject push if z_m too small
          if (z[m]<eps*current_sqr_r) {
            return false;
          }
          
          // update c, sqr_r
          value_t e = p.distance_squared_to(c[m-1])-sqr_r[m-1];
          f[m]=e/z[m];
          
          for (i=0; i<dimension; ++i)
            c[m][i] = c[m-1][i]+f[m]*v[m][i];
          sqr_r[m] = sqr_r[m-1] + e*f[m]/2;
        }
        current_c = &(c[m]);
        current_sqr_r = sqr_r[m];
        s = ++m;
        return true;
      }
      
      void pop() {
        --m;
      }
      
    public: // checking
      value_t slack() const {
        value_t l[dimension+1], min_l=0;
        l[0] = 1;
        for (int i=s-1; i>0; --i) {
          l[i] = f[i];
          for (int k=s-1; k>i; --k)
            l[i]-=a[k][i]*l[k];
          if (l[i] < min_l) min_l = l[i];
          l[0] -= l[i];
        }
        if (l[0] < min_l) min_l = l[0];
        return ( (min_l < 0) ? -min_l : 0);
      }
    };

  protected: // Data

    point_list_t points_;             //< the list of points
    std::vector<ball_t> balls_; //< the list of spheres
    
    mbp_basis             mbp_basis_; //< basis keeping the current ball
    point_list_iterator_t mbp_support_end_;   //< past-the-end iterator of support set
    
  public: // Creation & Destruction

    /// Default init
    inline ball_builder() {
    }

    inline virtual ~ball_builder() {
    }
  
  public: // Construction

    /// Start point/ball cloud
    virtual inline void begin_model() {
      super_t::begin_model();
      points_.clear();
      balls_.clear();
    }
    
    /// Add point to current cloud
    virtual inline void put_point(const point_t& p) {
      if (balls_.empty()) {
        points_.push_back(p);
      } else {
        balls_.push_back(ball_t(p));
      }
    }

    /// Add ball to current cloud
    virtual inline void put_ball(const ball_t& b) {
      if (!b.is_empty()) {
        if (b.is_point()) {
          put_point(b.center());
        } else {
          while (!points_.empty()) {
            balls_.push_back(ball_t(*points_.begin()));
            points_.pop_back();
          }
          balls_.push_back(b);
        }
      }
    }

    /// Finish point/ball cloud and build ball
    virtual inline void end_model() {
      construct();
      points_.clear();
      balls_.clear();
      super_t::end_model();
    }

  protected:

    virtual void construct() {
      (this->last_bounding_volume_).to_empty();
      if (!points_.empty()) {
        SL_REQUIRE("No balls", balls_.empty());
        construct_from_points();
      } else if (!balls_.empty()) {
        SL_REQUIRE("No points", points_.empty());
        construct_from_balls();

        // FIXME FIXME FIXME
        // The following code is a fallback for cases
        // where numerical instabilities break the
        // minball of balls code!! This is
        // a temporary hack, that should be removed ASAP
        // by fixing the restricted_mbb code
        
        bool ok = true;
        std::size_t N = balls_.size();
        for (std::size_t i=0; i<N && ok; ++i) {
          ok = (this->last_bounding_volume_).epsilon_contains(balls_[i], value_t(0.01));
        }
        
        ball_t b2 = quick_and_dirty_stable_miniball_of_balls();
        if (!ok || b2.radius() < (this->last_bounding_volume_).radius()) {
          SL_TRACE_OUT(0) << "Fixing miniball code: in: " << ok << "r: " << (this->last_bounding_volume_).radius() << " => " << b2.radius() << std::endl;
          (this->last_bounding_volume_) = b2;
        }
      }
    }

  protected: // Balls of Balls
    
    // Quick'n'dirty approximation algorithm
    ball_t quick_and_dirty_stable_miniball_of_balls() const {
      ball_t result;
      result.to_empty();
      
      std::size_t N = balls_.size();
      switch (N) {
      case 0: result.to_empty(); break;
      case 1: result = balls_[0]; break;
      case 2: result = balls_[0]; result.merge(balls_[1]); break;
      default:
        // Pick initial approximation
        std::size_t m=0;
        std::size_t r=0;
        value_t dmr = balls_[m].radius();
        for (std::size_t i=0; i<N; ++i) {
          value_t dmi = balls_[m].center().distance_to(balls_[i].center())+balls_[i].radius();
          if (dmi>dmr) { r = i; dmr = dmi; }
        }
        std::size_t l=r;
        value_t dlr = balls_[l].radius();
        for (std::size_t i=0; i<N; ++i) {
          value_t dri = balls_[r].center().distance_to(balls_[i].center())+balls_[i].radius();
          if (dri>dlr) { r = i; dlr = dri; }
        }
        result = balls_[l];
        result.merge(balls_[r]);

        // Grow to include all other balls
        for (std::size_t i=0; i<N; ++i) {
          if (i!= l && i!= r) {
            result.merge(balls_[i]);
          }
        }
      }

      return result;
    }

#if 0
    void mbb_lp(random::irng_marsaglia&   irng,
                std::vector<std::size_t>& untested,
                mbb_basis&                basis) {
      std::size_t N = untested.size();
      if (N == 0) {
        // done, current basis meets all constraints
      } else {
        // Remove random constraint
        std::size_t picked  = irng.value() % N;
        std::size_t removed = untested[picked];
        untested[picked] = untested[N-1];
        untested.pop_back();
        
        // Build new basis for smaller constraint set
        mbb_basis basis_prime = basis;
        mbb_lp(irng, untested, basis_prime);
        
        if (basis_prime.is_violated_by(removed)) {
          // Removed is part of the basis, add it
          basis = basis_prime;
          basis.expand(removed);
          mbb_lp(irng,untested,basis);
        } else {
          // Accept lower level solution
          basis = basis_prime;
        }

        untested.push_back(removed);
      }
    }
#else
    void mbb_lp(random::irng_marsaglia&   irng,
                std::vector<std::size_t>& untested,
                mbb_basis&                basis) {
      std::size_t N = untested.size();
      if (N>0) {
        // Remove random constraint
        std::size_t picked  = irng.value() % N;
        std::size_t removed = untested[picked];
        untested[picked] = untested[N-1];
        untested.pop_back();
        
        // Build new basis for smaller constraint set
        mbb_lp(irng, untested, basis);
        
        if (basis.is_violated_by(removed)) {
          // Removed is part of the basis, add it
          basis.expand(removed);
          mbb_lp(irng,untested,basis);
        } 
        untested.push_back(removed);
      }
    }
#endif
    
    virtual void construct_from_balls() {
      std::size_t N = balls_.size();
      switch (N) {
      case 0: (this->last_bounding_volume_).to_empty(); break;
      case 1: (this->last_bounding_volume_) = balls_[0]; break;
      case 2: (this->last_bounding_volume_) = balls_[0]; (this->last_bounding_volume_).merge(balls_[1]); break;
      default:
        std::vector<std::size_t> untested;
        for (std::size_t i=0; i<N; ++i) {
          untested.push_back(i);
        }
        random::irng_marsaglia   irng;
        mbb_basis                basis(&balls_);
        mbb_lp(irng, untested, basis);
        (this->last_bounding_volume_) = basis.solution();
      }
    }
    
  protected: // Balls of points

    // Exact algorithm, based on Bernd Gaertner's miniball
    virtual void construct_from_points() {
      // Construct from point samples
      mbp_basis_.reset();
      mbp_support_end_ = points_.begin();
      pivot_miniball (points_.end());
      (this->last_bounding_volume_) = ball_t(mbp_basis_.center(),
                                     std::sqrt(mbp_basis_.squared_radius()));
    }

    void move_to_front_miniball(point_list_iterator_t i) {
      mbp_support_end_ = points_.begin();
      if ((mbp_basis_.size())==dimension+1) return;
      for (point_list_iterator_t k=points_.begin(); k!=i;) {
        point_list_iterator_t j=k++;
        if (is_positive(mbp_basis_.excess(*j))) {
          if (mbp_basis_.push(*j)) {
            move_to_front_miniball (j);
            mbp_basis_.pop();
            move_to_front(j);
          }
        }
      }
    }

    void pivot_miniball(point_list_iterator_t i) {
      point_list_iterator_t t = ++points_.begin();
      move_to_front_miniball(t);
      value_t max_e, old_sqr_r = 0.0f;
      do {
        point_list_iterator_t pivot;
        max_e = max_excess (t, i, pivot);
        if (is_positive(max_e)) {
          t = mbp_support_end_;
          if (t==pivot) ++t;
          old_sqr_r = mbp_basis_.squared_radius();
          mbp_basis_.push (*pivot);
          move_to_front_miniball (mbp_support_end_);
          mbp_basis_.pop();
          move_to_front (pivot);
        }
      } while (is_positive(max_e) && (mbp_basis_.squared_radius() > old_sqr_r));
    }
    
    void move_to_front(point_list_iterator_t j) {
      if (mbp_support_end_ == j)
        ++mbp_support_end_;
      points_.splice (points_.begin(), points_, j);
    }
    
    value_t max_excess (point_list_iterator_t t, point_list_iterator_t i, point_list_iterator_t& pivot) const {
      const point_t& c = mbp_basis_.center();
      const value_t  sqr_r = mbp_basis_.squared_radius();
      value_t e, max_e = 0;
      for (point_list_iterator_t k=t; k!=i; ++k) {
        const value_t *p = (*k).begin();
        e = -sqr_r;
        for (int j=0; j<dimension; ++j)
          e += sqr(p[j]-c[j]);
        if (e > max_e) {
          max_e = e;
          pivot = k;
        }
      }
      return max_e;
    }

  }; // class ball_builder

}; // namespace sl

// ---------------------------------------------------------------------------
// Typedefs for common cases
// ---------------------------------------------------------------------------

namespace sl {
  // A 2D ball with single precision floating point components
  typedef ball<2,float> ball2f;
  // A 3D ball with single precision floating point components
  typedef ball<3,float> ball3f;
  // A 4D ball with single precision floating point components
  typedef ball<4,float> ball4f;

  // A 2D ball with double precision floating point components
  typedef ball<2,double> ball2d;
  // A 3D ball with double precision floating point components
  typedef ball<3,double> ball3d;
  // A 4D ball with double precision floating point components
  typedef ball<4,double> ball4d;
}

#endif
