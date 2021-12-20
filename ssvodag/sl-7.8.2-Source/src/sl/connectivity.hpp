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
#ifndef SL_CONNECTIVITY_HPP
#define SL_CONNECTIVITY_HPP

#include <sl/operators.hpp>
#include <sl/hash.hpp>
#include <cassert>

namespace sl {

  // -----------------------------------------------------------------------
  /// Directed edge connectivity
  class directed_edge_connectivity {
  public:
    typedef uint64_t                   vertex_t;
    typedef directed_edge_connectivity this_t;
  protected:
    vertex_t vidx_[2];
  public:
    inline directed_edge_connectivity() {}
    
    inline directed_edge_connectivity(vertex_t v0, vertex_t v1) {
      vidx_[0] = v0;
      vidx_[1] = v1;
    }

    inline vertex_t operator[] (std::size_t i) const {
      SL_REQUIRE("Good index", i<2);
      return vidx_[i];
    }
    
    inline bool operator< (const this_t& other) const {
      if ((*this)[0] < other[0]) return true;
      if ((*this)[0] > other[0]) return false;
      if ((*this)[1] < other[1]) return true;
      //if ((*this)[1] > other[1]) return false;
      return false;
    }
    
    inline bool operator== (const this_t& other) const {
      return
        ((*this)[0] == other[0]) &&
        ((*this)[1] == other[1]);
    }

    static inline vertex_t index_end() {
      return vertex_t(-1);
    }
    
    inline vertex_t index_of(vertex_t vidx) const {
      vertex_t result = index_end();
      if (vidx_[1] == vidx) result = 1;
      if (vidx_[0] == vidx) result = 0;
      return result;
    }

    inline bool has(vertex_t i) const {
      return index_of(i) != index_end();
    }

    inline this_t remapped(vertex_t old_vtx, vertex_t new_vtx) const {
      return this_t(vidx_[0] == old_vtx ? new_vtx : vidx_[0],
                    vidx_[1] == old_vtx ? new_vtx : vidx_[1]);
    }
    
    inline void remap(vertex_t old_vtx, vertex_t new_vtx) {
      (*this) = remapped(old_vtx, new_vtx);
    }
    
    inline vertex_t opposite(vertex_t vidx) const {
      SL_REQUIRE("Exists", has(vidx));
      
      if (vidx_[0] == vidx) {
        return vidx_[1];
      } else {
        SL_REQUIRE("Exists", vidx_[1] == vidx);
        return vidx_[0];
      }
    }
    
    inline bool is_valid() const {
      return
        (vidx_[0] != vidx_[1]);
    }

    SL_OP_COMPARABLE1(this_t);
    SL_OP_EQUALITY_COMPARABLE1(this_t);
  };

  inline std::size_t hash_value(const directed_edge_connectivity& x) {
    std::size_t seed = 0;
    hash_combine(seed, sl::uint64_t(x[0]));
    hash_combine(seed, sl::uint64_t(x[1]));
    return seed;
  }

  SL_HASH_SPECIALIZE_REF(directed_edge_connectivity)
}

namespace sl {
  // -----------------------------------------------------------------------
  /// Undirected edge connectivity
  class undirected_edge_connectivity {
  public:
    typedef uint64_t                     vertex_t;
    typedef undirected_edge_connectivity this_t;
  protected:
    vertex_t vidx_[2];
  public:
    inline undirected_edge_connectivity() {}
    
    inline undirected_edge_connectivity(vertex_t v0, vertex_t v1) {
      if (v0<v1) {
        vidx_[0] = v0;
        vidx_[1] = v1;
      } else {
        vidx_[0] = v1;
        vidx_[1] = v0;
      }
    }

    inline vertex_t operator[] (std::size_t i) const {
      SL_REQUIRE("Good index", i<2);
      return vidx_[i];
    }
    
    inline bool operator< (const this_t& other) const {
      if ((*this)[0] < other[0]) return true;
      if ((*this)[0] > other[0]) return false;
      if ((*this)[1] < other[1]) return true;
      //if ((*this)[1] > other[1]) return false;
      return false;
    }
    
    inline bool operator== (const this_t& other) const {
      return
        ((*this)[0] == other[0]) &&
        ((*this)[1] == other[1]);
    }

    static inline vertex_t index_end() {
      return vertex_t(-1);
    }
    
    inline vertex_t index_of(vertex_t vidx) const {
      vertex_t result = index_end();
      if (vidx_[0] == vidx) result = 0;
      if (vidx_[1] == vidx) result = 1;
      return result;
    }

    inline bool has(vertex_t i) const {
      return index_of(i) != index_end();
    }

    inline this_t remapped(vertex_t old_vtx, vertex_t new_vtx) const {
      return this_t(vidx_[0] == old_vtx ? new_vtx : vidx_[0],
                    vidx_[1] == old_vtx ? new_vtx : vidx_[1]);
    }
    
    inline void remap(vertex_t old_vtx, vertex_t new_vtx) {
      (*this) = remapped(old_vtx, new_vtx);
    }
    
    inline vertex_t opposite(vertex_t vidx) const {
      SL_REQUIRE("Exists", has(vidx));
      if (vidx_[0] == vidx) {
        return vidx_[1];
      } else {
        SL_REQUIRE("Exists", vidx_[1] == vidx);
        return vidx_[0];
      }
    }
    
    inline bool is_valid() const {
      return
        (vidx_[0] != vidx_[1]);
    }

    SL_OP_COMPARABLE1(this_t);
    SL_OP_EQUALITY_COMPARABLE1(this_t);
  };
  
  inline std::size_t hash_value(const undirected_edge_connectivity& x) {
    std::size_t seed = 0;
    hash_combine(seed, sl::uint64_t(x[0]));
    hash_combine(seed, sl::uint64_t(x[1]));
    return seed;
  }

  SL_HASH_SPECIALIZE_REF(undirected_edge_connectivity)
}

namespace sl {
  // -----------------------------------------------------------------------
  /// Triangle connectivity
  class triangle_connectivity {
  public:
    typedef uint64_t                     vertex_t;
    typedef directed_edge_connectivity   half_edge_t;
    typedef undirected_edge_connectivity edge_t;
    typedef triangle_connectivity        this_t;
  protected:
    vertex_t vidx_[3];
  public:
    inline triangle_connectivity() {}
    
    inline triangle_connectivity(vertex_t v0, vertex_t v1, vertex_t v2) {
      if (v0<v1) {
        if (v0<v2) {
          vidx_[0] = v0;
          vidx_[1] = v1;
          vidx_[2] = v2;
        } else {
          vidx_[0] = v2;
          vidx_[1] = v0;
          vidx_[2] = v1;
        }          
      } else {
        if (v1<v2) {
          vidx_[0] = v1;
          vidx_[1] = v2;
          vidx_[2] = v0;
        } else {
          vidx_[0] = v2;
          vidx_[1] = v0;
          vidx_[2] = v1;
        }
      }
    }

    inline vertex_t operator[] (std::size_t i) const {
      SL_REQUIRE("Good index", i<3);
      return vidx_[i];
    }

    inline edge_t undirected_edge(std::size_t i) const {
      SL_REQUIRE("Good index", i<3);
      return edge_t(vidx_[i], vidx_[(i+1)%3]);
    }

    inline half_edge_t directed_edge(std::size_t i) const {
      SL_REQUIRE("Good index", i<3);
      return half_edge_t(vidx_[i], vidx_[(i+1)%3]);
    }
    
    inline bool operator< (const this_t& other) const {
      if ((*this)[0] < other[0]) return true;
      if ((*this)[0] > other[0]) return false;
      if ((*this)[1] < other[1]) return true;
      if ((*this)[1] > other[1]) return false;
      if ((*this)[2] < other[2]) return true;
      //if ((*this)[2] > other[2]) return false;
      return false;
    }
    
    inline bool operator== (const this_t& other) const {
      return
        ((*this)[0] == other[0]) &&
        ((*this)[1] == other[1]) &&
        ((*this)[2] == other[2]);
    }

    inline this_t remapped(vertex_t old_vtx, vertex_t new_vtx) const {
      return this_t(vidx_[0] == old_vtx ? new_vtx : vidx_[0],
                    vidx_[1] == old_vtx ? new_vtx : vidx_[1],
                    vidx_[2] == old_vtx ? new_vtx : vidx_[2]);
    }
    
    inline void remap(vertex_t old_vtx, vertex_t new_vtx) {
      (*this) = remapped(old_vtx, new_vtx);
    }

    inline bool is_valid() const {
      return
        (vidx_[0] != vidx_[1]) &&
        (vidx_[0] != vidx_[2]) &&
        (vidx_[1] != vidx_[2]);
    }

    static inline vertex_t index_end() {
      return vertex_t(-1);
    }
    
    inline std::size_t index_of(vertex_t vidx) const {
      vertex_t result = index_end();
      if (vidx_[0] == vidx) result = 0;
      if (vidx_[1] == vidx) result = 1;
      if (vidx_[2] == vidx) result = 2;
      return result;
    }

    inline bool has(vertex_t i) const {
      return index_of(i) != index_end();
    }

    inline std::size_t index_of(const half_edge_t& e) const {
      std::size_t result = index_end();
      if (vidx_[0] == e[0] && vidx_[1] == e[1]) result = 0;
      if (vidx_[1] == e[0] && vidx_[2] == e[1]) result = 1;
      if (vidx_[2] == e[0] && vidx_[0] == e[1]) result = 2;
      return result;
    }

    inline bool has(const half_edge_t& e) const {
      return index_of(e) != index_end();
    }

    inline std::size_t index_of(const edge_t& e) const {
      std::size_t result = index_of(half_edge_t(e[0], e[1]));
      if (result == index_end()) result = index_of(half_edge_t(e[1], e[0]));
      return result;
    }

    inline bool has(const edge_t& e) const {
      return index_of(e) != index_end();
    }

    inline vertex_t opposite(const edge_t& e) const {
      SL_REQUIRE("exists", has(e));
      std::size_t i = index_of(e);
      return vidx_[(i+2)%3];
    }
    
    SL_OP_COMPARABLE1(this_t);
    SL_OP_EQUALITY_COMPARABLE1(this_t);
  };

  inline std::size_t hash_value(const triangle_connectivity& x) {
    std::size_t seed = 0;
    hash_combine(seed, sl::uint64_t(x[0]));
    hash_combine(seed, sl::uint64_t(x[1]));
    hash_combine(seed, sl::uint64_t(x[2]));
    return seed;
  }

  SL_HASH_SPECIALIZE_REF(triangle_connectivity)
}

#endif
