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
#ifndef SL_INDEXED_CONTAINER_HPP
#define SL_INDEXED_CONTAINER_HPP

#include <sl/indexed.hpp>
#include <sl/indexed_manifest_initializer.hpp>

namespace sl {

  /**
   *  N-dimensional assignable indexed entities
   */
  template <
    class  G_numtype, 
    size_t G_rank, 
    class  G_derived, 
    class  G_userdefined,
    class  G_discriminant
  >
  class indexed_container: 
    public indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant > {

  public:
    
    typedef indexed<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant>           super_t;
    typedef indexed_container<G_numtype, G_rank, G_derived, G_userdefined, G_discriminant> this_t;
    typedef typename super_t::derived_t                                                    derived_t;
    typedef typename super_t::user_defined_data_t                                          user_defined_data_t;
    typedef typename super_t::value_t                                                      value_t;
    typedef typename super_t::discriminant_t                                               discriminant_t;
    typedef typename super_t::subscript_t                                                  subscript_t;

    /// Consider all dimensions fixed inside an expression
    enum { xpr_dynamic_dimensions = false };

    /// The rank of the object
    enum { rank_c = super_t::rank_c };

  public:

    /// Default init (no operation)
    explicit inline indexed_container() {
    }
    
    /// Data init
    explicit inline indexed_container(const user_defined_data_t& data): super_t(data) {
    }
    
  public: // Dynamic expression bounds

    /**
     *  The upper bound for rank r.
     *  This is extent(r) if r has a fixed dimension,
     *  XPR_DYNAMIC_BOUND otherwise.
     */
    size_t xpr_ubound(size_t r) const {
      return r >= rank_c ? 0 : (this->derived_ref()).extent()(r);
    }

  public: // Access
    
    /// The element referenced by subscript idx
    inline value_t operator()(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->derived_ref()).item(idx);
    }

    /// The element referenced by subscript idx (alias of operator())
    inline value_t operator[](const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->derived_ref()).item(idx);
    }

  public: // Update

    /// Set element at index idx to zero
    //  (SUBCLASS RESPONSIBILITY)
    inline void put_zero(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->derived_ref()).put_zero(idx);
    }

    /// Set element at index idx to the value v, different than zero
    //  (SUBCLASS RESPONSIBILITY)
    inline void put_non_zero(const value_t& v, const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->derived_ref()).put_non_zero(v, idx);
    }

    /// Set element at index idx to the value v
    //  (SUBCLASS RESPONSIBILITY)
    inline void put(const value_t& v, const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->derived_ref()).put(v, idx);
    }

    /// Write access to element at subscript idx (alias of operator())
    // (SUBCLASS RESPONSIBILITY)
    inline value_t& item_ref(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->derived_ref()).item_ref(idx);
    }

    /// Write access to element at subscript idx
    inline value_t& operator()(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->derived_ref()).item_ref(idx);
    }

    /// Write access to element at subscript idx (alias of operator())
    inline value_t& operator[](const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->derived_ref()).item_ref(idx);
    }

    /// Is current resizable to size sz? 
    // (SUBCLASS RESPONSIBILITY)
    inline bool resizable(const subscript_t& sz) const {
      return 
	sz == (this->extent()) || 
	(this->derived_ref()).resizable(sz);
    }	

    /// If resizable, make this size equal to sz
    // (SUBCLASS RESPONSIBILITY)
    inline void resize(const subscript_t& sz) {
      (this->derived_ref()).resize(sz);
      SL_ENSURE("Good size", (this->extent()) == sz);
    }	

    /// If resizable, make this size equal to sz
    inline void resize(size_t i0) {
      resize(subscript_t(i0));
    }
    /// If resizable, make this size equal to sz
    inline void resize(size_t i0, size_t i1) {
      resize(subscript_t(i0,i1));
    }
    /// If resizable, make this size equal to sz
    inline void resize(size_t i0, size_t i1, size_t i2) {
      resize(subscript_t(i0,i1,i2));
    }
    /// If resizable, make this size equal to sz
    inline void resize(size_t i0, size_t i1, size_t i2, size_t i3) {
      resize(subscript_t(i0,i1,i2,i3));
    }
    /// If resizable, make this size equal to sz
    inline void resize(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) {
      resize(subscript_t(i0,i1,i2,i3,i4));
    }
                       
    /// Set everything to zero 
    // (SUBCLASS RESPONSIBILITY)
    inline void clear() {
      (this->derived_ref()).clear();
    }
    
  public: // Subscript remapping

    //    SL_DECLARE_SUBSCRIPT_REMAPPING(derived_t);

  public: // Helpers for indexing

    inline value_t operator()(size_t i0) const {
      return operator()(subscript_t(i0));
    }

    inline value_t operator()(size_t i0, size_t i1) const {
      return operator()(subscript_t(i0, i1));
    }

    inline value_t operator()(size_t i0, size_t i1, size_t i2) const {
      return operator()(subscript_t(i0, i1, i2));
    }

    inline value_t operator()(size_t i0, size_t i1, size_t i2, size_t i3) const {
      return operator()(subscript_t(i0, i1, i2, i3));
    }

    inline value_t& operator()(size_t i0) {
      return operator()(subscript_t(i0));
    }

    inline value_t& operator()(size_t i0, size_t i1) {
      return operator()(subscript_t(i0, i1));
    }

    inline value_t& operator()(size_t i0, size_t i1, size_t i2) {
      return operator()(subscript_t(i0, i1, i2));
    }

    inline value_t& operator()(size_t i0, size_t i1, size_t i2, size_t i3) {
      return operator()(subscript_t(i0, i1, i2, i3));
    }

  public: // Assignment
    
    /** 
     *  Assign from expression rhs (actual computation is performed here!).
     *  All other assignement routines should make use of this one,
     *  eventually converting their parameters to a suitable indexed.
     *  At the minimum, subclasses should define operator= for 
     *  indexed expressions and containers.
     */
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline derived_t& assign_from(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs);
    
  public: // Arithmetic

#if 0
    ////////////////////// REMOVED - NOT CURRENTLY IMPLEMENTED...

    /// Elementwise sum
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline derived_t& operator+=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
      SL_REQUIRE("Good size", length() == rhs.length());
      return assign_from((this->derived_ref()) + rhs);
    }
    
    /// Elementwise sum
    inline derived_t& operator+=(const G_numtype& rhs) {
      return assign_from((this->derived_ref()) + rhs);
    }
    
    /// Elementwise difference
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline derived_t& operator-=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
      SL_REQUIRE("Good size", length() == rhs.length());
      return assign_from((this->derived_ref()) - rhs);
    }
    
    /// Elementwise difference
    inline derived_t& operator-=(const G_numtype& rhs) {
      return assign_from((this->derived_ref()) - rhs);
    }
    
    /// Elementwise multiplication 
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline derived_t& operator*=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
      SL_REQUIRE("Good size", length() == rhs.length());
      return assign_from((this->derived_ref()) * rhs);
    }
    
    /// Elementwise multiplication
    inline derived_t& operator*=(const G_numtype& rhs) {
      return assign_from((this->derived_ref()) * rhs);
    }
    
    /// Elementwise division
    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline derived_t& operator/=(const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
      SL_REQUIRE("Good size", length() == rhs.length());
      return assign_from((this->derived_ref()) / rhs);
    }
    
    /// Elementwise division
    inline derived_t& operator/=(const G_numtype& rhs) {
      return assign_from((this->derived_ref()) / rhs);
    }
#endif
    
  };

}

//---------------------------------------------------------------------------
// Assignment

#include <sl/indexed_assign.hpp>

#endif

