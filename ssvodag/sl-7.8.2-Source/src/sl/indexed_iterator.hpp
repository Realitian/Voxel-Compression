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
#ifndef SL_INDEXED_ITERATOR_HPP
#define SL_INDEXED_ITERATOR_HPP

#include <sl/index.hpp>
#include <sl/utility.hpp>

namespace sl {

  //----------------------------------------------------------
  /**
   *  Helper for defining forward iterators that keep
   *  track of the current index.
   *  Subclasses should define all features marked
   *  SUBCLASS RESPONSIBILITY.
   */
  template <
    class  G_numtype, 
    size_t G_rank, 
    class  G_derived, 
    class  G_userdefined
  >
  class indexed_const_iterator_base {
  public:

    typedef indexed_const_iterator_base<G_numtype, G_rank, G_derived, G_userdefined> this_t;
    typedef G_derived                                                                derived_t;
    typedef G_userdefined                                                            user_defined_data_t;
    typedef G_numtype                                                                value_t;

    typedef ::sl::index<G_rank>                                                      subscript_t;

    enum { rank_c = G_rank, max_rank = G_rank-1 };

  protected:

    /// User defined data. Declared here to remove the space overhead associated to empty base classes.
    SL_DECLARE_GENERIC_SUPERCLASS_DATA(user_defined_data_t,data_);

  protected:

    /// Data init
    inline explicit indexed_const_iterator_base(const user_defined_data_t& data): data_(data) {
    }

    /// Default init
    inline explicit indexed_const_iterator_base() {
    }

  public: // Cast features

    /// Features to access this as a derived_t pointer
    SL_DECLARE_GENERIC_SUPERCLASS_FEATURES(derived_t);

  public: // Subclass responsibility features

    /// Element access 
    // (SUBCLASS RESPONSIBILITY)
    inline value_t value() const {
      return (this->derived_ref()).value();
    }

    /// Current index 
    // (SUBCLASS RESPONSIBILITY)
    inline subscript_t index() const {
      return (this->derived_ref()).index();
    }

    /// Pre-increment 
    // (SUBCLASS RESPONSIBILITY)
    inline derived_t& operator++() {
      return ++(this->derived_ref());
    }

    /// Is the iterator at the end of the sequence?
    // (SUBCLASS RESPONSIBILITY)
    inline bool off() const {
      return (this->derived_ref()).off();
    }

  public: // Features inherited without change

    /// Element access (SUBCLASS RESPONSIBILITY)
    inline value_t operator() () const {
      return (this->derived_ref()).value();
    }

    /// Pre-increment - this is much slower than post-increment, don't use it if not necessary
    inline derived_t operator++(int) { 
      derived_t tmp((this->derived_ref())); ++((this->derived_ref())); return tmp; 
    }

  };
  
  //----------------------------------------------------------
  /**
   *  Default constant iterator for indexed objects.
   *  The implementation explicitely maintains an index,
   *  which is incremented and used for accessing the
   *  indexed object through value()
   */
  template <class G_indexed>
  class indexed_default_const_iterator: 
    public indexed_const_iterator_base<typename G_indexed::value_t,
				       G_indexed::rank_c,
				       indexed_default_const_iterator<G_indexed>,
                                       const G_indexed * restrict> {
  public:
    typedef indexed_const_iterator_base<typename G_indexed::value_t,
				        G_indexed::rank_c,
				        indexed_default_const_iterator<G_indexed>,
                                        const G_indexed * restrict>                 super_t;
    typedef indexed_default_const_iterator<G_indexed>                               this_t;
    typedef typename super_t::derived_t                                             derived_t;
    typedef typename super_t::user_defined_data_t                                   user_defined_data_t;
    typedef typename super_t::value_t                                               value_t;
    typedef typename super_t::subscript_t                                           subscript_t;

  protected:

    subscript_t index_;
    subscript_t extent_;

  public:

    /// Explicit init
    inline indexed_default_const_iterator(const user_defined_data_t& wrapped,
					  const subscript_t& i0):
      super_t(wrapped), 
      index_(i0), 
      extent_(wrapped->extent()) {
    }

  public:

    /// Element access 
    inline value_t value() const {
      SL_REQUIRE("Not off", !off());
      return (*(this->data_))(index_);
    }

    /// Current index 
    inline subscript_t index() const {
      return index_;
    }

    /// Pre-increment 
    inline derived_t& operator++() {
      SL_REQUIRE("Not off", !off());
      index_.increment(extent_);
      return (this->derived_ref());
    }

    /// Is the iterator at the end of the sequence?
    inline bool off() const {
      return index_[0] >= extent_[0];
    }

  };
  
  //----------------------------------------------------------
  /**
   *  Default constant iterator wrapper for indexed objects.
   *  The implementation forwards all calls to the 
   *  wrapped object's const iterator.
   */
  template <class G_indexed>
  class indexed_const_iterator_wrapper: 
    public indexed_const_iterator_base<typename G_indexed::value_t,
				       G_indexed::rank_c,
				       indexed_const_iterator_wrapper<G_indexed>,
                                       typename G_indexed::const_iterator> {
  public:
    typedef indexed_const_iterator_base<typename G_indexed::value_t,
				        G_indexed::rank_c,
				        indexed_const_iterator_wrapper<G_indexed>,
                                        typename G_indexed::const_iterator>         super_t;
    typedef indexed_const_iterator_wrapper<G_indexed>                               this_t;
    typedef typename super_t::derived_t                                             derived_t;
    typedef typename super_t::user_defined_data_t                                   user_defined_data_t;
    typedef typename super_t::value_t                                               value_t;
    typedef typename super_t::subscript_t                                           subscript_t;

  public:

    inline indexed_const_iterator_wrapper() {
    }

    inline indexed_const_iterator_wrapper(const user_defined_data_t& wrapped): super_t(wrapped) {
    }

  public:

    /// Element access 
    inline value_t value() const {
      return (this->data_).value();
    }

    /// Current index 
    inline subscript_t index() const {
      return (this->data_).index();
    }

    /// Pre-increment 
    inline derived_t& operator++() {
      ++(this->data_);
      return (this->derived_ref());
    }

    /// Is the iterator at the end of the sequence?
    inline bool off() const {
      return (this->data_).off();
    }
  };

  //----------------------------------------------------------
  /**
   *  Default sparse iterator wrapper for indexed objects.
   *  The implementation forwards all calls to the 
   *  wrapped object's const iterator.
   */
  template <class G_indexed>
  class indexed_const_sparse_iterator_wrapper: 
    public indexed_const_iterator_base<typename G_indexed::value_t,
				       G_indexed::rank_c,
				       indexed_const_sparse_iterator_wrapper<G_indexed>,
                                       typename G_indexed::const_sparse_iterator> {
  public:
    typedef indexed_const_iterator_base<typename G_indexed::value_t,
				        G_indexed::rank_c,
				        indexed_const_sparse_iterator_wrapper<G_indexed>,
                                        typename G_indexed::const_sparse_iterator>super_t;
    typedef indexed_const_sparse_iterator_wrapper<G_indexed>                      this_t;
    typedef typename super_t::derived_t                                             derived_t;
    typedef typename super_t::user_defined_data_t                                   user_defined_data_t;
    typedef typename super_t::value_t                                               value_t;
    typedef typename super_t::subscript_t                                           subscript_t;

  public:
    inline indexed_const_sparse_iterator_wrapper() {
    }

    inline indexed_const_sparse_iterator_wrapper(const user_defined_data_t& wrapped): super_t(wrapped) {
    }

  public:

    /// Element access 
    inline value_t value() const {
      return (this->data_).value();
    }

    /// Current index 
    inline subscript_t index() const {
      return (this->data_).index();
    }

    /// Pre-increment 
    inline derived_t& operator++() {
      ++(this->data_);
      return (this->derived_ref());
    }
 
    /// Is the iterator at the end of the sequence?
    inline bool off() const {
      return (this->data_).off();
    }
  };

};

#endif
