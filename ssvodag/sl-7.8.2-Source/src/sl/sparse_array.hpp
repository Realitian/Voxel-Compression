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
#ifndef SL_SPARSE_ARRAY_HPP
#define SL_SPARSE_ARRAY_HPP

#include <sl/indexed_container.hpp>
#include <map> // std::map

namespace sl {

  /**
   *  Indexed iterator view over a std::map<index, value>
   */
  template <
    class  G_numtype, 
    size_t G_rank
  >
  class indexed_const_sparse_array_iterator:
    public indexed_const_iterator_base<G_numtype,
                                       G_rank,
                                       indexed_const_sparse_array_iterator<G_numtype, G_rank>,
                                       typename std::map<sl::index<G_rank>,G_numtype>::const_iterator> {
  public:
    typedef  indexed_const_iterator_base<G_numtype,
                                       G_rank,
                                       indexed_const_sparse_array_iterator<G_numtype, G_rank>,
                                       typename std::map<sl::index<G_rank>,G_numtype>::const_iterator> super_t;
    typedef indexed_const_sparse_array_iterator<G_numtype, G_rank>                  this_t;
    typedef typename super_t::derived_t                                             derived_t;
    typedef typename super_t::user_defined_data_t                                   user_defined_data_t;
    typedef typename super_t::value_t                                               value_t;
    typedef typename super_t::subscript_t                                           subscript_t;

    typedef typename std::map<subscript_t,value_t>::const_iterator                  wrapped_iter_t;

  protected:

    wrapped_iter_t it_end_;
    subscript_t index_end_mark_;

  public:

    indexed_const_sparse_array_iterator(const wrapped_iter_t& it_start,
					       const wrapped_iter_t& it_end,
					       const subscript_t& index_end_mark)
      : 
      super_t(it_start), 
      it_end_(it_end),
      index_end_mark_(index_end_mark)
    {
    }

  public:

    /// Element access 
    inline value_t operator() () const {
      SL_REQUIRE("Not off", (this->data_) != it_end_);
      return (this->data_)->second;
    }

    /// Current index 
    inline subscript_t index() const {
      return ((this->data_) == it_end_) ? index_end_mark_ : (this->data_)->first;
    }

    /// Pre-increment 
    inline derived_t& operator++() {
      SL_REQUIRE("Not off", (this->data_) != it_end_);
      ++(this->data_);
      return this->derived_ref();
    }
 
    /// Is other pointing to the same element?
    inline bool operator==(const this_t& other) const {
      return (this->data_) == (other.data_);
    }
  };



  /**
   *  N-dimensional sparse arrays
   */
  template <
    class  G_numtype, 
    size_t G_rank, 
    class  G_discriminant
  >
  class sparse_array: 
    public indexed_container<G_numtype, 
                             G_rank, 
                             sparse_array<G_numtype,G_rank,G_discriminant>, 
                             std::map< ::sl::index<G_rank>, G_numtype > , 
                             G_discriminant > {
  public:
    
    typedef indexed_container<G_numtype, 
                             G_rank, 
                             sparse_array<G_numtype,G_rank,G_discriminant>, 
                             std::map< ::sl::index<G_rank>, G_numtype > , 
                             G_discriminant >                                super_t;
    typedef sparse_array<G_numtype, G_rank,  G_discriminant >                this_t;
    typedef typename super_t::derived_t                                      derived_t;
    typedef typename super_t::user_defined_data_t                            user_defined_data_t;
    typedef typename super_t::value_t                                        value_t;
    typedef typename super_t::discriminant_t                                 discriminant_t;
    typedef typename super_t::subscript_t                                    subscript_t;
    
    typedef user_defined_data_t                                              map_t;

    typedef indexed_default_const_iterator<this_t>                           const_iterator;
    typedef indexed_const_sparse_array_iterator<G_numtype,G_rank>            const_sparse_iterator;
    
    typedef sparse_tag                                                       sparsity_t;

    /// Never return an object of this type on the stack!
    enum { xpr_pass_by_value = false };

  protected:

    subscript_t extent_;

  private: // disable copy, should pass through expressions

    sparse_array( const sparse_array& );
    const sparse_array& operator=( const sparse_array& );
	
  public:

    /// Default init (no operation)
    explicit inline sparse_array( subscript_t extent = subscript_t())
      : 
      super_t(map_t()), 
      extent_(extent)
    {
      SL_REQUIRE("Good extent", extent == subscript_t() || extent.element_count() > 0);
    }

    inline ~sparse_array() {
      resize(subscript_t());
    }

    explicit inline sparse_array(size_t i0)
      : 
      super_t(map_t()),
      extent_(subscript_t(i0))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0).element_count() > 0);
    }

    explicit inline sparse_array(size_t i0, size_t i1)
      : 
      super_t(map_t()),
      extent_(subscript_t(i0,i1))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1).element_count() > 0);
    }

    explicit inline sparse_array(size_t i0, size_t i1, size_t i2)
      : 
      super_t(map_t()),
      extent_(subscript_t(i0,i1,i2))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1,i2).element_count() > 0);
    }

    explicit inline sparse_array(size_t i0, size_t i1, size_t i2, size_t i3)
      : 
      super_t(map_t()),
      extent_(subscript_t(i0,i1,i2,i3))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1,i2,i3).element_count() > 0);
    }

    explicit inline sparse_array(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4)
      : 
      super_t(map_t()),
      extent_(subscript_t(i0,i1,i2,i3,i4))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1,i2,i3,i4).element_count() > 0);
    }
    
  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << extent_;
      super_t::store_to(s);
    }
    
    void retrieve_from(input_serializer& s) {
      s >> extent_;
      super_t::retrieve_from(s);
    }

  public: // Manifest init
      
   /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    indexed_manifest_initializer<this_t> operator= (value_t x) {      
      return indexed_manifest_initializer<this_t>(this->derived_ref(), x);
    }

  public: // Features from indexed

    /// The number of indexed elements for each of the ranks
    inline const subscript_t& extent() const {
      return extent_;
    }

    /// The element referenced by subscript idx
    inline value_t item(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      typename map_t::const_iterator it = (this->data_).find(idx);

      return 
	(it == (this->data_).end()) ?
	scalar_math<value_t>::zero() :
	it->second;
    }
    
    /// Iterator pointing to the beginning of the element sequence
    inline const_iterator begin() const {
      return const_iterator(this,subscript_t());
    }

    /// Iterator pointing to the end of the element sequence
    inline const_iterator end() const {
      return const_iterator(this,extent_.end_mark());
    }

    /// Iterator pointing to the beginning of the non zero element sequence
    inline const_sparse_iterator sparse_begin() const {
      return const_sparse_iterator((this->data_).begin(), (this->data_).end(), extent_.end_mark());
    }

    /// Iterator pointing to the end of the non zero element sequence
    inline const_sparse_iterator sparse_end() const {
      return const_sparse_iterator((this->data_).end(), (this->data_).end(), extent_.end_mark());
    }

  public: // Assignment

    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline this_t& operator= (const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
      resize(rhs.extent());
      return assign_from(rhs);
    }

    // TODO: Assign from manifest array

  public: // Update
    
    /// Set element at index idx to zero
    inline void put_zero(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->data_).erase(idx);
    }

    /// Set element at index idx to the value v, different than zero
    inline void put_non_zero(const value_t& v, const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->data_)[idx] = v;
    }

    /// Set element at index idx to the value v
    inline void put(const value_t& v, const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      if (is_zero(v)) {
	put_zero(idx);
      } else {
	put_non_zero(v,idx);
      }
    }

    /// Write access to element at subscript idx
    inline value_t& item_ref(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->data_)[idx];
    }

    /// Is current resizable to size sz? 
    inline bool resizable(const subscript_t&) const {
      return true; // Always resizable
    }	

    /// If resizable, make this size equal to sz 
    void resize(const subscript_t& sz) {
      SL_REQUIRE("Resizable", resizable(sz));
      if (extent_ != sz) {
	size_t n = sz.element_count();
	if (n == 0) {
	  /// Clear
	  (this->data_).clear();
	  extent_ = subscript_t();
	} else {
	  (this->data_).clear();
	  extent_ = sz;
	}
      }
      SL_ENSURE("Good size", sz == extent());
    }
    
    /// Set everything to zero 
    void clear() {
      (this->data_).clear();
    }
    
  };

}


  
#endif
