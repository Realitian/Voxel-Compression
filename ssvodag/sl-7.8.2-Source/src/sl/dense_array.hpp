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
#ifndef SL_DENSE_ARRAY_HPP
#define SL_DENSE_ARRAY_HPP

#include <sl/smart_pointer.hpp>
#include <sl/indexed_container.hpp>

namespace sl {


  /**
   *  N-dimensional dense arrays
   */
  template <
    class  G_numtype, 
    size_t G_rank, 
    class  G_discriminant
  >
  class dense_array: 
    public indexed_container<G_numtype, 
                             G_rank, 
                             dense_array<G_numtype,G_rank,G_discriminant>, 
                             G_numtype*, 
                             G_discriminant > {
  public:
    
    typedef indexed_container<G_numtype, 
                              G_rank, 
                              dense_array<G_numtype,G_rank,G_discriminant>, 
                              G_numtype*, 
                              G_discriminant >                               super_t;
    typedef dense_array<G_numtype, G_rank,  G_discriminant >                 this_t;
    typedef typename super_t::derived_t                                      derived_t;
    typedef typename super_t::user_defined_data_t                            user_defined_data_t;
    typedef typename super_t::value_t                                        value_t;
    typedef typename super_t::discriminant_t                                 discriminant_t;
    typedef typename super_t::subscript_t                                    subscript_t;
    
    typedef indexed_default_const_iterator<this_t>                           const_iterator;
    typedef indexed_default_const_iterator<this_t>                           const_sparse_iterator; // HACK
    
    typedef dense_tag                                                        sparsity_t;

    typedef G_numtype*                                                       data_pointer_t;
  protected:

    subscript_t extent_;
    subscript_t stride_;

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << extent_ << stride_;
      super_t::store_to(s);
    }
    
    void retrieve_from(input_serializer& s) {
      s >> extent_ >> stride_;
      super_t::retrieve_from(s);
    }

  protected: // Helper

    inline void compute_strides() {
      if (G_rank == 1) {
	stride_(0) = 1;
      } else {
	size_t s = 1;
	for (size_t n=0; n<G_rank; ++n) {
	  stride_(n) = s;
	  s *= extent_(n);
	}
      }
    }
	
  public:

    /// Default init (no operation)
    explicit inline dense_array(subscript_t extent = subscript_t())
      : 
      super_t(data_pointer_t(0)),
      extent_(extent)
    {
      SL_REQUIRE("Good extent", extent == subscript_t() || extent.element_count() > 0);
      compute_strides();
      if (this->count()) this->data_ = (new value_t[this->count()]);
      if (!detail::init_value<value_t>::init_value_is_default) clear(); 
    }

    inline ~dense_array() {
      // Smart pointer handles delete
      if (this->data_) delete[] this->data_;
      this->data_ = 0;
    }

    explicit inline dense_array(size_t i0)
      : 
      super_t(data_pointer_t(0)),
      extent_(subscript_t(i0))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0).element_count() > 0);
      compute_strides();
      if (this->count()) this->data_ = (new value_t[this->count()]);
      if (!detail::init_value<value_t>::init_value_is_default) clear(); 
    }

    explicit inline dense_array(size_t i0, size_t i1)
      : 
      super_t(data_pointer_t(0)),
      extent_(subscript_t(i0,i1))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1).element_count() > 0);
      compute_strides();
      if (this->count()) this->data_ = (new value_t[this->count()]);
      if (!detail::init_value<value_t>::init_value_is_default) clear(); 
    }

    explicit inline dense_array(size_t i0, size_t i1, size_t i2)
      : 
      super_t(data_pointer_t(0)),
      extent_(subscript_t(i0,i1,i2))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1,i2).element_count() > 0);
      compute_strides();
      if (this->count()) this->data_ = (new value_t[this->count()]);
      if (!detail::init_value<value_t>::init_value_is_default) clear(); 
    }

    explicit inline dense_array(size_t i0, size_t i1, size_t i2, size_t i3)
      : 
      super_t(data_pointer_t(0)),
      extent_(subscript_t(i0,i1,i2,i3))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1,i2,i3).element_count() > 0);
      compute_strides();
      if (this->count()) this->data_ = (new value_t[this->count()]);
      if (!detail::init_value<value_t>::init_value_is_default) clear(); 
     }

    explicit inline dense_array(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4)
      : 
      super_t(data_pointer_t(0)),
      extent_(subscript_t(i0,i1,i2,i3,i4))
    {
      SL_REQUIRE("Good extent", extent_ == subscript_t() || subscript_t(i0,i1,i2,i3,i4).element_count() > 0);
      compute_strides();
      if (this->count()) this->data_ = (new value_t[this->count()]);
      if (!detail::init_value<value_t>::init_value_is_default) clear(); 
    }

    
  public: // Copy 

#if 0
    /**
     *  Constructor making a shared view of other's data. After
     *  this constructor is called, both arrays point to the
     *  same data. 
     */
    inline dense_array(const this_t& other): 
      super_t(other.data_),
      extent_(other.extent_),
      stride_(other.stride_) {
    }
#else
    /**
     *  Constructor assigning from other. After
     *  this constructor is called, the arrays point to two
     *  different copies of the same data. 
     */
    inline dense_array(const this_t& other): 
      super_t(data_pointer_t(0)),
      extent_(subscript_t()),
      stride_(subscript_t()) {
      resize(other.extent()); // FIXME optimize
      this->assign_from(other);
    }
    
#endif
       
    inline this_t& operator=(const this_t& other) {
      resize(other.extent());
      return this->assign_from(other);
    }


    template <class G_numtype2, class G_derived2, class G_userdefined2> 
    inline this_t& operator= (const indexed<G_numtype2, G_rank, G_derived2, G_userdefined2, G_discriminant>& rhs) {
      resize(rhs.extent());
      return this->assign_from(rhs);
    }

    
    /**
     *  A copy of current array, not shared with anyone.
     */
    inline this_t copy() const {
      this_t result(extent());
      return result.assign_from(*this);
    }
     
  public:

    /// Offset into array
    inline size_t offset(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      size_t result = 0;
      for (size_t r=0; r<G_rank; ++r) {
	result += idx(r) * stride_(r);
      }
      return result;
    }
      
  public: // Features from indexed

    /// The number of indexed elements for each of the ranks
    inline const subscript_t& extent() const {
      return extent_;
    }

    /// The element referenced by subscript idx
    inline value_t item(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->data_)[offset(idx)];
    }
    
    /// Iterator pointing to the beginning of the element sequence
    inline const_iterator begin() const {
      return const_iterator(this,subscript_t());
    }

    /// Iterator pointing to the beginning of the non zero element sequence
    inline const_sparse_iterator sparse_begin() const {
      return const_sparse_iterator(this,subscript_t());
    }

  public: // Assignment

   /** 
     *  Initialize from manifest constant.
     *  This allows initializations such as v = 1.0; (fill) and 
     *  v = 1.0, 2.0, 3.0; (component init)
     */
    inline indexed_manifest_initializer<this_t> operator= (value_t x) {      
      return indexed_manifest_initializer<this_t>(this->derived_ref(), x);
    }
    
  public: // Update
    
    /// Set element at index idx to zero
    inline void put_zero(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->data_)[offset(idx)] = scalar_math<value_t>::zero();
    }

    /// Set element at index idx to the value v, different than zero
    inline void put_non_zero(const value_t& v, const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->data_)[offset(idx)] = v;
    }

    /// Set element at index idx to the value v
    inline void put(const value_t& v, const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      (this->data_)[offset(idx)] = v;
    }

    /// Write access to element at subscript idx
    inline value_t& item_ref(const subscript_t& idx) {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));
      return (this->data_)[offset(idx)];
    }

    /// Is current resizable to size sz? 
    inline bool resizable(const subscript_t&) const {
      return true; // Always resizable
    }	

    /// If resizable, make this size equal to sz 
    inline void resize(const subscript_t& sz) {
      SL_REQUIRE("Resizable", resizable(sz));
      if (extent_ != sz) {
	size_t n = sz.element_count();
	if (n == 0) {
	  /// Clear
          if (this->data_) delete[] this->data_;
	  this->data_ = NULL;
	  extent_ = subscript_t();
	  stride_ = subscript_t();
	} else {
          if (this->data_) delete[] this->data_;
	  this->data_ = new value_t[n];
	  extent_ = sz;
	  compute_strides();
	  clear();
	}
      }
      SL_ENSURE("Good size", sz == extent());
    }
    
    /// Set everything to zero 
    inline void clear() {
      size_t n = extent_.element_count();
      const value_t Zero = detail::init_value<value_t>::value();
      for (size_t i=0; i<n; ++i) {
	(this->data_)[i] = Zero;
      }
    }
    
  public: // Raw data access
    
    /// The raw array pointer - very dangerous!
    sized_raw_array_pointer<value_t> as_sized_raw_array_pointer() {
      return sized_raw_array_pointer<value_t>(this->data_, this->count());
    }

    /// The raw array pointer - very dangerous!
    sized_raw_array_pointer<const value_t> as_const_sized_raw_array_pointer() const {
      return sized_raw_array_pointer<const value_t>(this->data_, this->count());
    }
    
  };

}

#endif
