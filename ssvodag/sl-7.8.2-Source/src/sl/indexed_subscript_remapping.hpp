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
#ifndef SL_INDEXED_SUBSCRIPT_REMAPPING_HPP
#define SL_INDEXED_SUBSCRIPT_REMAPPING_HPP

#include <sl/indexed.hpp>

namespace sl {

  namespace detail {
    
    template <size_t G_rank>
    class subscript_remapper {
    public:
      template <class  G_numtype, size_t G_srcrank, class  G_derived, class  G_userdefined, class  G_discriminant>
      inline static G_numtype item(const indexed<G_numtype, G_srcrank, G_derived, G_userdefined, G_discriminant>& a,
				   const ::sl::index<G_rank>& idx,
				   size_t i0, size_t i1, size_t i2, size_t i3,
				   size_t i4, size_t i5, size_t i6, size_t i7,
				   size_t i8, size_t i9, size_t i10) {
	::sl::index<G_srcrank> a_idx;
	switch (G_srcrank) {
	case 0: 
	  break;
	case 1: 
	  a_idx[0] = idx[i0]; 
	  break;
	case 2: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  break;
	case 3: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  break;
	case 4: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  break;
	case 5: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  break;
	case 6: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  a_idx[5] = idx[i5]; 
	  break;
	case 7: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  a_idx[5] = idx[i5]; 
	  a_idx[6] = idx[i6]; 
	  break;
	case 8: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  a_idx[5] = idx[i5]; 
	  a_idx[6] = idx[i6]; 
	  a_idx[7] = idx[i7]; 
	  break;
	case 9: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  a_idx[5] = idx[i5]; 
	  a_idx[6] = idx[i6]; 
	  a_idx[7] = idx[i7]; 
	  a_idx[8] = idx[i8]; 
	  break;
	case 10: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  a_idx[5] = idx[i5]; 
	  a_idx[6] = idx[i6]; 
	  a_idx[7] = idx[i7]; 
	  a_idx[8] = idx[i8]; 
	  a_idx[9] = idx[i9]; 
	  break;
	case 11: 
	  a_idx[0] = idx[i0]; 
	  a_idx[1] = idx[i1]; 
	  a_idx[2] = idx[i2]; 
	  a_idx[3] = idx[i3]; 
	  a_idx[4] = idx[i4]; 
	  a_idx[5] = idx[i5]; 
	  a_idx[6] = idx[i6]; 
	  a_idx[7] = idx[i7]; 
	  a_idx[8] = idx[i8]; 
	  a_idx[9] = idx[i9]; 
	  a_idx[10] = idx[i10]; 
	  break;
	default:
	  SL_FAIL("Sorry, too many subscripts, not implemented!");
	  break;
	}
	return a(a_idx);
      }
    };
  }
		    
  //----------------------------------------------------------
  template <
    size_t N0,
    size_t N1 = 0,
    size_t N2 = 0,
    size_t N3 = 0,
    size_t N4 = 0,
    size_t N5 = 0,
    size_t N6 = 0,
    size_t N7 = 0,
    size_t N8 = 0,
    size_t N9 = 0,
    size_t N10 = 0
  >
  class gen_max_size_t {
  public:
   enum {
        tmp1 =  (N0 > N1  )  ? N0  : N1,
        tmp2 =  (N2 > tmp1)  ? N2  : tmp1,
        tmp3 =  (N3 > tmp2)  ? N3  : tmp2,
        tmp4 =  (N4 > tmp3)  ? N4  : tmp3,
        tmp5 =  (N5 > tmp4)  ? N5  : tmp4,
        tmp6 =  (N6 > tmp5)  ? N6  : tmp5,
        tmp7 =  (N7 > tmp6)  ? N7  : tmp6,
        tmp8 =  (N8 > tmp7)  ? N8  : tmp7,
        tmp9 =  (N9 > tmp8)  ? N9  : tmp8,
        tmp10 = (N10 > tmp9) ? N10 : tmp9
    };

    enum {
      value = tmp10+1
    };
  };

  //----------------------------------------------------------
  /**
   *  Index remapping expression
   */
  template <
    class  G_indexed, 
    size_t N0,
    size_t N1,
    size_t N2,
    size_t N3,
    size_t N4,
    size_t N5,
    size_t N6,
    size_t N7,
    size_t N8,
    size_t N9,
    size_t N10
  >
  class indexed_subscript_remapping: 
    public indexed<typename G_indexed::value_t,
                   gen_max_size_t<N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10>::value,
                   indexed_subscript_remapping<G_indexed,N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10>,
                   G_indexed,
                   typename G_indexed::discriminant_t> {
  public:
    typedef indexed_subscript_remapping<G_indexed,N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10> this_t;
    typedef indexed<typename G_indexed::value_t,
                    gen_max_size_t<N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10>::value,
                    this_t,
                    G_indexed,
                    typename G_indexed::discriminant_t>                    super_t;

    typedef typename super_t::derived_t                                    derived_t;
    typedef typename super_t::user_defined_data_t                          user_defined_data_t;
    typedef typename super_t::value_t                                      value_t;
    typedef typename super_t::discriminant_t                               discriminant_t;
    typedef typename super_t::subscript_t                                  subscript_t;

    typedef dense_tag                                                      sparsity_t;
    
    typedef indexed_default_const_iterator<this_t>                         const_iterator;
    typedef indexed_default_const_iterator<this_t>                         const_sparse_iterator;
    
    typedef G_indexed                                                      arg_t;

    /// True iff the dimensions are adaptable inside an expression
    enum { xpr_dynamic_dimensions = true };

    /// The rank of the object
    enum { rank_c = super_t::rank_c };

  protected:

    subscript_t extent_; 

  public:
    
    /// Data init
    explicit inline indexed_subscript_remapping(const arg_t& arg): super_t(user_defined_data_t(arg)) {
      size_t index_map[11] = {N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10 };
      const typename arg_t::subscript_t arg_extent = (this->data_).extent();
      for (size_t r=0; r < rank_c; r++) {
	extent_[r] = index_map[r] < arg_extent.count() ? arg_extent(index_map[r]) : 1;
	std::cerr << "Map: " << index_map[r] << " --> " << r << std::endl;
      }
    }

    /// The number of indexed elements for each of the ranks
    inline subscript_t extent() const {
      return extent_;
    }

    /// The element referenced by subscript idx
    inline value_t item(const subscript_t& idx) const {
      SL_REQUIRE("Good subscript", this->good_subscript(idx));

      return detail::subscript_remapper<this_t::rank_c>::item((this->data_), idx,
                                                              N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10);
    }
    
    /// Iterator pointing to the beginning of the element sequence
    inline const_iterator begin() const {
      return const_iterator(this,subscript_t());
    }

    /// Iterator pointing to the beginning of the non zero element sequence
    inline const_sparse_iterator sparse_begin() const {
      return begin();
    }

  public: // Dynamic expression bounds

    /**
     *  The upper bound for rank r when this is part of an 
     *  expresssion. This is extent(r) if r has a fixed dimension,
     *  XPR_DYNAMIC_BOUND otherwise.
     */
    inline size_t xpr_ubound(size_t r) const {
      if (N0 == r) {
	return (this->data_).xpr_ubound(0);
      } else if ((N1 == r) && (rank_c > 1)) {
	return (this->data_).xpr_ubound(1);
      } else if ((N2 == r) && (rank_c > 2)) {
	return (this->data_).xpr_ubound(2);
      } else if ((N3 == r) && (rank_c > 3)) {
	return (this->data_).xpr_ubound(3);
      } else if ((N4 == r) && (rank_c > 4)) {
	return (this->data_).xpr_ubound(4);
      } else if ((N5 == r) && (rank_c > 5)) {
	return (this->data_).xpr_ubound(5);
      } else if ((N6 == r) && (rank_c > 6)) {
	return (this->data_).xpr_ubound(6);
      } else if ((N7 == r) && (rank_c > 7)) {
	return (this->data_).xpr_ubound(7);
      } else if ((N8 == r) && (rank_c > 8)) {
	return (this->data_).xpr_ubound(8);
      } else if ((N9 == r) && (rank_c > 9)) {
	return (this->data_).xpr_ubound(9);
      } else if ((N10 == r) && (rank_c > 10)) {
	return (this->data_).xpr_ubound(10);
      } else {
	return XPR_DYNAMIC_BOUND;
      }
    }

  };

}


#endif
