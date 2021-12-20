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
/////// ALWAYS TEST IN DEBUG MODE, UNLESS IT BREAKS THE COMPILER
# if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
# endif 
///////

// Indexed 

//=============================================================

#include <sl/dense_array.hpp>
#include <sl/sparse_array.hpp>

template <class G_numtype, size_t G_rank, class G_derived, class G_userdefined, class G_discriminant >
static void dump(std::string s,
		 const sl::indexed<G_numtype,G_rank,G_derived,G_userdefined,G_discriminant>& u) {
  std::cerr 
    << s 
    << " ( " 
    << (sl::is_same< typename G_derived::sparsity_t, sl::dense_tag >::value ? "dense" : "sparse")
    << " ) = " << std::endl
    << u <<
    std::endl;
}

template <bool t1_sparse, bool t2_sparse> 
class basic_array_test {
public:

  typedef typename sl::gen_if<
    t1_sparse,
    sl::sparse_array<float, 2, void>,
    sl::dense_array<float, 2, void>
  >::type array1_t;

  typedef typename sl::gen_if<
    t2_sparse,
    sl::sparse_array<float, 2, void>,
    sl::dense_array<float, 2, void>
  >::type array2_t;

    
  static void execute() {
    array1_t u(2, 3);
    array2_t v(2, 3);

    std::cerr << "===================" << std::endl;
    dump("u", u);
    u(1,2) = 12;
    dump("u(1,2) = 12; u", u);
    dump("v", v);
    v(1,1) = 11;
    dump("v(1,1) = 11; v", v);
#if 0
    dump("sl::sqr(std::sqrt(u))", sl::sqr(std::sqrt(u)));
    dump("u + v", u+v);
    u = u+v;
    dump("u = u + v", u);
    dump("2 * u", 2.0f * u);
    dump("u * 2", 2.0f * u);
    std::cerr << "-------------------" << std::endl;
    dump("u(i,j)", u(sl::tensor::i,sl::tensor::j));
    dump("u(j,i)", u(sl::tensor::j,sl::tensor::i));
#endif    
    std::cerr << "===================" << std::endl;
  }

};
 
int main() {
  basic_array_test<false,false>::execute();
  basic_array_test<false,true>::execute();
  basic_array_test<true,false>::execute();
  basic_array_test<true,true>::execute();

  return 0;
}
