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

// Template expression tests

#if 0
//=============================================================

#include <sl/xpr/dense_array.hpp>
#include <sl/xpr/sparse_array.hpp>

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
    dump("sl::sqr(std::sqrt(u))", sl::sqr(std::sqrt(u)));
    dump("u + v", u+v);
    u = u+v;
    dump("u = u + v", u);
    dump("2 * u", 2.0f * u);
    dump("u * 2", 2.0f * u);
    std::cerr << "-------------------" << std::endl;
    dump("u(i,j)", u(sl::tensor::i,sl::tensor::j));
    dump("u(j,i)", u(sl::tensor::j,sl::tensor::i));
    
    std::cerr << "===================" << std::endl;
  }

};

void tensor_test() {
  sl::dense_array<float, 1, void> x(4), y(4);
  sl::dense_array<float, 2, void> A(4,4);

  x = 1,2,3,4;
  y = 1,0,0,1;

  //  A = x(sl::tensor::i) * y(sl::tensor::j);

  dump("x_i * y_j", A);

};
 
int main() {
  basic_array_test<false,false>::execute();
  basic_array_test<false,true>::execute();
  basic_array_test<true,false>::execute();
  basic_array_test<true,true>::execute();
  tensor_test();

  return 0;
}

//=============================================================
#else


#include <sl/dense_vector.hpp>

///////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------//
// Function to print out the contents of a dense_vector.
//---------------------------------------------------------------------------//


template<class P, class E, class O, class D> 
static void dump( const sl::xpr::vec_expression<P,E,O,D>& u ) {
  for( int i=0; i < u.size(); i++ ) {
    std::cout << "u[" << i << "]=" << u(i) << std::endl;
  }
}

template<class P,class O,class D> 
static void dump( const sl::dense_vector<P,O,D>& u ) {
  dump(u.as_expression());
}


//---------------------------------------------------------------------------//
// Check that we can index dense_vector the usual way.
//---------------------------------------------------------------------------//

void t1() {
    std::cout << "t1: start\n";

    sl::dense_vector<float> u( 5);

    for( int i=0; i < 5; i++ )
	u(i) = 2.*i;

    dump(u);

    std::cout << "t1: end\n";
}

//---------------------------------------------------------------------------//
// Check a simple binary operation, with an assignment.
//---------------------------------------------------------------------------//

void t2() {
    std::cout << "t2: start\n";

    sl::dense_vector<float> u( 5), v( 5), w( 5);

    for( int i=0; i < 5; i++ ) {
	v(i) = 2.*i;
	w(i) = 5.-i;
    }
    u = 0;
    u += v + w;
    u = v + w;

    dump(u);

    std::cout << "t2: end\n";
}


//---------------------------------------------------------------------------//
// Mixed mode
//---------------------------------------------------------------------------//

void t2_mix()
{
    std::cout << "t2: start\n";

    sl::dense_vector<float>  u( 5);
    sl::dense_vector<double> v( 5), w( 5);

    for( int i=0; i < 5; i++ ) {
	v(i) = 2.*i;
	w(i) = 5.-i;
    }
    u = 0;
    u += v + w;
    u = v + w;

    dump(u);

    std::cout << "t2: end\n";
}
  
//---------------------------------------------------------------------------//
// Now a more involved test of binary operations.
//---------------------------------------------------------------------------//

void t3()
{
    std::cout << "t3: start\n";

    sl::dense_vector<float> a(5), b(5), c(5), d(5), e(5), f(5);

    a = 4.;
    b = 2.;
    c = 2.;
    d = 5.;
    e = 2.;

    //f = (a+b)*(c+d);
    //f = (a+b)*c/(d-e);

    dump(f);

    std::cout << "t3: end\n";
}

//---------------------------------------------------------------------------//
// Now test distributing function application
//---------------------------------------------------------------------------//

void t4() {
    std::cout << "t4: start\n";

    sl::dense_vector<float> a(5), b(5), c(5), d(5);

    a = 4.;
    b = 2.;

    c = std::pow( a, 2.f );

    dump(c);

    c = 2.;
    // d = std::pow( c, 2 );

    dump(d);

    // d = 1.f + std::pow( c, 2 );

    dump(d);

    std::cout << "t4: end\n";
}

//---------------------------------------------------------------------------//
// Check incompatible participation.
//---------------------------------------------------------------------------//

class foobar_domain {};

void t5()
{
    std::cout << "t5: start\n";

    sl::dense_vector<float> a(5), b(5);
    sl::dense_vector<float,sl::xpr::vec_column_orientation_tag,foobar_domain> c(5);

    a = -99.;

    b = 2.; c = 3.;

// This statement should be illegal!
//     a = b + c;

    dump(a);

    std::cout << "t5: end\n";
}

//---------------------------------------------------------------------------//
// Main program, just run through each test in turn.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    t1();
    t2();
    t2_mix();
    t3();
    t4();
    t5();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstUV.cc
//---------------------------------------------------------------------------//



#endif
