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
/////// ALWAYS TEST IN DEBUG MODE
#if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////
#include <sl/floating_cone.hpp>//mbr: order matters in clang, which is not able to find operator<<
#include <sl/fixed_size_point.hpp>
#include <sl/serializer.hpp>
#include <sl/buffer_serializer.hpp>
#include <sl/std_serializer.hpp>
#include <sl/tester.hpp>
#include <iterator>

// ======= Debug IO
namespace std {

template <class T1, class T2>
inline std::ostream& operator<<(std::ostream& s, const std::pair<T1,T2>& x) {
  return s << "[" << x.first << x.second << "]";
}

template <class T, class Allocator>
inline std::ostream& operator<<(std::ostream& s, const std::list<T,Allocator>& x)  {
  s << "{ ";
  std::copy(x.begin(), x.end(), std::ostream_iterator<T>(s, " "));
  s << "} ";
  return s;
}	  

template <class T, class Allocator>
inline std::ostream& operator << (std::ostream& s, const std:: vector<T,Allocator>& x) {
  s << "{ ";
  std::copy(x.begin(), x.end(), std::ostream_iterator<T>(s, " "));
  s << "} ";
  return s;
}

template <class Key, class Compare, class Allocator>
inline std::ostream& operator << (std::ostream& s,
                                  const std::set<Key,Compare,Allocator>& x) {
  s << "{ ";
  std::copy(x.begin(), x.end(), std::ostream_iterator<Key>(s, " "));
  s << "} ";
  return s;
}	  

template <class Key, class Compare, class Allocator>
inline std::ostream& operator << (std::ostream& s,
                                  const std::multiset<Key,Compare,Allocator>& x) {
  s << "{ ";
  std::copy(x.begin(), x.end(), std::ostream_iterator<Key>(s, " "));
  s << "} ";
  return s;
}	  

template <class Key, class T, class Compare, class Allocator>
inline std::ostream& operator << (std::ostream& s,
                                  const std::map<Key,T,Compare,Allocator>& x) {
  s << "{ ";
  std::copy(x.begin(), x.end(), std::ostream_iterator< std::pair<Key,T> >(s, " "));
  s << "} ";
  return s;
}	  

template <class Key, class T, class Compare, class Allocator>
inline std::ostream& operator << (std::ostream& s,
                                  const std::multimap<Key,T,Compare,Allocator>& x) {
  s << "{ ";
  std::copy(x.begin(), x.end(), std::ostream_iterator< std::pair<Key,T> >(s, " "));
  s << "} ";
  return s;
}	  

}
//=======================

static std::size_t failed_test_count = 0;

//--------------------

template <class T>
void test_serialization(sl::tester& tester,
                        const std::string& msg,
                        const T& x) {
  sl::input_buffer_serializer ibuf;
  sl::output_buffer_serializer obuf;

  signed char z = 'z';
  
  obuf << x << z << x;

  ibuf.buffer() = obuf.buffer(); ibuf.reset();
  tester.test(msg+" -- off at begin", !ibuf.off());

  T x_prime, x_second; signed char z_prime = 'y';
  ibuf >> x_prime >> z_prime >> x_second;

  tester.test(msg+" -- at begin", x_prime, x);
  tester.test(msg+" -- separator", z_prime, z);
  tester.test(msg+" -- at end", x_second, x);
  tester.test(msg+" -- off at end", ibuf.off());
}

void test_serializer() {
  sl::tester tester(std::string() + "serializer");

  test_serialization(tester, "float serialization", 1.0f);
  test_serialization(tester, "double serialization", 1.0);
  test_serialization(tester, "int serialization", 1);
  test_serialization(tester, "bool serialization", true);

  test_serialization(tester, "std::string serialization", std::string("A quick brown..."));

  std::list<double> a_lst; a_lst.push_back(1.0); a_lst.push_back(2.0);
  test_serialization(tester, "std::list<double> serialization", a_lst);

  std::vector<double> a_vec; a_vec.push_back(1.0); a_vec.push_back(2.0);
  test_serialization(tester, "std::vector<double> serialization", a_vec);

  std::set<double> a_set; a_set.insert(1.0); a_set.insert(2.0);
  test_serialization(tester, "std::set<double> serialization", a_set);

  std::map<double,double> a_map; a_map[1.0] = 2.0; a_map[2.0] = 3.0;
  test_serialization(tester, "std::map<double,double> serialization", a_map);

  std::vector< std::set<double> > a_vec_set; a_vec_set.push_back(a_set); a_vec_set.push_back(a_set);
  test_serialization(tester, "std::vector< std::set<double> > serialization", a_vec_set);

  std::vector<sl::uint64_t> vec64; vec64.push_back(1 << 31); vec64.push_back(1LL << 48);
  test_serialization(tester, "std::vector<int64_t> serialization", vec64);

  test_serialization(tester, "point3f serialization", sl::point3f(1.0f,2.0f,3.0f));
  test_serialization(tester, "floating_cone3f serialization", sl::floating_cone3f(sl::vector3f(1.0f,0.0f,0.0f), 0.1f));
  test_serialization(tester, "ball3f serialization", sl::ball3f(sl::point3f(1.0f,0.0f,0.0f), 0.1f));

  failed_test_count += tester.failed_test_count();
}

int main() {
  test_serializer();
  return (int)failed_test_count;
}
