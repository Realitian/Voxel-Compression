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
# if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif

#include <sl/smart_pointer.hpp>
#include <sl/tester.hpp>


//--------------------

static std::size_t failed_test_count = 0;

//--------------------

class test_not_intrusive {
public:
  static std::size_t live_instance_count;

  test_not_intrusive() {
    ++live_instance_count;
  }
  ~test_not_intrusive() {
    --live_instance_count;
  }
};

std::size_t test_not_intrusive::live_instance_count = 0;

class test_intrusive: public sl::reference_counted {
public:
  static std::size_t live_instance_count;

  test_intrusive() {
    ++live_instance_count;
  }
  ~test_intrusive() {
    --live_instance_count;
  }
};

std::size_t test_intrusive::live_instance_count = 0;

//--------------------

template<class T>
static void test_shared_pointer(const char *id, 
				const T&,
				const std::size_t expected_pointer_size) {

  typedef sl::shared_pointer<T> pointer_t;
  
  sl::tester tester(id);
  tester.test("sizeof(pointer)", sizeof(pointer_t), expected_pointer_size);
  
  std::size_t cnt = T::live_instance_count;
  pointer_t ptr1(new T);
  
  tester.test("use count at creation", ptr1.use_count(), 1L);
  tester.test("instance count at creation", T::live_instance_count, cnt+1);
  
  pointer_t ptr2 = ptr1;
  
  tester.test("use count after sharing", ptr1.use_count(), 2L);
  tester.test("instance count after sharing", T::live_instance_count, cnt+1);
  ptr2.reset();
  tester.test("use count after reset2", ptr1.use_count(), 1L);
  tester.test("instance count after reset1", T::live_instance_count, cnt+1);
  ptr1.reset();
  tester.test("use count after reset1", ptr1.use_count(), 0L);
  tester.test("instance count after reset2", T::live_instance_count, cnt);
    
  failed_test_count += tester.failed_test_count();
}

int main() {

  test_shared_pointer("shared_pointer<not_intrusive>", test_not_intrusive(), sizeof(void*)+sizeof(long*));
  test_shared_pointer("shared_pointer<intrusive>", test_intrusive(), sizeof(void*));

  return (int)failed_test_count;
}
