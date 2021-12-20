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
#include <sl/tester.hpp>

namespace sl {

  tester::tester()
    : 
    name_("Anonymous"), last_test_failed_(false), failed_test_count_(0), passed_test_count_(0), silent_(false) {
  }

  tester::tester(const std::string& n): 
    name_(n),  last_test_failed_(false), failed_test_count_(0), passed_test_count_(0), silent_(false) {
  }

  tester::~tester() {
    if (is_testing()) {
      report();
      reset();
    }
  }

  bool tester::is_silent() const {
    return silent_;
  }

  void tester::set_silent(bool b) {
    silent_ = b;
  }

  bool tester::has_last_test_failed() const {
    return last_test_failed_;
  }

  size_t tester::failed_test_count() const {
    return failed_test_count_;
  }
  size_t tester::passed_test_count() const {
    return passed_test_count_;
  }
  
  size_t tester::current_test_count() const {
    return failed_test_count() + passed_test_count();
  }

  bool tester::is_testing() const {
    return current_test_count() > 0;
  }

  const std::string& tester::name() const {
    return name_;
  }

  void tester::reset() {
    failed_test_count_ = 0;
    passed_test_count_ = 0;
    last_test_failed_ = false;
  }

  void tester::set_name(const std::string& n) {
    SL_REQUIRE("Not testing", ! is_testing());
    name_ = n;
    SL_ENSURE("New name", name_ == n);
  }

  void tester::test(const std::string& id, 
		    bool actual) {
    const bool expected = true;
    test(id, actual, expected);
  }

  void tester::report() {
    if (!is_silent()) {
      if (current_test_count() > 0) { 
	if (failed_test_count() == 0) {
	  std::cerr << "### [Passed " << passed_test_count() << "/" << passed_test_count() << "]";
	} else {
	  std::cerr << "### [FAILED " << failed_test_count() << "/" << (failed_test_count() + passed_test_count()) << "]";
	}
	std::cerr << ": Test suite ``" << name() << "''" << std::endl;
	std::cerr << "### ------------------------------------------------------" << std::endl;
      }
    }
  }

} // namespace sl
