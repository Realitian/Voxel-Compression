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
#ifndef SL_TESTER_HPP
#define SL_TESTER_HPP

#include <sl/interval.hpp>
#include <string>
#include <iostream>
#include <functional>


namespace sl {

  /**
   *  A class for implementing unit test suites.
   *  Individual tests are implemented by creating
   *  a 'tester' objects and by applying to it
   *  the feature 'test' for each of the
   *  individual tests. 
   */
  class tester {
  protected: 
    
    std::string name_;

    bool   last_test_failed_;
    size_t failed_test_count_;
    size_t passed_test_count_;

    bool silent_;

  public: // Creation / Copy / Destruction

    /// Default creation.
    tester();

    /// Initialize a tester objects named \e n
    tester(const std::string& n);

    /// Destruct test, eventually printing a report if needed.
    ~tester();

  public: // Queries

    /// Is the tester not dumping results on std::cerr?
    bool is_silent() const;

    /**
     *  If b is false, the tester dumps result on std::cerr, otherwise it will
     *  not write anything on output
     */
    void set_silent(bool b);

    /// The number of tests failed since the last reset.
    size_t failed_test_count() const;

    /// The number of tests passed since the last reset.
    size_t passed_test_count() const;

    /// The number of tests executed since tha last reset.
    size_t current_test_count() const;

    /// True if at least one test has been executed.
    bool is_testing() const;

    /// The name of the current test suite.
    const std::string& name() const;

    /// Did the last test fail?
    bool has_last_test_failed() const;

  public: // Commands
    
    /// Reset to zero the number of tests passed/failed.
    void reset();

    /// Set the name of the current test suite to \e n.
    void set_name(const std::string& n);

  public: // Testing
    
    /** Check result of a test, and update test statistics.
     *  @param id       the test name.
     *  @param actual   the actual result.
     *  @param expected the expected result.
     *  @param eq       the function for comparing actual to expected.
     */
    template <class T>
    void test(const std::string& id, 
	      const T& actual, 
	      const T& expected, 
	      const std::equal_to<T>& eq);

    /** Check result of a test, and update test statistics.
     *  @param id       the test name.
     *  @param actual   the actual result.
     *  @param expected the expected result.
     *  Actual and expected are compared with std::equal_to<T>
     */
    template <class T>
    void test(const std::string& id, 
	      const T& actual, 
	      const T& expected);

    /** Check result of a test, and update test statistics.
     *  @param id       the test name.
     *  @param actual   the actual result.
     *  @param expected the expected result range.
     */
    template <class T>
    void test(const std::string& id, 
	      const T& actual, 
	      const interval<T>& expected);

    /** Consider a test passed, and update test statistics.
     *  @param id       the test name.
     *  @param actual   the actual result.
     */
    template <class T>
    void unchecked_test(const std::string& id, 
			const T& actual);

    /** Check result of a test, and update test statistics.
     *  @param id       the test name.
     *  @param actual   the actual result (compared to true)
     */
     void test(const std::string& id, 
	       bool actual);

  public: // Reporting
    
    /// Report to std::cerr the results of last test.
    template <class T_ACTUAL, class T_EXPECTED>
    void report_last_test_result(const std::string& id, 
				 const T_ACTUAL& actual, 
				 const T_EXPECTED& expected,
				 bool passed);

    /// Summarize to std::cerr the results of the tests executed so far.
    void report();

  }; // class tester

} // namespace sl

//-----------------------------------------------------------------------
// --- Template implementation
//-----------------------------------------------------------------------

namespace sl {

  namespace detail {

    template <class T>
    inline void test_dispatcher(sl::tester& tst,
			 const std::string& id,
			 T const& actual, 
			 T const& expected) {
      std::equal_to<T> eq;
      tst.test(id, actual, expected, eq);
    }

    template <>
    inline void test_dispatcher(sl::tester& tst,
				const std::string& id,
				const char* const& actual, 
				const char* const& expected) {
      std::equal_to<std::string> eq;
      tst.test(id, 
	       std::string(actual), 
	       std::string(expected), eq);
    }

  } // namespace detail

  /// Report to std::cerr the results of last test.
  template <class T_ACTUAL, class T_EXPECTED>
  void tester::report_last_test_result(const std::string& id, 
				       const T_ACTUAL& actual, 
				       const T_EXPECTED& expected,
				       bool passed) {
    if (!is_silent()) {
      if (current_test_count() == 1) {
	std::cerr << "### ------------------------------------------------------" << std::endl;
	std::cerr << "### Test suite ``" << name() << "'' started." << std::endl;
	std::cerr << "###" << std::endl;
      }
      
      std::cerr <<
	(passed ? "---" : "!!!") <<
	" [" << current_test_count() << "]: " <<
	(passed ? "Passed" : "FAILED") << " " <<
	id << std::endl;
      if (!passed) {
	std::cerr << "  Actual = " << std::endl << actual << std::endl;
	std::cerr << "  Expected = " << std::endl << expected << std::endl;
      }
    }
  }

  template <class T>
  void tester::test(const std::string& id, 
		    const T& actual, 
		    const T& expected, 
		    const std::equal_to<T>& eq) {
    bool passed = eq(actual,expected); 
    
    if (passed) {
      passed_test_count_++;
    } else {
      failed_test_count_++;
    }
    last_test_failed_ = !passed;

    report_last_test_result(id, actual, expected, passed);
  }

  template <class T>
  void tester::test(const std::string& id, 
		    const T& actual, 
		    const T& expected) {
    detail::test_dispatcher(*this, id, actual, expected);
  }

  template <class T>
  void tester::unchecked_test(const std::string& id, 
			      const T& actual) {
    bool passed = true;
    passed_test_count_++;
    last_test_failed_ = !passed;

    report_last_test_result(id, actual, actual, passed);
  }

  template <class T>
  void tester::test(const std::string& id, 
		    const T& actual, 
		    const interval<T>& expected) {
    bool passed = expected.contains(actual);
    
    if (passed) {
      passed_test_count_++;
    } else {
      failed_test_count_++;
    }
    last_test_failed_ = !passed;
   
    report_last_test_result(id, actual, expected, passed);
  }
  
} // namespace sl

#endif
