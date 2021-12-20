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
#ifndef SL_ASSERT_HPP
#define SL_ASSERT_HPP

#include <sl/config.hpp>
#include <iostream>

// ---------------------------------------------------------------------------
// Trace levels
// ---------------------------------------------------------------------------

/** 
 * \def SL_TRACE_LEVEL
 * trace level (integer): 0 = no trace, 1 = minimal trace, 2+ = more verbosity.
 * The defaule is 0 for no debugging and 1 for debugging mode
 */

#ifndef SL_TRACE_LEVEL
#ifdef NDEBUG
#  define   SL_TRACE_LEVEL 0
#else 
#  define   SL_TRACE_LEVEL 1
#endif
#endif

/** 
 * \def SL_TRACE_OUT(X)
 * The command SL_TRACE_OUT(X) << msg outputs msg to std::cerr if X is 
 * lower than the current trace level.
 */

#define SL_TRACE_OUT(X) if ((X) < (SL_TRACE_LEVEL)) \
  std::cerr << CXX_CURRENT_FUNCTION << " - TRACE(" << (X) << "): " 

// ---------------------------------------------------------------------------
// Assertion levels
// ---------------------------------------------------------------------------

/** 
 * \def SL_CHECK_NO
 * assertion check level: no assertion checking
 */
/** 
 * \def SL_CHECK_MINIMAL
 * assertion check level: only method preconditions, no stack frames
 */
/** 
 * \def SL_CHECK_REQUIRE
 * assertion check level: only method preconditions
 */
/** 
 * \def SL_CHECK_ENSURE
 * assertion check level: require+ensure, i.e. method pre/postconditions
 */
/** 
 * \def SL_CHECK_INVARIANT
 * assertion check level: as ensure plus check class invariant on entry/exit
 */
/** 
 * \def SL_CHECK_LOOP
 * assertion check level: as invariant plus loop invariants + variants
 */
/** 
 * \def SL_CHECK_ALL
 * assertion check level: everything including check statements
 */

#define  SL_CHECK_NO        0
#define  SL_CHECK_MINIMAL   1
#define  SL_CHECK_REQUIRE   2
#define  SL_CHECK_ENSURE    3
#define  SL_CHECK_INVARIANT 4
#define  SL_CHECK_LOOP      5
#define  SL_CHECK_ALL       6

/** 
 * \def SL_ASSERT_LEVEL
 * assertion check level
 */

#ifndef SL_ASSERT_LEVEL
#  ifdef  NDEBUG
#    define SL_ASSERT_LEVEL     SL_CHECK_NO
#  else
#    define SL_ASSERT_LEVEL     SL_CHECK_ALL
#  endif
#endif

// ---------------------------------------------------------------------------
// Behavior upon failure
// ---------------------------------------------------------------------------


#define SL_ABORT    true
                             // Abort
#define SL_CONTINUE false 
                             // Log and continue

/** 
 * \def SL_BEHAVIOR_UPON_FAILURE
 * assertion behavior upon failure: if true, abort, if false log and continue
 */
#ifndef SL_BEHAVIOR_UPON_FAILURE
#define SL_BEHAVIOR_UPON_FAILURE SL_ABORT
#endif
 

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

namespace sl {
  namespace detail {
    extern void failure(const char* tp,
		        const char* tag,
		        const char* expr,
		        const char* file,
		        int line,
		        const char* fname,
			bool do_abort);
  }
};


#define SL_ASSERT_IMPL(tp,tag, e, abrt) \
  ((void) ((e) ? 0 : (sl::detail::failure((tp), (tag), (#e), __FILE__, __LINE__, CXX_CURRENT_FUNCTION,(abrt)), 0)))

// ---------------------------------------------------------------------------
// Eiffel-style Assertion Macros
// ---------------------------------------------------------------------------

/** 
 * \def SL_COMPILE_TIME_CHECK(tag, e)
 * fail if compile-time constant e evaluates to false
 */
#define SL_COMPILE_TIME_CHECK(tag, e) \
  SL_COMPILE_TIME_CHECK_IMPL(tag, e,__INCLUDE_LEVEL__,__LINE__)
 
#define SL_COMPILE_TIME_CHECK_IMPL(tag, e, il, ll) \
  SL_COMPILE_TIME_CHECK_IMPL2(tag, e, il, ll)


#ifdef USE_TYPEDEF_FOR_COMPILE_TIME_CHECK

#define SL_COMPILE_TIME_CHECK_IMPL2(tag, e, il, ll) \
  typedef bool SL_COMPILE_TIME_CHECK_ ## il  ## _ ## ll  [(e) ? 1 : -1]

#else 

#define SL_COMPILE_TIME_CHECK_IMPL2(tag, e, il, ll) \
   enum { SL_COMPILE_TIME_CHECK_ ## il  ## _ ## ll \
      = sizeof(::sl::detail::STATIC_ASSERTION_FAILURE< ( e ) >) }

namespace sl{
  namespace detail {
    template <bool> struct STATIC_ASSERTION_FAILURE;
    template <> struct STATIC_ASSERTION_FAILURE<true>{};
  }
}

#endif

/** 
 * \def SL_REQUIRE(tag, e)
 * precondition check, fail if e evaluates to false
 */
#if (SL_ASSERT_LEVEL >= SL_CHECK_MINIMAL) 
#define SL_REQUIRE(tag, e) \
    SL_ASSERT_IMPL("Precondition", tag, e, SL_BEHAVIOR_UPON_FAILURE)
#else 
#define SL_REQUIRE(tag, e)
#endif

/** 
 * \def SL_ENSURE(tag, e)
 * postcondition check, fail if e evaluates to false
 */
#if (SL_ASSERT_LEVEL >= SL_CHECK_ENSURE)
#define SL_ENSURE(tag, e) \
    SL_ASSERT_IMPL("Postcondition", tag, e, SL_BEHAVIOR_UPON_FAILURE)
#else 
#define SL_ENSURE(tag, e)
#endif

/** 
 * \def SL_INVARIANT(e)
 * invariant check, fail if e evaluates to false
 */
#if (SL_ASSERT_LEVEL >= SL_CHECK_INVARIANT)
#define SL_INVARIANT(e) \
    SL_ASSERT_IMPL("Invariant", "INV", e, SL_BEHAVIOR_UPON_FAILURE)
#else 
#define SL_INVARIANT(e)
#endif

/** 
 * \def SL_LOOP_INVARIANT(tag, e)
 * loop invariant check, fail if e evaluates to false
 */
#if (SL_ASSERT_LEVEL >= SL_CHECK_LOOP)
#define SL_LOOP_INVARIANT(tag, e) \
    SL_ASSERT_IMPL("Loop invariant",tag, e, SL_BEHAVIOR_UPON_FAILURE)
#else 
#define SL_LOOP_INVARIANT(tag, e)
#endif

/** 
 * \def SL_CHECK(tag, e)
 * assertion check, fail if e evaluates to false
 */
#if (SL_ASSERT_LEVEL >= SL_CHECK_ALL)
#define SL_CHECK(tag, e) \
    SL_ASSERT_IMPL("Assertion",tag, e, SL_BEHAVIOR_UPON_FAILURE)
#else 
#define SL_CHECK(tag, e)
#endif

#if (SL_ASSERT_LEVEL >= SL_CHECK_ALL)
#define SL_CHECK_WARNING(tag, e) \
    SL_ASSERT_IMPL("Warning", tag, e, false)
#else 
#define SL_CHECK_WARNING(tag, e)
#endif

/** 
 * \def SL_FAIL(tag)
 * alwais fails
 */
#define SL_FAIL(tag) \
    sl::detail::failure("Failure", (tag), "fail request", __FILE__, __LINE__, CXX_CURRENT_FUNCTION, true)

#endif
