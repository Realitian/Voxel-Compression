#+++HDR+++
#======================================================================
#  This file is part of the SL software library.
#
#  Copyright (C) 1993-2014 by Enrico Gobbetti (gobbetti@crs4.it)
#  Copyright (C) 1996-2014 by CRS4 Visual Computing Group, Pula, Italy
#
#  For more information, visit the CRS4 Visual Computing Group 
#  web pages at http://www.crs4.it/vvr/.
#
#  This file may be used under the terms of the GNU General Public
#  License as published by the Free Software Foundation and appearing
#  in the file LICENSE included in the packaging of this file.
#
#  CRS4 reserves all rights not expressly granted herein.
#  
#  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
#  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
#  FOR A PARTICULAR PURPOSE.
#
#======================================================================
#---HDR---#
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sl)

#################################################
# The version number

SET(SL_VERSION_MAJOR 7)
SET(SL_VERSION_MINOR 8)
SET(SL_VERSION_MICRO 2)
SET(SL_VERSION "${SL_VERSION_MAJOR}.${SL_VERSION_MINOR}.${SL_VERSION_MICRO}")

#################################################
# Default settings
IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF (NOT CMAKE_BUILD_TYPE)

MESSAGE( "-----------------------------------------------------------------" )
MESSAGE( STATUS "Configuring ${PROJECT_NAME} ${SL_VERSION} in ${CMAKE_BUILD_TYPE} mode" )
MESSAGE( "-----------------------------------------------------------------" )

#################################################
# Options

OPTION(SL_DOC "Build HTML doc (via Doxygen)" OFF)
OPTION(SL_TEST "Build testsuite & benchmark apps" ON)
OPTION(SL_INSTALL_SLFIND_MODULE "Install FindSL.cmake into cmake modules" ON)

#################################################
# Include

INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckTypeSize)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckLibraryExists)

#################################################
# Macro

MACRO (cxx_check_include_file arg result)
  CHECK_CXX_SOURCE_COMPILES("#include <${arg}>
                             int main(void)
                             {
                                return 0;
                             }
                            " ${result})
ENDMACRO (cxx_check_include_file)

MACRO (check_function arg result)
  CHECK_CXX_SOURCE_COMPILES("#include <iostream>
                             int main(void)
                             {
                                std::cout << ${arg};
                             }
                            " ${result})
ENDMACRO (check_function)

MACRO (check_restrict arg result)
  CHECK_CXX_SOURCE_COMPILES("void f (int * ${arg} x, int * ${arg} y) { *x = *y; }
                             int main(void)
                             {
                               int x=0;
                               int y=1;
                               f(&x, &y);
                             }
                            " ${result})
ENDMACRO (check_restrict)

MACRO (sizeof arg have_result have_result_bits)
  CHECK_TYPE_SIZE(${arg} ${have_result_bits})
  IF (HAVE_${have_result_bits})
    MATH( EXPR ${have_result_bits} "${${have_result_bits}}*8" )
    SET(${have_result} 1)
  ENDIF (HAVE_${have_result_bits})
ENDMACRO (sizeof)

#################################################
# Configuration 

INCLUDE(TestForSTDNamespace) 
IF (NOT CMAKE_NO_STD_NAMESPACE)
  SET(HAVE_NAMESPACE_STD 1)
ENDIF (NOT CMAKE_NO_STD_NAMESPACE)

INCLUDE(TestForANSIForScope) 
IF (NOT CMAKE_NO_ANSI_FOR_SCOPE)
  SET(HAVE_ANSI_FOR_SCOPE 1)
ENDIF (NOT CMAKE_NO_ANSI_FOR_SCOPE)

sizeof("short" HAVE_SHORT HAVE_SHORT_BITS)
sizeof("int" HAVE_INT HAVE_INT_BITS)
sizeof("long" HAVE_LONG HAVE_LONG_BITS)
sizeof("long long" HAVE_LONG_LONG HAVE_LONG_LONG_BITS)
sizeof("long double" HAVE_LONG_DOUBLE HAVE_LONG_DOUBLE_BITS)
sizeof("size_t" HAVE_SIZE_T HAVE_SIZE_T_BITS)

check_function(__FUNCTION__ HAVE_FUNCTION)
check_function(__func__ HAVE_PRETTY_FUNC)
check_function(__PRETTY_FUNCTION__ HAVE_PRETTY_FUNCTION)

check_restrict(restrict HAVE_RESTRICT)
check_restrict(__restrict HAVE_PRIVATE_RESTRICT)
check_restrict(__restrict__ HAVE_PRIVATE_SGI_RESTRICT)

CHECK_CXX_SOURCE_COMPILES("#include <unordered_map>
                           #include <unordered_set>
			     int main(void)
                             {
                                std::unordered_set<int> my_set;
                                std::unordered_map<int,int> my_map;
				return 0;
                             }
                             " HAVE_STL_UNORDERED_CONTAINERS)
CHECK_CXX_SOURCE_COMPILES("#include <tr1/unordered_map>
                           #include <tr1/unordered_set>
			     int main(void)
                             {
                                std::tr1::unordered_set<int> my_set;
                                std::tr1::unordered_map<int,int> my_map;
				return 0;
                             }
                             " HAVE_STL_TR1_UNORDERED_CONTAINERS)
CHECK_CXX_SOURCE_COMPILES("#include <unordered_map>
                           #include <unordered_set>
                             int main(void)
                             {
                                std::tr1::unordered_set<int> my_set;
                                std::tr1::unordered_map<int,int> my_map;
				return 0;
                             }
                             " HAVE_STL_TR1_UNORDERED_CONTAINERS_MSVC)

CHECK_LIBRARY_EXISTS(m random "" HAVE_LIBM) 
CHECK_LIBRARY_EXISTS(m lrint "" HAVE_LRINT) 
CHECK_LIBRARY_EXISTS(m lrintf "" HAVE_LRINTF) 

# Threads
FIND_PACKAGE(Threads)
IF(CMAKE_USE_WIN32_THREADS_INIT) 
  SET(HAVE_WIN32_THREADS 1)
  MESSAGE("-- Using win32 threads")
ELSE(CMAKE_USE_WIN32_THREADS_INIT)
  IF(CMAKE_USE_PTHREADS_INIT)
    SET(HAVE_POSIX_THREADS 1)
    MESSAGE("-- Using POSIX threads")
  ELSE(CMAKE_USE_PTHREADS_INIT)
    MESSAGE("-- Unsupported thread library - threads are disabled")
  ENDIF(CMAKE_USE_PTHREADS_INIT)
ENDIF(CMAKE_USE_WIN32_THREADS_INIT)

# ZLIB
FIND_PACKAGE(ZLIB)
IF (ZLIB_FOUND)
    SET(HAVE_ZLIB 1)
    MESSAGE("-- Using ZLIB")
ELSE(ZLIB_FOUND)
    MESSAGE("-- No zlib")
ENDIF(ZLIB_FOUND)

SET(conf_dir ${CMAKE_CURRENT_SOURCE_DIR}/src/sl/config)
SET(conf_file_in ${conf_dir}/config.h.in)
SET(conf_file ${conf_dir}/config.h)

MESSAGE("-- Generating Config file: ${conf_file}")
CONFIGURE_FILE(${conf_file_in}
               ${conf_file})

# Adding a dummy target is to allow ADDITIONAL_MAKE_CLEAN_FILES
ADD_CUSTOM_TARGET(
 config
 COMMAND ${CMAKE_COMMAND} -E echo "Generating Config file..."
 WORKING ${conf_dir}
)
SET_DIRECTORY_PROPERTIES(PROPERTIES
                         ADDITIONAL_MAKE_CLEAN_FILES
                         "${conf_file};")
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall-cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#################################################
# Compiler dependent definitions


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
    MESSAGE("-- Current build compiler is Clang")
    ADD_DEFINITIONS(-D_ISOC9X_SOURCE) 
    ADD_DEFINITIONS(-D_ISOC99_SOURCE)
    ADD_DEFINITIONS(-D__USE_ISOC9X)
    ADD_DEFINITIONS(-D__USE_ISOC99)
    # Compiler flags
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall"  )
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall"  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
    MESSAGE("-- Current build compiler is g++")

    ##define _ISOC9X_SOURCE 1
    ADD_DEFINITIONS(-D_ISOC9X_SOURCE) 
    ADD_DEFINITIONS(-D_ISOC99_SOURCE)
    ADD_DEFINITIONS(-D__USE_ISOC9X)
    ADD_DEFINITIONS(-D__USE_ISOC99)

    # Compiler flags
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall"  )
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall"  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
    MESSAGE("-- Current build compiler is Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
    MESSAGE("-- Current build compiler is Microsoft Visual Studio")
    ADD_DEFINITIONS(-D_ISOC9X_SOURCE) 
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
endif()

#################################################
# Enable testing

IF (SL_TEST)
  ENABLE_TESTING()
ENDIF (SL_TEST)

#################################################
# make package package_source

SET(CPACK_GENERATOR "TGZ")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES  "CMakeCache.txt;/CMakeFiles/;/Makefile/;/Testing/;
                                /_CPack_Packages/;/CVS/;/\\\\.svn/;\\\\.swp$;
                                \\\\.#;/#;.*~;cscope.*;.*\\\\.tar\\\\.gz;
                                .*\\\\.tar\\\\.bz2;.*\\\\.tar\\\\.Z")
SET(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
SET(CPACK_PACKAGE_VENDOR "CRS4/ViC")
SET(CPACK_PACKAGE_VERSION_MAJOR ${SL_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${SL_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${SL_VERSION_MICRO})
SET(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/COPYING)
SET(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README)
INCLUDE(CPack)

#################################################
# Print status
MESSAGE( "-- CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )

#################################################
# Subdirs

MESSAGE( "-----------------------------------------------------------------" )

ADD_SUBDIRECTORY(src)

IF (SL_INSTALL_SLFIND_MODULE)
   MESSAGE( "-- SL_INSTALL_SLFIND_MODULE = ON" )
   ADD_SUBDIRECTORY(cmake_modules)
ENDIF (SL_INSTALL_SLFIND_MODULE)

IF (SL_TEST)
  ADD_SUBDIRECTORY(testsuite)
  ADD_SUBDIRECTORY(benchmark)
ENDIF (SL_TEST)

IF (SL_DOC)
   ADD_SUBDIRECTORY(doc)
ENDIF (SL_DOC)

MESSAGE( "-----------------------------------------------------------------" )

