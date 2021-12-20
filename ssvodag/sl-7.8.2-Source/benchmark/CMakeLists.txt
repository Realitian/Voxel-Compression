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
PROJECT(benchmark)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

SET(BENCHMARK_APPS bench_hull bench_obox 
                   fixed_size_daxpy bench_knn
                   bench_compression bench_external_array
		   bench_geometry bench_containers bench_minla)
#                   dense_daxpy)

FOREACH(APP ${BENCHMARK_APPS})
  ADD_EXECUTABLE(sl_${APP} ${APP}.cpp)
  ADD_DEPENDENCIES(sl_${APP} sl)
  TARGET_LINK_LIBRARIES(sl_${APP} sl)
  INSTALL(TARGETS sl_${APP} DESTINATION bin)
ENDFOREACH(APP ${BENCHMARK_APPS})

