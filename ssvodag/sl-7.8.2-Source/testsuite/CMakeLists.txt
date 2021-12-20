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
PROJECT(TestSuite)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

ADD_EXECUTABLE(test_argument_parser test_argument_parser.cpp)
ADD_DEPENDENCIES(test_argument_parser sl)
TARGET_LINK_LIBRARIES(test_argument_parser sl)
ADD_TEST(test_argument_parser test_argument_parser)

ADD_EXECUTABLE(test_bitops test_bitops.cpp)
ADD_DEPENDENCIES(test_bitops sl)
TARGET_LINK_LIBRARIES(test_bitops sl)
ADD_TEST(test_bitops test_bitops)

ADD_EXECUTABLE(test_bounded_scalar test_bounded_scalar.cpp)
ADD_DEPENDENCIES(test_bounded_scalar sl)
TARGET_LINK_LIBRARIES(test_bounded_scalar sl)
ADD_TEST(test_bounded_scalar test_bounded_scalar)

ADD_EXECUTABLE(test_color test_color.cpp)
ADD_DEPENDENCIES(test_color sl)
TARGET_LINK_LIBRARIES(test_color sl)
ADD_TEST(test_color test_color)

ADD_EXECUTABLE(test_compression test_compression.cpp)
ADD_DEPENDENCIES(test_compression sl)
TARGET_LINK_LIBRARIES(test_compression sl)
ADD_TEST(test_compression test_compression)

ADD_EXECUTABLE(test_cubature test_cubature.cpp)
ADD_DEPENDENCIES(test_cubature sl)
TARGET_LINK_LIBRARIES(test_cubature sl)
ADD_TEST(test_cubature test_cubature)

ADD_EXECUTABLE(test_external_array test_external_array.cpp)
ADD_DEPENDENCIES(test_external_array sl)
TARGET_LINK_LIBRARIES(test_external_array sl)
ADD_TEST(test_external_array test_external_array)

ADD_EXECUTABLE(test_fem test_fem.cpp)
ADD_DEPENDENCIES(test_fem sl)
TARGET_LINK_LIBRARIES(test_fem sl)
ADD_TEST(test_fem test_fem)

ADD_EXECUTABLE(test_fixed_unit_real test_fixed_unit_real.cpp)
ADD_DEPENDENCIES(test_fixed_unit_real sl)
TARGET_LINK_LIBRARIES(test_fixed_unit_real sl)
ADD_TEST(test_fixed_unit_real test_fixed_unit_real)

ADD_EXECUTABLE(test_floating_cone test_floating_cone.cpp)
ADD_DEPENDENCIES(test_floating_cone sl)
TARGET_LINK_LIBRARIES(test_floating_cone sl)
ADD_TEST(test_floating_cone test_floating_cone)

ADD_EXECUTABLE(test_geometry test_geometry.cpp)
ADD_DEPENDENCIES(test_geometry sl)
TARGET_LINK_LIBRARIES(test_geometry sl)
ADD_TEST(test_geometry test_geometry)

ADD_EXECUTABLE(test_graph_minimum_linear_arranger test_graph_minimum_linear_arranger.cpp)
ADD_DEPENDENCIES(test_graph_minimum_linear_arranger sl)
TARGET_LINK_LIBRARIES(test_graph_minimum_linear_arranger sl)
ADD_TEST(test_graph_minimum_linear_arranger test_graph_minimum_linear_arranger)

ADD_EXECUTABLE(test_absolute_orientation test_absolute_orientation.cpp)
ADD_DEPENDENCIES(test_absolute_orientation sl)
TARGET_LINK_LIBRARIES(test_absolute_orientation sl)
ADD_TEST(test_absolute_orientation test_absolute_orientation)

ADD_EXECUTABLE(test_indexed test_indexed.cpp)
ADD_DEPENDENCIES(test_indexed sl)
TARGET_LINK_LIBRARIES(test_indexed sl)
ADD_TEST(test_indexed test_indexed)

ADD_EXECUTABLE(test_interval test_interval.cpp)
ADD_DEPENDENCIES(test_interval sl)
TARGET_LINK_LIBRARIES(test_interval sl)
ADD_TEST(test_interval test_interval)

ADD_EXECUTABLE(test_knn test_knn.cpp)
ADD_DEPENDENCIES(test_knn sl)
TARGET_LINK_LIBRARIES(test_knn sl)
ADD_TEST(test_knn test_knn)

ADD_EXECUTABLE(test_matrix test_matrix.cpp)
ADD_DEPENDENCIES(test_matrix sl)
TARGET_LINK_LIBRARIES(test_matrix sl)
ADD_TEST(test_matrix test_matrix)

ADD_EXECUTABLE(test_miniball test_miniball.cpp)
ADD_DEPENDENCIES(test_miniball sl)
TARGET_LINK_LIBRARIES(test_miniball sl)
ADD_TEST(test_miniball test_miniball)

ADD_EXECUTABLE(test_normal_compressor test_normal_compressor.cpp)
ADD_DEPENDENCIES(test_normal_compressor sl)
TARGET_LINK_LIBRARIES(test_normal_compressor sl)
ADD_TEST(test_normal_compressor test_normal_compressor)

ADD_EXECUTABLE(test_random test_random.cpp)
ADD_DEPENDENCIES(test_random sl)
TARGET_LINK_LIBRARIES(test_random sl)
ADD_TEST(test_random test_random)

ADD_EXECUTABLE(test_serializer test_serializer.cpp)
ADD_DEPENDENCIES(test_serializer sl)
TARGET_LINK_LIBRARIES(test_serializer sl)
ADD_TEST(test_serializer test_serializer)

ADD_EXECUTABLE(test_smart_pointer test_smart_pointer.cpp)
ADD_DEPENDENCIES(test_smart_pointer sl)
TARGET_LINK_LIBRARIES(test_smart_pointer sl)
ADD_TEST(test_smart_pointer test_smart_pointer)

ADD_EXECUTABLE(test_tester test_tester.cpp)
ADD_DEPENDENCIES(test_tester sl)
TARGET_LINK_LIBRARIES(test_tester sl)
ADD_TEST(test_tester test_tester)

ADD_EXECUTABLE(test_threads test_threads.cpp)
ADD_DEPENDENCIES(test_threads sl)
TARGET_LINK_LIBRARIES(test_threads sl)
ADD_TEST(test_threads test_threads)

ADD_EXECUTABLE(test_time test_time.cpp)
ADD_DEPENDENCIES(test_time sl)
TARGET_LINK_LIBRARIES(test_time sl)
ADD_TEST(test_time test_time)

ADD_EXECUTABLE(test_triangle_mesh test_triangle_mesh.cpp)
ADD_DEPENDENCIES(test_triangle_mesh sl)
TARGET_LINK_LIBRARIES(test_triangle_mesh sl)
ADD_TEST(test_triangle_mesh test_triangle_mesh)

ADD_EXECUTABLE(test_containers test_containers.cpp)
ADD_DEPENDENCIES(test_containers sl)
TARGET_LINK_LIBRARIES(test_containers sl)
ADD_TEST(test_containers test_containers)
