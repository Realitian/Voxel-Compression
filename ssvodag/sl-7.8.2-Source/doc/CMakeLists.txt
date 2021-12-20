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
#################################################
# Doc

INCLUDE(FindDoxygen)
IF (${DOXYGEN_FOUND} STREQUAL "YES")
  SET(PACKAGE_NAME ${CMAKE_PROJECT_NAME})
  SET(VERSION ${SL_VERSION})
  SET(TOP_SRCDIR "${CMAKE_SOURCE_DIR}")
  CONFIGURE_FILE(Doxyfile.in "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")

  SET(HTML_TARGET "doc" ) 
  ADD_CUSTOM_TARGET(${HTML_TARGET}
                    ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
                    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile") 
  INSTALL( DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html/" DESTINATION "share/doc/${PACKAGE_NAME}" )
  SET_DIRECTORY_PROPERTIES(PROPERTIES
                           ADDITIONAL_MAKE_CLEAN_FILES
                           "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile;${CMAKE_CURRENT_BINARY_DIR}/html;${CMAKE_CURRENT_BINARY_DIR}/latex")

ENDIF (${DOXYGEN_FOUND} STREQUAL "YES")


