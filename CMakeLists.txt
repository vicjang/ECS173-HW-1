#This software is property of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  Do not distribute or copy this software without the consent of Owen Carmichael, Chris Schwarz, and the Regents of the University of California.  All rights reserved.  Copyright ©  2007-2013  Owen Carmichael, Chris Schwarz, and the Regents of the University of California.
cmake_minimum_required(VERSION 2.4)

# This project is designed to be built outside the Insight source tree.
PROJECT(eigen_faces)

# include information to allow running this in a debugger
# change "Debug" to "Release" or remove this line entirely to remove
SET(CMAKE_BUILD_TYPE Debug)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(eigen_faces eigen_faces.cxx)

TARGET_LINK_LIBRARIES(eigen_faces ${ITK_LIBRARIES})
