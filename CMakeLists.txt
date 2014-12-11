cmake_minimum_required(VERSION 2.8)
project(VShape)

set( VSHAPE_VERSION_MAJOR 0 )
set( VSHAPE_VERSION_MINOR 1 )
set( VSHAPE_VERSION_BUILD 0 )

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type specified. Using 'release'")
	SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/vshape-${VSHAPE_VERSION_MAJOR}.${VSHAPE_VERSION_MINOR}.${VSHAPE_VERSION_BUILD}/)

include_directories (libs)
include_directories (libs/bamtools-2.3.0/src)
include_directories (libs/zlib-1.2.7)

add_subdirectory(libs/zlib-1.2.7)
add_subdirectory(libs/bamtools-2.3.0)
add_subdirectory(src)