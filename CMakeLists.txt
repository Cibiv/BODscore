cmake_minimum_required(VERSION 2.8)
project(VShape)

set( VSHAPE_VERSION_MAJOR 1 )
set( VSHAPE_VERSION_MINOR 0 )
set( VSHAPE_VERSION_BUILD 0 )

# Set a default build type for single-configuration
# CMake generators if no build type is set.

#IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
#	message(STATUS "No build type specified. Using 'release'")
#	SET(CMAKE_BUILD_TYPE Release)
# ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

#####################################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
#####################################

SET(CMAKE_BUILD_TYPE DEBUG)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/vshape-${VSHAPE_VERSION_MAJOR}.${VSHAPE_VERSION_MINOR}.${VSHAPE_VERSION_BUILD}/)

include_directories (libs)

# set(sqlite_path "${CMAKE_CURRENT_SOURCE_DIR}/libs/sqlite-amalgamation-3080802")
# include_directories(sqlite_path)
# add_subdirectory( sqlite_path )
# add_library(sqlite3 STATIC "${CMAKE_CURRENT_SOURCE_DIR}/libs/sqlite-amalgamation-3080802/sqlite3.c")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/libs/sqlite3pp/src")
add_subdirectory("libs/sqlite3pp/src")
# add_library(sqlite3pp SHARED "${CMAKE_CURRENT_SOURCE_DIR}/libs/sqlite3pp/src")

# TARGET_LINK_LIBRARIES(sqlite3pp "${CMAKE_CURRENT_SOURCE_DIR}/libs/sqlite-amalgamation-3080802/sqlite3.c" )
# target_include_directories( sqlite3pp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libs/sqlite-amalgamation-3080802> )
#######################################

include_directories (libs/bamtools-2.3.0/src)
add_subdirectory(libs/bamtools-2.3.0)

include_directories (libs/zlib-1.2.7)
add_subdirectory(libs/zlib-1.2.7)

#######################################

add_subdirectory(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0 -Wall ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g3 -O0 -Wall ")

SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}-g3 -O0 -Wall -std=c++0x -pedantic ")
SET (CMAKE_C_FLAGS_DEBUG "-g3 -O0 -Wall ")

#######################################
# install (TARGETS VShape DESTINATION bin)
# install (FILES sqlite.c DESTINATION include)
