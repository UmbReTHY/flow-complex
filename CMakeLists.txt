project(flow_complex)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

# TODO detect clang and gcc and restrict flags to these two
set (CMAKE_CXX_FLAGS "-std=c++11 -fno-exceptions -fno-rtti")
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -O0")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")

option(BUILD_TESTS
       "builds the tests and provide a test target for Makefiles"
       OFF)

option(BUILD_EXAMPLES "" OFF)

find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
     
if (BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif (BUILD_TESTS)

if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif (BUILD_EXAMPLES)

