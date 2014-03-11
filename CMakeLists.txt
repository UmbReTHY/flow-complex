project(flow_complex)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

# TODO detect clang and gcc and restrict flags to these two
set (CMAKE_CXX_FLAGS "--std=c++11")
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -O0")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")

option(ENABLE_TESTS
       "if set, builds the tests and provides a test target for Makefiles"
       OFF)

find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
     
if (ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif (ENABLE_TESTS)
