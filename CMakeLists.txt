project(flow_complex)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
include(cotire)

set(BASE_FLAGS "-std=c++11 -fcolor-diagnostics")
set(DIAG_FLAGS "-Wall -Wextra -g")
set(REL_FLAGS  "-O2 -DNDEBUG -pipe -fomit-frame-pointer")
set (CMAKE_CXX_FLAGS ${BASE_FLAGS})
set (CMAKE_CXX_FLAGS_DEBUG "${DIAG_FLAGS} -O0")
set (CMAKE_CXX_FLAGS_RELEASE ${REL_FLAGS})
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO ${REL_FLAGS} ${DIAG_FLAGS})

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

