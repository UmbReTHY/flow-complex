project(flow_complex)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (CMAKE_CXX_FLAGS "-fcolor-diagnostics")
  endif()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
  set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
  set (CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-parameter -DNDEBUG -O3")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_definitions(-DLOGGING)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
  set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
  set (CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-parameter -DNDEBUG -Ofast -mavx")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_definitions(-DLOGGING)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

option(BUILD_TESTS
       "builds the tests and provide a test target for Makefiles"
       OFF)

option(BUILD_EXAMPLES "" OFF)

find_package (glog 0.3.4 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(GooglePerfTools)

include_directories(SYSTEM
                    ${EIGEN3_INCLUDE_DIR}
                    "${CMAKE_SOURCE_DIR}/nanoflann"
                    ${TBB_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR})

if (BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif (BUILD_TESTS)

if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif (BUILD_EXAMPLES)

add_subdirectory(tools)

