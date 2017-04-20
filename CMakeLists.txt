project(flow_complex)
cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

set(CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (CMAKE_CXX_FLAGS "-fcolor-diagnostics")
  endif()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
  set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
  set (CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-parameter -DNDEBUG -O3")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # using Intel C++
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
  set (CMAKE_CXX_FLAGS_RELEASE "-Wno-unused-parameter -DNDEBUG -Ofast -mavx")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "[Mm]SVC")
  # using Visual Studio C++
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
  add_definitions(-DNOMINMAX)
endif()

option(BUILD_TESTS
       "builds the tests and provide a test target for Makefiles"
       OFF)

option(BUILD_EXAMPLES "" OFF)

find_package (glog 0.3.4 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(GooglePerfTools)
cmake_policy(SET CMP0012 NEW)
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

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

