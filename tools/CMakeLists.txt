find_package(gflags REQUIRED)

add_executable(fc fc.cpp)
target_link_libraries(fc ${TBB_LIBRARIES} glog::glog gflags)

add_executable(perturb perturb.cpp)
target_link_libraries(perturb glog::glog)

add_executable(simplify simplify.cpp)
target_link_libraries(simplify ${TBB_LIBRARIES} glog::glog)

add_executable(fc2pc fc2pc.cpp)
target_link_libraries(fc2pc ${TBB_LIBRARIES} glog::glog)

add_executable(inc-pts inc_pts.cpp)
target_link_libraries(inc-pts ${TBB_LIBRARIES} glog::glog)

if (GOOGLE_PERFTOOLS_FOUND)
  target_link_libraries(fc ${TCMALLOC_LIBRARIES})
endif()

