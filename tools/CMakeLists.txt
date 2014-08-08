add_executable(fc fc.cpp)
target_link_libraries(fc ${TBB_LIBRARIES})

add_executable(perturb perturb.cpp)

add_executable(simplify simplify.cpp)
target_link_libraries(simplify ${TBB_LIBRARIES})

add_executable(fc2pc fc2pc.cpp)
target_link_libraries(fc2pc ${TBB_LIBRARIES})

if (GOOGLE_PERFTOOLS_FOUND)
  target_link_libraries(fc ${TCMALLOC_LIBRARIES})
endif()

