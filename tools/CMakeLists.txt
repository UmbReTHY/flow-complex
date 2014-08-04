add_executable(fc fc.cpp)
target_link_libraries(fc ${TBB_LIBRARIES})

add_executable(perturb perturb.cpp)

if (GOOGLE_PERFTOOLS_FOUND)
  target_link_libraries(fc ${TCMALLOC_LIBRARIES})
endif()

