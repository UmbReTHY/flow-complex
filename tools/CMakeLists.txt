add_executable(fc fc.cpp)
target_link_libraries(fc ${TBB_LIBRARIES})

if (GOOGLE_PERFTOOLS_FOUND)
  target_link_libraries(fc ${TCMALLOC_LIBRARIES})
endif()

