find_package(GooglePerfTools)

add_executable(simple simple.cpp)
target_link_libraries(simple ${TBB_LIBRARIES})

if (GOOGLE_PERFTOOLS_FOUND)
  target_link_libraries(simple ${TCMALLOC_LIBRARIES})
endif()

