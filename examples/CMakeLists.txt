find_package(GooglePerfTools)

add_executable(fc_single fc_single.cpp)
target_link_libraries(fc_single ${TBB_LIBRARIES})

add_executable(fc_tbb_hl fc_tbb_hl.cpp)
target_link_libraries(fc_tbb_hl ${TBB_LIBRARIES})

add_executable(fc_tbb_ll fc_tbb_ll.cpp)
target_link_libraries(fc_tbb_ll ${TBB_LIBRARIES})

add_executable(fc_openmp fc_openmp.cpp)
target_link_libraries(fc_openmp ${TBB_LIBRARIES})

if (GOOGLE_PERFTOOLS_FOUND)
  target_link_libraries(fc_single ${TCMALLOC_LIBRARIES})
  target_link_libraries(fc_tbb_hl ${TCMALLOC_LIBRARIES})
  target_link_libraries(fc_tbb_ll ${TCMALLOC_LIBRARIES})
  target_link_libraries(fc_openmp ${TCMALLOC_LIBRARIES})
endif()

