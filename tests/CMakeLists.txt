include_directories(${CMAKE_SOURCE_DIR})

# this macro simplifies adding test executables
macro (do_test arg)
  add_executable(${arg} ${arg}.cpp)
  add_test (${arg}_test ${arg})
endmacro (do_test)

do_test(point_cloud)
do_test(dynamic_qr)
do_test(affine_hull)
do_test(nn_along_ray)
