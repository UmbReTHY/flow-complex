find_package(GMPxx REQUIRED)
find_package(Boost 1.53.0)

include_directories(SYSTEM
                          ${GMP_INCLUDES}
                          ${Boost_INCLUDE_DIRS})

add_executable(simple simple.cpp)
target_link_libraries(simple ${GMP_LIBRARIES})

