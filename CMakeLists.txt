cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)


# use linalg from NGSolve, otherwise comment it out
#find_package(NGSolve CONFIG REQUIRED
#  HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
#)

set (CXXFLAGS "-g -O0")
include_directories(src)
add_subdirectory (src)

include_directories(BLA/src)

include_directories(BLA/HPC/src)

add_executable (test_ode demos/test_ode.cc)
target_link_libraries (test_ode PUBLIC)

add_executable (test_newmark demos/test_newmark.cc)
target_link_libraries (test_newmark PUBLIC)

add_executable (test_alpha demos/test_alpha.cc)
target_link_libraries (test_alpha PUBLIC)

add_executable (test_cube demos/test_cube.cc)
target_link_libraries (test_cube PUBLIC)

add_subdirectory (mass_spring)
