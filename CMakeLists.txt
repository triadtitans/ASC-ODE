cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)


find_package(NGSolve CONFIG REQUIRED
  HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)


include_directories(src)
# include_directories(/Applications/Netgen.app/Contents/Resources/include)
# link_directories(/Applications/Netgen.app/Contents/MacOS)

add_subdirectory (src)




add_executable (test_ode demos/test_ode.cc)
target_sources (test_ode PUBLIC src/ode.h src/nonlinfunc.h)
target_link_libraries (test_ode PUBLIC ngbla)


add_executable (test_verlet demos/test_verlet.cc)
target_sources (test_verlet PUBLIC src/ode.h src/nonlinfunc.h)
target_link_libraries (test_verlet PUBLIC ngbla)

add_executable (test_newmark demos/test_newmark.cc)
target_sources (test_newmark PUBLIC src/ode.h src/nonlinfunc.h)
target_link_libraries (test_newmark PUBLIC ngbla)

add_executable (test_shake demos/test_shake.cc)
target_sources (test_shake PUBLIC src/ode.h src/nonlinfunc.h)
target_link_libraries (test_shake PUBLIC ngbla)

add_executable (test_alpha demos/test_alpha.cc)
target_sources (test_alpha PUBLIC src/ode.h src/nonlinfunc.h)
target_link_libraries (test_alpha PUBLIC ngbla)


add_subdirectory (mass_spring)
