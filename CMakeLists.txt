cmake_minimum_required(VERSION 3.15...3.27)
project(ASC-ode)

set (CMAKE_CXX_STANDARD 17)

include_directories(src)
include_directories(BLA/src)
include_directories(BLA/HPC/src)
include_directories(BLA/HPC/concurrentqueue)

# if invoked using scikit-build-core
if(DEFINED SKBUILD)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)

    # rigid_body
    python_add_library(rigid_body MODULE rigid_body/bind_rigid_body.cc WITH_SOABI)
    target_link_libraries(rigid_body PRIVATE pybind11::headers)

    target_compile_definitions(rigid_body PRIVATE VERSION_INFO=${PROJECT_VERSION})

    #rigid_body_fem
    python_add_library(rigid_body_FEM MODULE rigid_body_FEM/bind_rigid_body_FEM.cc WITH_SOABI)
    target_link_libraries(rigid_body_FEM PRIVATE pybind11::headers)

    target_compile_definitions(rigid_body_FEM PRIVATE VERSION_INFO=${PROJECT_VERSION})


    install(TARGETS rigid_body rigid_body_FEM DESTINATION lib_rigid_body)

else()
    set (CXXFLAGS "-g -O0")

    add_subdirectory (src)

    add_executable (test_ode demos/test_ode.cc)
    target_link_libraries (test_ode PUBLIC)

    add_executable (test_newmark demos/test_newmark.cc)
    target_link_libraries (test_newmark PUBLIC)

    add_executable (test_alpha demos/test_alpha.cc)
    target_link_libraries (test_alpha PUBLIC)

    add_executable (test_autodiff demos/test_autodiff.cc)
    target_link_libraries (test_autodiff PUBLIC)

    add_executable (test_cube demos/test_cube.cc)
    target_link_libraries (test_cube PUBLIC)

    add_subdirectory (mass_spring)
    add_subdirectory (rigid_body)
    add_subdirectory (rigid_body_FEM)
endif()