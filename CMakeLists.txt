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

    #rigid_body_fem
    python_add_library(rigid_body_FEM MODULE rigid_body_FEM/bind_rigid_body_FEM.cc WITH_SOABI)
    target_link_libraries(rigid_body_FEM PRIVATE pybind11::headers)

    target_compile_definitions(rigid_body_FEM PRIVATE VERSION_INFO=${PROJECT_VERSION})


    install(TARGETS rigid_body_FEM DESTINATION lib_rigid_body)

else()
    set (CXXFLAGS "-g -O0")

    add_subdirectory (src)
    
    add_subdirectory (rigid_body_FEM)

    add_subdirectory (rigid_body_FEM_GGL)

    add_subdirectory (rigid_body_FEM_SR)
endif()