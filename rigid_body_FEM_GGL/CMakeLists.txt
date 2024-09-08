find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

find_package(LAPACK REQUIRED)

pybind11_add_module(rigid_body_FEM_GGL bind_rigid_body_FEM_GGL.cc)
target_link_libraries (rigid_body_FEM_GGL PUBLIC)
target_link_libraries (rigid_body_FEM_GGL PUBLIC LAPACK::LAPACK)

add_executable (test_rb_fem_GGL rbs_fem_GGL.cc)
target_link_libraries (test_rb_fem_GGL PUBLIC)
target_link_libraries (test_rb_fem_GGL PUBLIC LAPACK::LAPACK)

target_compile_options(rigid_body_FEM_clean PRIVATE "-flto=auto")