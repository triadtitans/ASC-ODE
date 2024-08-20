find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

find_package(LAPACK REQUIRED)

pybind11_add_module(rigid_body_FEM bind_rigid_body_FEM.cc)
pybind11_add_module(tigid_body_FEM_clean bind_rigid_body_FEM_clean.cc)
target_link_libraries (rigid_body_FEM PUBLIC)
target_link_libraries (rigid_body_FEM PUBLIC LAPACK::LAPACK)
# add_custom_command(TARGET rigid_body_FEM POST_BUILD 
#   COMMAND "${CMAKE_COMMAND}" -E copy 
#      "$<TARGET_FILE:rigid_body_FEM>"
#      "rigid_body_FEM.so" 
#   COMMENT "Copying to output directory")

# fix ltrans warning: https://stackoverflow.com/a/72222512
target_compile_options(rigid_body_FEM PRIVATE "-flto=auto")

add_executable (test_rb_fem rigid_body_fem.cc)
target_link_libraries (test_rb_fem PUBLIC)
target_link_libraries (test_rb_fem PUBLIC LAPACK::LAPACK)

target_compile_options(rigid_body_FEM PRIVATE "-flto=auto")

add_executable (test_rb_fem_clean rbs_fem_clean.cc)
target_link_libraries (test_rb_fem_clean PUBLIC)
target_link_libraries (test_rb_fem_clean PUBLIC LAPACK::LAPACK)

#target_compile_options(rigid_body_FEM_clean PRIVATE "-flto=auto")