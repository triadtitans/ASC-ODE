add_executable (test_rigid_body rigid_body.cc)
target_link_libraries (test_rigid_body PUBLIC)

add_executable (test_rb_system rb_system.cc)
target_link_libraries (test_rb_system PUBLIC)

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(rigid_body bind_rigid_body.cc)
target_link_libraries (rigid_body PUBLIC)
add_custom_command(TARGET rigid_body POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
     "$<TARGET_FILE:rigid_body>"
     "rigid_body.so" 
  COMMENT "Copying to output directory")