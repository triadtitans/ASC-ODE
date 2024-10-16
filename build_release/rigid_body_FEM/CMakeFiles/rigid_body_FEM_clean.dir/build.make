# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/linus/ghq/github.com/triadtitans/ASC-ODE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release

# Include any dependencies generated for this target.
include rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/compiler_depend.make

# Include the progress variables for this target.
include rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/progress.make

# Include the compile flags for this target's objects.
include rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/flags.make

rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o: rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/flags.make
rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o: /home/linus/ghq/github.com/triadtitans/ASC-ODE/rigid_body_FEM/bind_rigid_body_FEM_clean.cc
rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o: rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o"
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o -MF CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o.d -o CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o -c /home/linus/ghq/github.com/triadtitans/ASC-ODE/rigid_body_FEM/bind_rigid_body_FEM_clean.cc

rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.i"
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/linus/ghq/github.com/triadtitans/ASC-ODE/rigid_body_FEM/bind_rigid_body_FEM_clean.cc > CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.i

rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.s"
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/linus/ghq/github.com/triadtitans/ASC-ODE/rigid_body_FEM/bind_rigid_body_FEM_clean.cc -o CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.s

# Object files for target rigid_body_FEM_clean
rigid_body_FEM_clean_OBJECTS = \
"CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o"

# External object files for target rigid_body_FEM_clean
rigid_body_FEM_clean_EXTERNAL_OBJECTS =

rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so: rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/bind_rigid_body_FEM_clean.cc.o
rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so: rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/build.make
rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so: /opt/intel/oneapi/mkl/latest/lib/libmkl_intel_lp64.so
rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so: /opt/intel/oneapi/mkl/latest/lib/libmkl_sequential.so
rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so: /opt/intel/oneapi/mkl/latest/lib/libmkl_core.so
rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so: rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so"
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rigid_body_FEM_clean.dir/link.txt --verbose=$(VERBOSE)
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM && /usr/bin/strip /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so

# Rule to build all files generated by this target.
rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/build: rigid_body_FEM/rigid_body_FEM_clean.cpython-312-x86_64-linux-gnu.so
.PHONY : rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/build

rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/clean:
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM && $(CMAKE_COMMAND) -P CMakeFiles/rigid_body_FEM_clean.dir/cmake_clean.cmake
.PHONY : rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/clean

rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/depend:
	cd /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/linus/ghq/github.com/triadtitans/ASC-ODE /home/linus/ghq/github.com/triadtitans/ASC-ODE/rigid_body_FEM /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM /home/linus/ghq/github.com/triadtitans/ASC-ODE/build_release/rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : rigid_body_FEM/CMakeFiles/rigid_body_FEM_clean.dir/depend

