# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /home/jh/code/JHPC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jh/code/JHPC/build

# Include any dependencies generated for this target.
include CMakeFiles/dpo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/dpo.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/dpo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dpo.dir/flags.make

CMakeFiles/dpo.dir/src/bops.cpp.o: CMakeFiles/dpo.dir/flags.make
CMakeFiles/dpo.dir/src/bops.cpp.o: ../src/bops.cpp
CMakeFiles/dpo.dir/src/bops.cpp.o: CMakeFiles/dpo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jh/code/JHPC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dpo.dir/src/bops.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dpo.dir/src/bops.cpp.o -MF CMakeFiles/dpo.dir/src/bops.cpp.o.d -o CMakeFiles/dpo.dir/src/bops.cpp.o -c /home/jh/code/JHPC/src/bops.cpp

CMakeFiles/dpo.dir/src/bops.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dpo.dir/src/bops.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jh/code/JHPC/src/bops.cpp > CMakeFiles/dpo.dir/src/bops.cpp.i

CMakeFiles/dpo.dir/src/bops.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dpo.dir/src/bops.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jh/code/JHPC/src/bops.cpp -o CMakeFiles/dpo.dir/src/bops.cpp.s

# Object files for target dpo
dpo_OBJECTS = \
"CMakeFiles/dpo.dir/src/bops.cpp.o"

# External object files for target dpo
dpo_EXTERNAL_OBJECTS =

../bin/dpo: CMakeFiles/dpo.dir/src/bops.cpp.o
../bin/dpo: CMakeFiles/dpo.dir/build.make
../bin/dpo: CMakeFiles/dpo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jh/code/JHPC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/dpo"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dpo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dpo.dir/build: ../bin/dpo
.PHONY : CMakeFiles/dpo.dir/build

CMakeFiles/dpo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dpo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dpo.dir/clean

CMakeFiles/dpo.dir/depend:
	cd /home/jh/code/JHPC/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jh/code/JHPC /home/jh/code/JHPC /home/jh/code/JHPC/build /home/jh/code/JHPC/build /home/jh/code/JHPC/build/CMakeFiles/dpo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dpo.dir/depend
