# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/janos/libraries/cmake-3.13.0-rc2/bin/cmake

# The command to remove a file.
RM = /home/janos/libraries/cmake-3.13.0-rc2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/janos/sc1/conjugate_gradient

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/janos/sc1/conjugate_gradient/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/conjugate_gradient.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/conjugate_gradient.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/conjugate_gradient.dir/flags.make

CMakeFiles/conjugate_gradient.dir/src/main.cpp.o: CMakeFiles/conjugate_gradient.dir/flags.make
CMakeFiles/conjugate_gradient.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/janos/sc1/conjugate_gradient/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/conjugate_gradient.dir/src/main.cpp.o"
	/usr/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/conjugate_gradient.dir/src/main.cpp.o -c /home/janos/sc1/conjugate_gradient/src/main.cpp

CMakeFiles/conjugate_gradient.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/conjugate_gradient.dir/src/main.cpp.i"
	/usr/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/janos/sc1/conjugate_gradient/src/main.cpp > CMakeFiles/conjugate_gradient.dir/src/main.cpp.i

CMakeFiles/conjugate_gradient.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/conjugate_gradient.dir/src/main.cpp.s"
	/usr/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/janos/sc1/conjugate_gradient/src/main.cpp -o CMakeFiles/conjugate_gradient.dir/src/main.cpp.s

# Object files for target conjugate_gradient
conjugate_gradient_OBJECTS = \
"CMakeFiles/conjugate_gradient.dir/src/main.cpp.o"

# External object files for target conjugate_gradient
conjugate_gradient_EXTERNAL_OBJECTS =

conjugate_gradient: CMakeFiles/conjugate_gradient.dir/src/main.cpp.o
conjugate_gradient: CMakeFiles/conjugate_gradient.dir/build.make
conjugate_gradient: CMakeFiles/conjugate_gradient.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/janos/sc1/conjugate_gradient/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable conjugate_gradient"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/conjugate_gradient.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/conjugate_gradient.dir/build: conjugate_gradient

.PHONY : CMakeFiles/conjugate_gradient.dir/build

CMakeFiles/conjugate_gradient.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/conjugate_gradient.dir/cmake_clean.cmake
.PHONY : CMakeFiles/conjugate_gradient.dir/clean

CMakeFiles/conjugate_gradient.dir/depend:
	cd /home/janos/sc1/conjugate_gradient/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/janos/sc1/conjugate_gradient /home/janos/sc1/conjugate_gradient /home/janos/sc1/conjugate_gradient/cmake-build-debug /home/janos/sc1/conjugate_gradient/cmake-build-debug /home/janos/sc1/conjugate_gradient/cmake-build-debug/CMakeFiles/conjugate_gradient.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/conjugate_gradient.dir/depend
