cmake_minimum_required(VERSION 3.10)
project(finite_elemente)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

add_executable(finite_elemente
        main.cpp
        constant_size_matrix.hpp
        assemble_problem.hpp)

set_target_properties(finite_elemente PROPERTIES
        CXX_STANDARD 17
        )

target_link_libraries(finite_elemente PUBLIC
        CONAN_PKG::boost
        CONAN_PKG::fmt
        CONAN_PKG::range-v3
        CONAN_PKG::eigen
        CONAN_PKG::openblas)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(finite_elemente PUBLIC OpenMP::OpenMP_CXX)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -march=native")


find_package(PCL REQUIRED)

target_include_directories(finite_elemente PRIVATE ${PCL_INCLUDE_DIRS})
target_link_libraries(finite_elemente PRIVATE ${PCL_LIBRARIES} )
