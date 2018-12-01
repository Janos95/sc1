cmake_minimum_required(VERSION 3.10)

project(FiniteDifferences LANGUAGES CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

add_executable(FiniteDifferences src/main.cpp)

set_target_properties(FiniteDifferences PROPERTIES
    CXX_STANDARD 17
)

find_package(VTK REQUIRED)

# find mkl and link it to project
set(MKLROOT /home/janos/intel/compilers_and_libraries_2019.0.117/linux/mkl)
set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
set(MKL_LINK_LIBRARIES  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group tbb stdc++ pthread m dl)

target_include_directories(FiniteDifferences PUBLIC
        ${MKL_INCLUDE_DIRS}
        ${VTK_INCLUDE_DIRS}

)
target_link_libraries(FiniteDifferences
        ${MKL_LINK_LIBRARIES}
        ${VTK_LIBRARIES})


add_subdirectory(plot)

target_link_libraries(FiniteDifferences
        CONAN_PKG::fmt
        CONAN_PKG::eigen
        CONAN_PKG::Boost::timer)


