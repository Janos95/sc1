cmake_minimum_required(VERSION 3.10)

project(FiniteDifferences LANGUAGES CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

set(CMAKE_CXX_COMPILER /usr/bin/g++-8)

add_executable(FiniteDifferences src/main.cpp)

set_target_properties(FiniteDifferences PROPERTIES
    CXX_STANDARD 17
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -lpthread")

find_package(PCL 1.2 REQUIRED)

find_package(Boost REQUIRED COMPONENTS
    timer 
)

# find mkl!
set(MKLROOT /home/janos/intel/compilers_and_libraries_2019.0.117/linux/mkl)
set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
set(MKL_LINK_LIBRARIES  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group tbb stdc++ pthread m dl)

target_include_directories(FiniteDifferences PUBLIC 
#     Eigen3::eigen
    ${MKL_INCLUDE_DIRS}
#     ${PCL_INCLUDE_DIRS}
)

target_link_libraries(FiniteDifferences  ${MKL_LINK_LIBRARIES})


target_link_libraries(FiniteDifferences 
    CONAN_PKG::fmt 
    CONAN_PKG::eigen
    Boost::timer)
    

# ${PCL_LIBRARIES}
