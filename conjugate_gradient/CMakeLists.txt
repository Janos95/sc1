cmake_minimum_required(VERSION 3.10)

project(FiniteDifferences LANGUAGES CXX)

add_executable(FiniteDifferences src/main.cpp)

set_target_properties(FiniteDifferences PROPERTIES
    CXX_STANDARD 17
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -lpthread")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Eigen3)
find_package(MKL)
find_package(PCL 1.2 REQUIRED)

find_package(Boost REQUIRED COMPONENTS
    timer 
)

target_include_directories(FiniteDifferences PUBLIC 
    Eigen3::eigen
    ${MKL_INCLUDE_DIRS}
    ${PCL_INCLUDE_DIRS}
)

set(MKL_LINKING -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group)


target_link_libraries(FiniteDifferences Boost::timer ${PCL_LIBRARIES} ${MKL_LINKING} gomp pthread m dl)

