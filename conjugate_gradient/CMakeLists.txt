cmake_minimum_required(VERSION 3.8)

project(FiniteDifferences LANGUAGES CXX)

add_executable(FiniteDifferences src/main.cpp)

set_target_properties(FiniteDifferences PROPERTIES
    CXX_STANDARD 17
)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lboost_system")


find_package (Eigen3 REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS

)

target_include_directories(FiniteDifferences PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
