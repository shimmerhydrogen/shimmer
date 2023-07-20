cmake_minimum_required(VERSION 3.18)
project(shimmer++)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#######################################################################
# Qt stuff
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt6 COMPONENTS Core Widgets REQUIRED)

#######################################################################
# OpenXLSX stuff
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "")
set(OPENXLSX_BUILD_BENCHMARKS OFF CACHE BOOL "")
set(OPENXLSX_BUILD_SAMPLES OFF CACHE BOOL "")
set(OPENXLSX_BUILD_TESTS OFF CACHE BOOL "")
set(OPENXLSX_CREATE_DOCS OFF CACHE BOOL "")
add_subdirectory(deps/OpenXLSX)

find_package(Boost REQUIRED COMPONENTS graph)

add_executable(shimmer src/main.cpp)
target_link_libraries(shimmer OpenXLSX::OpenXLSX)
target_link_libraries(shimmer Qt6::Widgets)

