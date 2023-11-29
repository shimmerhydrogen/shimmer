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
set(LINK_LIBS ${LINK_LIBS} Qt6::Core Qt6::Widgets)


#######################################################################
# OpenXLSX stuff
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "")
set(OPENXLSX_BUILD_BENCHMARKS OFF CACHE BOOL "")
set(OPENXLSX_BUILD_SAMPLES OFF CACHE BOOL "")
set(OPENXLSX_BUILD_TESTS OFF CACHE BOOL "")
set(OPENXLSX_CREATE_DOCS OFF CACHE BOOL "")
add_subdirectory(deps/OpenXLSX)
set(LINK_LIBS ${LINK_LIBS} OpenXLSX::OpenXLSX)

#######################################################################
# Boost stuff
find_package(Boost REQUIRED COMPONENTS graph)

#######################################################################
# Lua stuff
find_package(Lua REQUIRED)
include_directories("${LUA_INCLUDE_DIR}")
set(LINK_LIBS ${LINK_LIBS} ${LUA_LIBRARIES})
set(BUILD_TESTING OFF CACHE BOOL "")
set(SOL2_BUILD_LUA OFF CACHE BOOL "")
add_subdirectory(deps/sol2)
set(LINK_LIBS ${LINK_LIBS} sol2::sol2)
######################################################################
## Eigen3
#add_subdirectory(contrib/eigen)
find_package(Eigen3 REQUIRED)
set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

######################################################################

add_executable(shimmer src/main.cpp src/xlsx_io.cpp src/infrastructure_graph.cpp)
target_link_libraries(shimmer ${LINK_LIBS})

add_executable(test_incidence src/test_incidence_matrix.cpp)
target_link_libraries(test_incidence ${LINK_LIBS})
