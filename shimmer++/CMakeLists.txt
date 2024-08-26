cmake_minimum_required(VERSION 3.22)
project(shimmer++)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUGSANITIZED "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASENATIVE "-O3 -DNDEBUG -march=native")

option(BUILD_SHIMMER_UNIT_TESTS "Build SHIMMER++ unitary tests" ON)

#######################################################################
# GERG
add_subdirectory(deps/gerg)
set(LINK_LIBS ${LINK_LIBS} GERG::GERG)

#######################################################################
# SQLite
find_package(SQLite3 REQUIRED)
set(LINK_LIBS ${LINK_LIBS} SQLite::SQLite3)

#######################################################################
# Threads library
find_package(Threads REQUIRED)
set(LINK_LIBS ${LINK_LIBS} Threads::Threads)

#######################################################################
# Qt stuff
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTORCC ON)
#set(CMAKE_AUTOUIC ON)
#find_package(Qt6 COMPONENTS Core Widgets REQUIRED)
#set(LINK_LIBS ${LINK_LIBS} Qt6::Core Qt6::Widgets)

#######################################################################
# Boost stuff
find_package(Boost REQUIRED COMPONENTS graph)

#######################################################################
# Lua stuff
#find_package(Lua REQUIRED)
#include_directories("${LUA_INCLUDE_DIR}")
#set(LINK_LIBS ${LINK_LIBS} ${LUA_LIBRARIES})
#set(BUILD_TESTING OFF CACHE BOOL "")
#set(SOL2_BUILD_LUA OFF CACHE BOOL "")
#add_subdirectory(deps/sol2)
#set(LINK_LIBS ${LINK_LIBS} sol2::sol2)
######################################################################
## Eigen3
#add_subdirectory(contrib/eigen)
find_package(Eigen3 REQUIRED)
set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

######################################################################

set(SHIMMER_SOURCES
    src/infrastructure_graph.cpp src/geometry_properties.cpp src/gas_law.cpp
    src/conservation_matrices.cpp src/incidence_matrix.cpp
    src/matrix_manipulations.cpp src/pipe_calculator.cpp src/matlab_manip.cpp
    src/fluid_solver.cpp  ../src/boundary.cpp ../src/nonpipe_over_edges.cpp
    src/variable.cpp src/sqlite.cpp)

#add_executable(shimmer src/main.cpp ${SHIMMER_SOURCES})
#target_link_libraries(shimmer ${LINK_LIBS})

add_library(shimmer SHARED ${SHIMMER_SOURCES})
add_library(shimmer::shimmer ALIAS shimmer)
target_link_libraries(shimmer ${LINK_LIBS})
target_include_directories(shimmer PUBLIC src)

if (BUILD_SHIMMER_UNIT_TESTS)
    include(CTest)
    add_subdirectory(unit_tests)
endif()
