cmake_minimum_required(VERSION 3.25)
project(shimmer++)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wshadow")
set(CMAKE_CXX_FLAGS_DEBUGSANITIZED "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASENATIVE "-O3 -DNDEBUG -march=native")

option(BUILD_SHIMMER_UNIT_TESTS "Build SHIMMER++ unitary tests" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

######################################################################
## Compiler reporting: disable deprecation warnings
option(OPT_DISABLE_DEPRECATION_WARNINGS "Disable deprecation warnings" OFF)
if(OPT_DISABLE_DEPRECATION_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()
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

set(SHIMMER_DBIO_SOURCES
    src/sqlite/sqlite.cpp src/sqlite/sqlite_entry_p_reg.cpp
    src/sqlite/sqlite_entry_l_reg.cpp src/sqlite/sqlite_exit_l_reg.cpp
    src/sqlite/sqlite_initial_conditions.cpp
    src/sqlite/sqlite_compressor.cpp src/sqlite/sqlite_pipe.cpp
    src/sqlite/sqlite_reduction.cpp src/sqlite/sqlite_outlet.cpp
    src/sqlite/sqlite_gases.cpp src/sqlite/sqlite_infra.cpp
)

add_library(shimmer_dbio SHARED ${SHIMMER_DBIO_SOURCES})
add_library(shimmer::shimmer_dbio ALIAS shimmer_dbio)
target_link_libraries(shimmer_dbio ${LINK_LIBS})
target_include_directories(shimmer_dbio PUBLIC src)

set(SHIMMER_SOURCES
    src/infrastructure_graph.cpp src/solver/geometry_properties.cpp src/solver/gas_law.cpp
    src/solver/conservation_matrices.cpp src/solver/incidence_matrix.cpp
    src/solver/matrix_manipulations.cpp src/solver/pipe_calculator.cpp src/solver/matlab_manip.cpp
    src/solver/fluid_solver.cpp src/solver/nonpipe_over_edges.cpp
    src/solver/variable.cpp src/solver/boundary.cpp
    src/infra/infrastructure.cpp
)

add_library(shimmer SHARED ${SHIMMER_SOURCES})
add_library(shimmer::shimmer ALIAS shimmer)
target_link_libraries(shimmer shimmer_dbio ${LINK_LIBS})
target_include_directories(shimmer PUBLIC src)

add_executable(ptd unit_tests/populate_test_db.cpp)
target_link_libraries(ptd shimmer ${LINK_LIBS})
target_include_directories(ptd PUBLIC src)

add_executable(ldb unit_tests/load_database.cpp)
target_link_libraries(ldb shimmer ${LINK_LIBS})
target_include_directories(ldb PUBLIC src)

add_executable(solver src/solver_main.cpp)
target_link_libraries(solver shimmer ${LINK_LIBS})
target_include_directories(solver PUBLIC src)

if (BUILD_SHIMMER_UNIT_TESTS)
    include(CTest)
    add_subdirectory(unit_tests)
endif()
