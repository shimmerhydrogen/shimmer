cmake_minimum_required(VERSION 3.22)
project(libgerg)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matlab REQUIRED)
list(APPEND LINK_LIBS Matlab::MatlabEngine)
list(APPEND LINK_LIBS Matlab::MatlabDataArray)

set(GERG_SOURCES
    src/MATLAB_GERG_functions.cpp
    src/Matlab_interface.cpp)

add_library(gerg SHARED ${GERG_SOURCES})
add_library(GERG::GERG ALIAS gerg)
target_link_libraries(gerg ${LINK_LIBS})
target_include_directories(gerg PUBLIC include)
