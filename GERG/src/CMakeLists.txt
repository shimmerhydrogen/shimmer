list(APPEND gerg_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MATLAB_GERG_interface.hpp)
list(APPEND gerg_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MATLAB_GERG_functions.hpp)
list(APPEND gerg_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/GERG_functions.hpp)

list(APPEND gerg_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MATLAB_GERG_interface.cpp)
list(APPEND gerg_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MATLAB_GERG_functions.cpp)

list(APPEND gerg_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gerg_src_sources ${gerg_src_sources} PARENT_SCOPE)
set(gerg_src_headers ${gerg_src_headers} PARENT_SCOPE)
set(gerg_src_includes ${gerg_src_includes} PARENT_SCOPE)
