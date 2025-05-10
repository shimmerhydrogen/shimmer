add_subdirectory(GERG)

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/shimmer_gerg_utilities.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/shimmer_gerg_concepts.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/shimmer_gerg_data.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/shimmer_gerg_functions.hpp)

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/shimmer_gerg_utilities.cpp)

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)
