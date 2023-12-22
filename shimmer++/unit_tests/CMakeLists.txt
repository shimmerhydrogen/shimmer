set(SHIMMER_SOURCES_T ../src/infrastructure_graph.cpp 
        ../src/geometry_properties.cpp 
        ../src/conservation_matrices.cpp 
        ../src/incidence_matrix.cpp 
        ../src/assemble.cpp
        verify_test.cpp
        ../src/matrix_manipulations.cpp 
        ../src/pipe_calculator.cpp)

#set(SHIMMER_SOURCES_T ../${SHIMMER_SOURCES} verify_test.cpp)

add_executable(test_incidence test_incidence_matrix.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_incidence ${LINK_LIBS})
add_test(NAME test_incidence COMMAND test_incidence)

add_executable(test_geometry test_geometry.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_geometry ${LINK_LIBS})
add_test(NAME test_geometry COMMAND test_geometry)

add_executable(test_conservation test_conservation_matrices.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_conservation ${LINK_LIBS})
add_test(NAME test_conservation COMMAND test_conservation)

add_executable(test_resistances test_resistances.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_resistances ${LINK_LIBS})
add_test(NAME test_resistances COMMAND test_resistances)

add_executable(test_mean_pressure test_mean_pressure.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_mean_pressure ${LINK_LIBS})
add_test(NAME test_mean_pressure COMMAND test_mean_pressure)

add_executable(test_assembler test_assembler.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_assembler ${LINK_LIBS})
add_test(NAME test_assembler COMMAND test_assembler)

add_executable(test_friction_factor test_friction_factor.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_friction_factor ${LINK_LIBS})
add_test(NAME test_friction_factor COMMAND test_friction_factor)