set(SHIMMER_SOURCES_T ../src/infrastructure_graph.cpp 
        ../src/geometry_properties.cpp 
        ../src/conservation_matrices.cpp 
        ../src/incidence_matrix.cpp 
        ../src/assemble.cpp
        verify_test.cpp
        ../src/matrix_manipulations.cpp
        ../src/matlab_manip.cpp 
        ../src/pipe_calculator.cpp
        ../src/gas_law.cpp)

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

add_executable(test_conservation_system test_conservation_system.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_conservation_system ${LINK_LIBS})
add_test(NAME test_conservation_system COMMAND test_conservation_system)

add_executable(test_friction_factor test_friction_factor.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_friction_factor ${LINK_LIBS})
add_test(NAME test_friction_factor COMMAND test_friction_factor)

add_executable(test_speed_of_sound test_speed_of_sound.cpp ${SHIMMER_SOURCES_T})
target_link_libraries(test_speed_of_sound ${LINK_LIBS})
add_test(NAME test_speed_of_sound COMMAND test_speed_of_sound)
