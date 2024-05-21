
add_executable(test_incidence test_incidence_matrix.cpp verify_test.cpp)
target_link_libraries(test_incidence shimmer::shimmer)
add_test(NAME test_incidence COMMAND test_incidence)

add_executable(test_geometry test_geometry.cpp verify_test.cpp)
target_link_libraries(test_geometry shimmer::shimmer)
add_test(NAME test_geometry COMMAND test_geometry)

add_executable(test_conservation test_conservation_matrices.cpp verify_test.cpp)
target_link_libraries(test_conservation shimmer::shimmer)
add_test(NAME test_conservation COMMAND test_conservation)

add_executable(test_resistances test_resistances.cpp verify_test.cpp)
target_link_libraries(test_resistances shimmer::shimmer)
add_test(NAME test_resistances COMMAND test_resistances)

add_executable(test_mean_pressure test_mean_pressure.cpp verify_test.cpp)
target_link_libraries(test_mean_pressure shimmer::shimmer)
add_test(NAME test_mean_pressure COMMAND test_mean_pressure)

add_executable(test_conservation_system test_conservation_system.cpp verify_test.cpp)
target_link_libraries(test_conservation_system shimmer::shimmer)
add_test(NAME test_conservation_system COMMAND test_conservation_system)

add_executable(test_friction_factor test_friction_factor.cpp verify_test.cpp)
target_link_libraries(test_friction_factor shimmer::shimmer)
add_test(NAME test_friction_factor COMMAND test_friction_factor)

#add_executable(test_speed_of_sound test_speed_of_sound.cpp ${SHIMMER_SOURCES_T})
#target_link_libraries(test_speed_of_sound ${LINK_LIBS})
#add_test(NAME test_speed_of_sound COMMAND test_speed_of_sound)

add_executable(test_fluid_solver test_fluid_solver.cpp verify_test.cpp)
target_link_libraries(test_fluid_solver shimmer::shimmer)
add_test(NAME test_fluid_solver COMMAND test_fluid_solver)

add_executable(test_time_cycle test_time_cycle.cpp verify_test.cpp)
target_link_libraries(test_time_cycle shimmer::shimmer)
add_test(NAME test_time_cycle COMMAND test_time_cycle)

add_executable(test_stations test_stations.cpp verify_test.cpp)
target_link_libraries(test_stations shimmer::shimmer)
add_test(NAME test_stations COMMAND test_stations)

add_executable(test_non-pipe_at_nodes test_non-pipe_at_nodes.cpp verify_test.cpp)
target_link_libraries(test_non-pipe_at_nodes shimmer::shimmer)
add_test(NAME test_non-pipe_at_nodes COMMAND test_stations)