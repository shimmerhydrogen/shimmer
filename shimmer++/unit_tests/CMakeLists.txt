
# verify_test should be called test_common or something like this
add_library(verify_test STATIC verify_test.cpp)
target_link_libraries(verify_test shimmer) #for eigen include paths

add_executable(test_incidence test_incidence_matrix.cpp)
target_link_libraries(test_incidence shimmer verify_test)
add_test(NAME test_incidence COMMAND test_incidence)

add_executable(test_geometry test_geometry.cpp)
target_link_libraries(test_geometry shimmer verify_test)
add_test(NAME test_geometry COMMAND test_geometry)

add_executable(test_conservation test_conservation_matrices.cpp)
target_link_libraries(test_conservation shimmer verify_test)
add_test(NAME test_conservation COMMAND test_conservation)

add_executable(test_resistances test_resistances.cpp)
target_link_libraries(test_resistances shimmer verify_test)
add_test(NAME test_resistances COMMAND test_resistances)

add_executable(test_mean_pressure test_mean_pressure.cpp)
target_link_libraries(test_mean_pressure shimmer verify_test)
add_test(NAME test_mean_pressure COMMAND test_mean_pressure)

add_executable(test_conservation_system test_conservation_system.cpp)
target_link_libraries(test_conservation_system shimmer verify_test)
add_test(NAME test_conservation_system COMMAND test_conservation_system)

add_executable(test_friction_factor test_friction_factor.cpp)
target_link_libraries(test_friction_factor shimmer verify_test)
add_test(NAME test_friction_factor COMMAND test_friction_factor)

#add_executable(test_speed_of_sound test_speed_of_sound.cpp ${SHIMMER_SOURCES_T})
#target_link_libraries(test_speed_of_sound ${LINK_LIBS})
#add_test(NAME test_speed_of_sound COMMAND test_speed_of_sound)

add_executable(test_fluid_solver test_fluid_solver.cpp)
target_link_libraries(test_fluid_solver shimmer verify_test)
add_test(NAME test_fluid_solver COMMAND test_fluid_solver)

add_executable(test_time_cycle test_time_cycle.cpp)
target_link_libraries(test_time_cycle shimmer verify_test)
add_test(NAME test_time_cycle COMMAND test_time_cycle)

add_executable(test_stations test_stations.cpp)
target_link_libraries(test_stations shimmer verify_test)
add_test(NAME test_stations COMMAND test_stations)

add_executable(test_non-pipe_at_nodes test_non-pipe_at_nodes.cpp)
target_link_libraries(test_non-pipe_at_nodes shimmer verify_test)
add_test(NAME test_non-pipe_at_nodes COMMAND test_stations)

add_executable(test_sqlite_import test_sqlite_import.cpp)
target_link_libraries(test_sqlite_import shimmer)
add_test(NAME test_sqlite_import COMMAND test_sqlite_import)

add_executable(test_compressor test_compressor.cpp)
target_link_libraries(test_compressor shimmer verify_test)
add_test(NAME test_compressor COMMAND test_compressor)




add_executable(test_non-pipe_at_nodes_withdb test_non-pipe_at_nodes_withdb.cpp)
target_link_libraries(test_non-pipe_at_nodes_withdb shimmer verify_test)
add_test(NAME test_non-pipe_at_nodes_withdb COMMAND test_non-pipe_at_nodes_withdb)
