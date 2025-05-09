function [num_stations_converted,numpipes_converted] = sql_populate_from_graph(db_path, graph)

%% Stations
[num_nodes, num_stations_converted] = sql_populate_stations(db_path, graph);
assert(num_nodes == num_stations_converted);

%% Pipes
[num_pipes, numpipes_converted] = sql_populate_pipes(db_path, graph);
assert(num_pipes == numpipes_converted);

%% Boundary conditions
sql_populate_boundary_conditions(db_path, graph);

%% Initial Conditions
sql_populate_initial_conditions(db_path, graph);

end