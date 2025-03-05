function [num_stations_converted] = sql_populate_from_graph(db_path, graph)

%% Stations
[num_nodes, num_stations_converted] = sql_populate_stations(db_path, graph);
assert(num_nodes == num_stations_converted);


end