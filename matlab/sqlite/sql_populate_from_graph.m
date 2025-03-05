function [num_stations_converted] = sql_populate_from_graph(db_path, graph)

%% Stations
num_stations_converted = sql_populate_stations(db_path, graph);
assert(num_stations_converted == size(graph.Nodes.Nodes_ID, 1));


end