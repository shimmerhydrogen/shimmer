function [graph] = graph_populate_stations(db_path, graph)

db = sqlite(db_path, 'readonly');
station_tab = sqlread(db, "stations");
close(db);

num_nodes = size(station_tab, 1);

if (num_nodes < 2)
    return;
end

stations_type = map_matlab_stations_type(station_tab);

graph.Nodes.Nodes_ID = zeros(num_nodes, 1);
graph.Nodes.Type = zeros(num_nodes, 1);
graph.Nodes.altitude = zeros(num_nodes, 1);
graph.Nodes.coordinates_XY = zeros(num_nodes, 2);

for s = 1:num_nodes
    graph.Nodes.Nodes_ID(s) = station_tab.s_number(s);
    graph.Nodes.Type(s) = stations_type(s);
    graph.Nodes.altitude(s) = station_tab.s_height(s);
    graph.Nodes.coordinates_XY(s, 1) = station_tab.s_latitude(s);
    graph.Nodes.coordinates_XY(s, 2) = station_tab.s_longitude(s);
end

end