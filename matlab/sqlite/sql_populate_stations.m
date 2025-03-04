function [num_nodes_converted] = sql_populate_stations(db_path, graph)

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes == 0)
    num_nodes_converted = 0;
    return;
end

stations = cell(num_nodes, 3);
for s = 1:num_nodes
    stations{s, 1} = graph.Nodes.Nodes_ID(s);
    stations{s, 2} = "s01_entry";
    stations{s, 3} = 0;
end

stations_tab = cell2table(stations, ...
    "VariableNames", ["s_number", "s_name", "t_type"] );

conn = sqlite(db_path, 'connect');
sqlwrite(conn, "stations", stations_tab);
close(conn);

num_nodes_converted = size(stations_tab, 1);

end