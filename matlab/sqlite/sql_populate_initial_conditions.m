function [] = sql_populate_initial_conditions(db_path, graph)

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 2)
    return;
end

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    return;
end

%% Initial conditions
%
% graph.Nodes.PRESSURES (p)
% graph.Nodes.G_EXE (L)
% graph.Edges.FLOWRATES (G)
% Ask if it is the guess or the last value

stations_initial_condition = cell(num_nodes, 3);
for s = 1:num_nodes
    stations_initial_condition{s, 1} = graph.Nodes.Nodes_ID(s);
    stations_initial_condition{s, 2} = graph.Nodes.PRESSURES(s);
    stations_initial_condition{s, 3} = graph.Nodes.G_EXE(s);
end

stations_initial_condition_tab_variables_name = ["s_number", "init_P", "init_L"];
stations_initial_condition_tab = cell2table(stations_initial_condition, ...
    "VariableNames", stations_initial_condition_tab_variables_name );

conn = sqlite(db_path, 'connect');
sqlwrite(conn, "station_initial_conditions", stations_initial_condition_tab);
close(conn);

end