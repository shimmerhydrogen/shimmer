function [] = sql_populate_initial_conditions(db_path, graph)

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 2)
    return;
end

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    return;
end

%% Nodes Initial conditions
%
% graph.Nodes.PRESSURES (p)
% graph.Nodes.G_EXE (L)

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

%% Pipes Initial conditions
%
% graph.Edges.FLOWRATES (G)

pipelines_initial_condition = cell(num_nodes, 4);
for p = 1:num_pipes
    pipelines_initial_condition{p, 1} = [num2str(p)];
    pipelines_initial_condition{p, 2} = graph.Edges.EndNodes(p, 1);
    pipelines_initial_condition{p, 3} = graph.Edges.EndNodes(p, 2);
    pipelines_initial_condition{p, 4} = graph.Edges.FLOWRATES(p);
end

pipelines_initial_condition_tab_variables_name = ["p_name", "s_from", "s_to", "init_G"];
pipelines_initial_condition_tab = cell2table(pipelines_initial_condition, ...
    "VariableNames", pipelines_initial_condition_tab_variables_name );

conn = sqlite(db_path, 'connect');
sqlwrite(conn, "pipe_initial_conditions", pipelines_initial_condition_tab);
close(conn);

end