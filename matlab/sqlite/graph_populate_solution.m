function [graph] = graph_populate_solution(db_path, graph)

%% FOR REAL SOLUTION
% to count total time intervals
% SELECT s_number, count(DISTINCT time) AS unique_count FROM 'table_solution' GROUP BY s_number LIMIT 0,30

num_time_intervals = 1;

%% Nodes Initial conditions
%
% graph.Nodes.PRESSURES (p)
% graph.Nodes.G_EXE (L)

db = sqlite(db_path, 'readonly');
station_solution_tab = sqlread(db, "station_initial_conditions");
close(db);

num_nodes = size(station_solution_tab, 1);

if (num_nodes < 1)
    return;
end

graph.Nodes.PRESSURES = zeros(num_nodes, num_time_intervals);
graph.Nodes.G_EXE = zeros(num_nodes, num_time_intervals);

for s = 1:num_nodes
    s_index = station_solution_tab.s_number(s);
    graph.Nodes.PRESSURES(s_index) = station_solution_tab.init_P(s);
    graph.Nodes.G_EXE(s_index) = station_solution_tab.init_L(s);
end

%% Pipes Initial conditions
%
% graph.Edges.FLOWRATES (G)

db = sqlite(db_path, 'readonly');
pipe_solution_tab = sqlread(db, "pipe_initial_conditions");
close(db);

num_pipes = size(pipe_solution_tab, 1);

if (num_pipes < 1)
    return;
end

graph.Edges.FLOWRATES = zeros(num_pipes, num_time_intervals);

for p = 1:num_pipes
    pipe_extremes = [pipe_solution_tab.s_from(p), pipe_solution_tab.s_to(p)];
    pipe_index = find(all(graph.Edges.EndNodes == pipe_extremes, 2));
    graph.Edges.FLOWRATES(pipe_index) = pipe_solution_tab.init_G(p);
end

end