function [graph] = graph_populate_solution(db_path, graph)

%% Count time intervals for nodes
db = sqlite(db_path, 'readonly');
query = 'SELECT s_number, count(DISTINCT timestep) AS num_time_intervals FROM solution_station_pressures GROUP BY s_number';
station_time_interval_tab = fetch(db, query);
close(db);

num_nodes = size(station_time_interval_tab, 1);

if (num_nodes < 1)
    return;
end

nodes_time_intervals = zeros(num_nodes, 1);
for s = 1:num_nodes
    s_index = station_time_interval_tab.s_number(s);
    nodes_time_intervals(s_index) = station_time_interval_tab.num_time_intervals(s);
end

num_node_time_intervals = max(nodes_time_intervals);

if (num_node_time_intervals < 1)
    return;
end

%% Nodes solutions
%
% graph.Nodes.PRESSURES (p)
% graph.Nodes.G_EXE (L)

db = sqlite(db_path, 'readonly');
station_solution_tab = sqlread(db, "solution_station_pressures");
close(db);

graph.Nodes.PRESSURES = zeros(num_nodes, num_node_time_intervals);
graph.Nodes.G_EXE = zeros(num_nodes, num_node_time_intervals);

for s = 1:size(station_solution_tab,1)
    s_index = station_solution_tab.s_number(s);
    t_index = station_solution_tab.timestep(s) + 1;
    graph.Nodes.PRESSURES(s_index, t_index) = station_solution_tab.pressure(s);
    %graph.Nodes.G_EXE(s_index, t_index) = station_solution_tab.init_L(s);
end

%% Count time intervals for pipes
db = sqlite(db_path, 'readonly');
query = 'SELECT p_name, s_from, s_to, count(DISTINCT timestep) AS num_time_intervals FROM solution_pipe_flowrates GROUP BY p_name, s_from, s_to';
pipe_time_interval_tab = fetch(db, query);
close(db);

num_pipes = size(pipe_time_interval_tab, 1);

if (num_pipes < 1)
    return;
end

pipes_time_intervals = zeros(num_nodes, 1);
for p = 1:num_pipes
    pipe_extremes = [pipe_time_interval_tab.s_from(p), pipe_time_interval_tab.s_to(p)];
    pipe_index = find(all(graph.Edges.EndNodes == pipe_extremes, 2));
    pipes_time_intervals(pipe_index) = pipe_time_interval_tab.num_time_intervals(p);
end

num_pipe_time_intervals = max(pipes_time_intervals);

if (num_pipe_time_intervals < 1)
    return;
end

%% Pipes solutions
%
% graph.Edges.FLOWRATES (G)

db = sqlite(db_path, 'readonly');
pipe_solution_tab = sqlread(db, "solution_pipe_flowrates");
close(db);

graph.Edges.FLOWRATES = zeros(num_pipes, num_node_time_intervals);

for p = 1:size(pipe_solution_tab, 1)
    pipe_extremes = [pipe_solution_tab.s_from(p), pipe_solution_tab.s_to(p)];
    pipe_index = find(all(graph.Edges.EndNodes == pipe_extremes, 2));
    t_index = pipe_solution_tab.timestep(p) + 1;
    graph.Edges.FLOWRATES(pipe_index, t_index) = pipe_solution_tab.flowrate(p);
end


end