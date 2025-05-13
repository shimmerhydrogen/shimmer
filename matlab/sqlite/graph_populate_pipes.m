function [graph] = graph_populate_pipes(db_path, graph)

db = sqlite(db_path, 'readonly');
pipe_tab = sqlread(db, "pipelines");
close(db);

num_pipes = size(pipe_tab, 1);

if (num_pipes < 2)
    return;
end

pipes_type = map_matlab_pipes_type(pipe_tab);

graph.Edges.EndNodes = zeros(num_pipes, 2);
graph.Edges.PIPE = zeros(num_pipes, 1);
graph.Edges.COMP = zeros(num_pipes, 1);
graph.Edges.VALV = zeros(num_pipes, 1);
graph.Edges.REDST = zeros(num_pipes, 1);
graph.Edges.REG = zeros(num_pipes, 1);

for p = 1:num_pipes
    graph.Edges.EndNodes(p, 1) = pipe_tab.s_from(p);
    graph.Edges.EndNodes(p, 2) = pipe_tab.s_to(p);

    if pipes_type(p) == 0 % pipe
        graph.Edges.PIPE(p) = 1;
    elseif pipes_type(p) == 1 % compressor
        graph.Edges.COMP(p) = 1;
    elseif pipes_type(p) == 3 % valve
        graph.Edges.VALV(p) = 1;
    elseif pipes_type(p) == 2 % reduction and regulation station
        graph.Edges.REDST(p) = 1;
        graph.Edges.REG(p) = 1;
    else
        graph.Edges.PIPE(p) = 1; % DEFAULT is pipe
    end
end

db = sqlite(db_path, 'readonly');
pipe_parameters_tab = sqlread(db, "pipe_parameters");
close(db);

if (size(pipe_parameters_tab, 1) < 1)
    return;
end

graph.Edges.Length = zeros(num_pipes, 1);
graph.Edges.Diameter = zeros(num_pipes, 1);
graph.Edges.Epsi = zeros(num_pipes, 1);

for p = 1:size(pipe_parameters_tab, 1)
    pipe_extremes = [pipe_parameters_tab.s_from(p), pipe_parameters_tab.s_to(p)];
    pipe_index = find(all(graph.Edges.EndNodes == pipe_extremes, 2));
    graph.Edges.Length(pipe_index) = pipe_parameters_tab.length(p);
    graph.Edges.Diameter(pipe_index) = pipe_parameters_tab.diameter(p);
    graph.Edges.Epsi(pipe_index) = pipe_parameters_tab.roughness(p);
end

end