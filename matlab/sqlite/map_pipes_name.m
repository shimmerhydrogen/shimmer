function [pipes_name] = map_pipes_name(graph)

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    pipes_name = {};
    return;
end

pipes_name = cell(1, num_pipes);

for p = 1:num_pipes
    pipes_name{p} = [num2str(p)];
end

end