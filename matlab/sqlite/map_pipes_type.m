function [pipes_type] = map_pipes_type(graph)

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    pipes_type = [];
    return;
end

pipes_type = zeros(1, num_pipes);

% pipes_type(p) = 1; Not implemented - No, resistor not supported yet

for p = 1:num_pipes
    if isfield(graph.Edges, 'PIPE') && size(graph.Edges.PIPE, 1) == num_pipes && graph.Edges.PIPE(p)
        pipes_type(p) = 0;
    elseif isfield(graph.Edges, 'COMP') && size(graph.Edges.COMP, 1) == num_pipes && graph.Edges.COMP(p) % compressor
        pipes_type(p) = 1;
    elseif isfield(graph.Edges, 'VALV') && size(graph.Edges.VALV, 1) == num_pipes && graph.Edges.VALV(p) % valve
        pipes_type(p) = 3;
    elseif isfield(graph.Edges, 'REDST') && size(graph.Edges.REDST, 1) == num_pipes && graph.Edges.REDST(p) % reduction station
    elseif isfield(graph.Edges, 'REG') && size(graph.Edges.REG, 1) == num_pipes && graph.Edges.REG(p) % regulation station
        pipes_type(p) = 2;
    else
        pipes_type(p) = 0; % DEFAULT is 0
    end
end

end