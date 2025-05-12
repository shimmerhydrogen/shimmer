function [stations_type] = map_stations_type(graph)
        
num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 1)
    stations_type = [];
    return;
end

%% Station type
% MATLAB 1 is SQL ReMi station w/o backflow (set pressure)
% MATLAB 2 is SQL Injection station w/ pressure control (set flow)
% MATLAB 3 is SQL Outlet station / Consumption point w/o pressure control (set flow)
% MATLAB 0 is SQL Junction (default)

stations_type = zeros(1, num_nodes);
for s = 1:num_nodes
    if graph.Nodes.Type(s) == 1 % ReMi station w/o backflow (set pressure)
        stations_type(s) = 1;
    elseif graph.Nodes.Type(s) == 2 % Injection station w/ pressure control (set flow)
        stations_type(s) = 2;
    elseif graph.Nodes.Type(s) == 3 % Outlet station / Consumption point w/o pressure control  (set flow)
        stations_type(s) = 3;
    elseif graph.Nodes.Type(s) == 0 % Junction
        stations_type(s) = 4;
    else
        stations_type(s) = 4; % DEFAULT is Junction
    end
end


end