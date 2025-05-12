function [sql_stations_type] = map_sql_stations_type(graph)
        
num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 1)
    sql_stations_type = [];
    return;
end

%% Station type
% MATLAB 1 is SQL ReMi station w/o backflow (set pressure)
% MATLAB 2 is SQL Injection station w/ pressure control (set flow)
% MATLAB 3 is SQL Outlet station / Consumption point w/o pressure control (set flow)
% MATLAB 0 is SQL Junction (default)

sql_stations_type = zeros(1, num_nodes);
for s = 1:num_nodes
    if graph.Nodes.Type(s) == 1 % ReMi station w/o backflow (set pressure)
        sql_stations_type(s) = 1;
    elseif graph.Nodes.Type(s) == 2 % Injection station w/ pressure control (set flow)
        sql_stations_type(s) = 2;
    elseif graph.Nodes.Type(s) == 3 % Outlet station / Consumption point w/o pressure control  (set flow)
        sql_stations_type(s) = 3;
    elseif graph.Nodes.Type(s) == 0 % Junction
        sql_stations_type(s) = 4;
    else
        sql_stations_type(s) = 4; % DEFAULT is Junction
    end
end


end