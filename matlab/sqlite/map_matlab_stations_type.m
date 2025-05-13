function [matlab_stations_type] = map_matlab_stations_type(station_tab)
        
num_nodes = size(station_tab, 1);

if (num_nodes < 1)
    return;
end

%% Station type
% MATLAB 1 is SQL ReMi station w/o backflow (set pressure)
% MATLAB 2 is SQL Injection station w/ pressure control (set flow)
% MATLAB 3 is SQL Outlet station / Consumption point w/o pressure control (set flow)
% MATLAB 0 is SQL Junction (default)

matlab_stations_type = zeros(1, num_nodes);
for s = 1:num_nodes
    if station_tab.t_type(s) == 1 % ReMi station w/o backflow (set pressure)
        matlab_stations_type(s) = 1;
    elseif station_tab.t_type(s) == 2 % Injection station w/ pressure control (set flow)
        matlab_stations_type(s) = 2;
    elseif station_tab.t_type(s) == 3 % Outlet station / Consumption point w/o pressure control  (set flow)
        matlab_stations_type(s) = 3;
    elseif station_tab.t_type(s) == 4 % Junction
        matlab_stations_type(s) = 0;
    else
        matlab_stations_type(s) = 0; % DEFAULT is Junction
    end
end


end