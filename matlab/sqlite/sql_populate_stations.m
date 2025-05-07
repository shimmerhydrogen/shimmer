function [num_nodes, num_nodes_converted] = sql_populate_stations(db_path, graph)

num_nodes = max(size(graph.Nodes.Nodes_ID, 1), size(graph.Nodes.Nodes_ID, 2));

if (num_nodes == 0)
    num_nodes_converted = 0;
    return;
end

num_variables = 3;
stations_tab_variables_name = ["s_number", "s_name", "t_type"];
altitude_index = 0;
latitude_index = 0;
longitude_index = 0;

if isfield(graph.Nodes, 'altitude') && size(graph.Nodes.altitude, 1) == num_nodes
    num_variables = num_variables + 1;
    altitude_index = num_variables;
    stations_tab_variables_name(altitude_index) = 's_height';
end

if isfield(graph.Nodes, 'coordinates_XY') && size(graph.Nodes.coordinates_XY, 1) == num_nodes
    num_variables = num_variables + 2;
    latitude_index = num_variables - 1;
    longitude_index = num_variables;
    stations_tab_variables_name(latitude_index) = 's_latitude';
    stations_tab_variables_name(longitude_index) = 's_longitude';
end

%% Check Station type
% MATLAB 1 is SQL ReMi station w/o backflow (set pressure)
% MATLAB 2 is SQL Injection station w/ pressure control (set flow)
% MATLAB 3 is SQL Outlet station / Consumption point w/o pressure control
% MATLAB 0 is SQL Junction

%% Boundary conditions
% MATLAB 1 trovi il valore in graph.Nodes.PRESSURES
% MATLAB 2/3 trovi il valore in graph.Nodes.G_EXE (G is flow exchanged)
% Metti questo nelle tabelle profiles con time 0
% Invece per i limits non so dove andarli a prendere, ma mi sa che sarà un
% valore settato a manina (da chiedere) -> vedere file small_network_disma_nodes_nonpipe.xlsx

%% Initial conditions
%
% graph.Nodes.PRESSURES (p)
% graph.Nodes.G_EXE (L)
% graph.Nodes.FLOWRATES (G)
% Ask if it is the guess or the last value


stations_type = zeros(num_nodes);
for s = 1:num_nodes
    if graph.Nodes.Type(s) == 1 % ReMi station w/o backflow (set pressure)
        stations_type(s) = 0;
    elseif graph.Nodes.Type(s) == 2 % Injection station w/ pressure control (set flow)
        stations_type(s) = 1;
    elseif graph.Nodes.Type(s) == 3 % Outlet station / Consumption point w/o pressure control
        stations_type(s) = 4; % Outlet station not used
    elseif graph.Nodes.Type(s) == 0 % Junction
        stations_type(s) = 3;
    else
        stations_type(s) = 3; % DEFAULT is Junction
    end
end


stations = cell(num_nodes, num_variables);
for s = 1:num_nodes
    stations{s, 1} = graph.Nodes.Nodes_ID(s);
    stations{s, 2} =  ['Station_', num2str(s)];
    stations{s, 3} = stations_type(s);
end

if altitude_index > 0
    for s = 1:num_nodes
        stations{s, altitude_index} = graph.Nodes.altitude(s);
    end
end

if latitude_index > 0 && longitude_index > 0
    for s = 1:num_nodes
        stations{s, latitude_index} = graph.Nodes.coordinates_XY(s, 1);
        stations{s, longitude_index} = graph.Nodes.coordinates_XY(s, 2);
    end
end

stations_tab = cell2table(stations, ...
    "VariableNames", stations_tab_variables_name );

conn = sqlite(db_path, 'connect');
sqlwrite(conn, "stations", stations_tab);
close(conn);

num_nodes_converted = size(stations_tab, 1);

end