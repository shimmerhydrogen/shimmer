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

if isfield(graph.Nodes, 'altitude')
    num_variables = num_variables + 1;
    altitude_index = num_variables;
    stations_tab_variables_name(altitude_index) = 's_height';
end

if isfield(graph.Nodes, 'coordinates_XY')
    num_variables = num_variables + 2;
    latitude_index = num_variables - 1;
    longitude_index = num_variables;
    stations_tab_variables_name(latitude_index) = 's_latitude';
    stations_tab_variables_name(longitude_index) = 's_longitude';
end

stations = cell(num_nodes, num_variables);
for s = 1:num_nodes
    stations{s, 1} = graph.Nodes.Nodes_ID(s);
    stations{s, 2} =  ['Station_', num2str(s)];
    stations{s, 3} = graph.Nodes.Type(s);
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