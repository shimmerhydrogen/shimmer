function [num_pipes, num_pipes_converted] = sql_populate_pipes(db_path, graph)

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    num_pipes_converted = 0;
    return;
end

pipes_tab_num_variables = 4;
pipes_tab_variables_name = ["p_name", "s_from", "s_to", "p_type"];

pipes_name = map_pipes_name(graph);
pipes_type = map_pipes_type(graph);

pipes = cell(num_pipes, pipes_tab_num_variables);
for p = 1:num_pipes
    pipes{p, 1} = pipes_name{p};
    pipes{p, 2} = graph.Edges.EndNodes(p, 1);
    pipes{p, 3} = graph.Edges.EndNodes(p, 2);
    pipes{p, 4} = pipes_type(p);
end

pipes_tab = cell2table(pipes, ...
    "VariableNames", pipes_tab_variables_name );

conn = sqlite(db_path, 'connect');
sqlwrite(conn, "pipelines", pipes_tab);
close(conn);


pipes_parameters_tab_num_variables = 3;
pipes_parameters_tab_variables_name = ["p_name", "s_from", "s_to"];

length_index = 0;
diameter_index = 0;
roughness_index = 0; % EPSI

if isfield(graph.Edges, 'Length') && size(graph.Edges.Length, 1) == num_pipes
    pipes_parameters_tab_num_variables = pipes_parameters_tab_num_variables + 1;
    length_index = pipes_parameters_tab_num_variables;
    pipes_parameters_tab_variables_name(length_index) = 'length';
end

if isfield(graph.Edges, 'Diameter') && size(graph.Edges.Diameter, 1) == num_pipes
    pipes_parameters_tab_num_variables = pipes_parameters_tab_num_variables + 1;
    diameter_index = pipes_parameters_tab_num_variables;
    pipes_parameters_tab_variables_name(diameter_index) = 'diameter';
end

if isfield(graph.Edges, 'Epsi') && size(graph.Edges.Epsi, 1) == num_pipes
    pipes_parameters_tab_num_variables = pipes_parameters_tab_num_variables + 1;
    roughness_index = pipes_parameters_tab_num_variables;
    pipes_parameters_tab_variables_name(roughness_index) = 'roughness';
end

pipes_parameters = cell(num_pipes, pipes_parameters_tab_num_variables);
for p = 1:num_pipes
    pipes_parameters{p, 1} = ['Pipe_', num2str(p)];
    pipes_parameters{p, 2} = graph.Edges.EndNodes(p, 1);
    pipes_parameters{p, 3} = graph.Edges.EndNodes(p, 2);
end

if length_index > 0
    for p = 1:num_pipes
        pipes_parameters{p, length_index} = graph.Edges.Length(p);
    end
end

if diameter_index > 0
    for p = 1:num_pipes
        pipes_parameters{p, diameter_index} = graph.Edges.Diameter(p);
    end
end

if roughness_index > 0
    for p = 1:num_pipes
        pipes_parameters{p, roughness_index} = graph.Edges.Epsi(p);
    end
end

pipes_parameters_tab = cell2table(pipes_parameters, ...
    "VariableNames", pipes_parameters_tab_variables_name );


conn = sqlite(db_path, 'connect');
sqlwrite(conn, "pipe_parameters", pipes_parameters_tab);
close(conn);

num_pipes_converted = size(pipes_tab, 1);

end