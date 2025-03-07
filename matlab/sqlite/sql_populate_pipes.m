function [num_pipes, num_pipes_converted] = sql_populate_pipes(db_path, graph)

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes == 0)
    num_pipes_converted = 0;
    return;
end

pipes_tab_num_variables = 4;
pipes_tab_variables_name = ["p_name", "s_from", "s_to", "p_type"];

pipes_type = zeros(num_pipes);

for p = 1:num_pipes
    if isfield(graph.Edges, 'PIPE') && size(graph.Edges.PIPE, 1) == num_pipes && graph.Edges.PIPE(p)
        pipes_type(p) = 0;
    elseif isfield(graph.Edges, 'COMP') && size(graph.Edges.COMP, 1) == num_pipes && graph.Edges.COMP(p)
        pipes_type(p) = 2;
    elseif isfield(graph.Edges, 'REDST') && size(graph.Edges.REDST, 1) == num_pipes && graph.Edges.REDST(p)
        pipes_type(p) = 1;
    elseif isfield(graph.Edges, 'VALV') && size(graph.Edges.VALV, 1) == num_pipes && graph.Edges.VALV(p)
        pipes_type(p) = 4;
    elseif isfield(graph.Edges, 'REG') && size(graph.Edges.REG, 1) == num_pipes && graph.Edges.REG(p)
        pipes_type(p) = 3; % TO ASK
    else
        pipes_type(p) = 0; % DEFAULT is 0
    end
end

pipes = cell(num_pipes, pipes_tab_num_variables);
for p = 1:num_pipes
    pipes{p, 1} = ['Pipe_', num2str(p)];
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
roug = 0; % EPSI

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

pipes_parameters_tab = cell2table(pipes_parameters, ...
    "VariableNames", pipes_parameters_tab_variables_name );


conn = sqlite(db_path, 'connect');
sqlwrite(conn, "pipe_parameters", pipes_parameters_tab);
close(conn);

num_pipes_converted = size(pipes_tab, 1);

end