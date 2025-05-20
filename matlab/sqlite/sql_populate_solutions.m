function [] = sql_populate_solutions(db_path, graph)

if isfield(graph.Nodes, 'PRESSURES')
    num_nodes = size(graph.Nodes.Nodes_ID, 1);
    if (num_nodes < 1)
        return;
    end

    num_variables = 3;
    stations_tab_variables_name = ["s_number", "timestep", "pressure"];

    stations_solution = cell(num_nodes, num_variables);
    for s = 1:num_nodes
        stations_solution{s, 1} = graph.Nodes.Nodes_ID(s);
        for t = 1:size(graph.Nodes.PRESSURES, 2)
            stations_solution{s, 2} = t - 1;
            stations_solution{s, 3} = graph.Nodes.PRESSURES(s, t);
        end
    end

    stations_solution_tab = cell2table(stations_solution, ...
        "VariableNames", stations_tab_variables_name );

    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, "solution_station_pressures", stations_solution_tab);
    close(conn);
end

if isfield(graph.Edges, 'FLOWRATES')
    num_pipes = size(graph.Edges.EndNodes, 1);
    if (num_pipes < 1)
        return;
    end

    pipes_name = map_pipes_name(graph);

    num_variables = 5;
    pipes_tab_variables_name = ["p_name", "s_from", "s_to", "timestep", "flowrate"];

    pipes_solution = cell(num_nodes, num_variables);
    for p = 1:num_pipes
        pipes_solution{p, 1} = pipes_name{p};
        pipes_solution{p, 2} = graph.Edges.EndNodes(p, 1);
        pipes_solution{p, 3} = graph.Edges.EndNodes(p, 2);
        for t = 1:size(graph.Edges.FLOWRATES, 2)
            pipes_solution{p, 4} = t - 1;
            pipes_solution{p, 5} = graph.Edges.FLOWRATES(p, t);
        end
    end

    pipes_solution_tab = cell2table(pipes_solution, ...
        "VariableNames", pipes_tab_variables_name );

    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, "solution_pipe_flowrates", pipes_solution_tab);
    close(conn);
end

end