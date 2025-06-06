function [] = sql_populate_solutions(db_path, graph)

if isfield(graph.Nodes, 'PRESSURES')
    num_nodes = size(graph.Nodes.Nodes_ID, 1);
    if (num_nodes < 1)
        return;
    end

    num_variables = 3;
    stations_tab_variables_name = ["s_number", "timestep", "pressure"];

    num_time_intervals = size(graph.Nodes.PRESSURES, 2);
    stations_solution = cell(num_nodes * num_time_intervals, num_variables);
    ti = 1;
    for s = 1:num_nodes
        for t = 1:size(graph.Nodes.PRESSURES, 2)
            stations_solution{ti, 1} = graph.Nodes.Nodes_ID(s);
            stations_solution{ti, 2} = t - 1;
            stations_solution{ti, 3} = graph.Nodes.PRESSURES(s, t);
            ti = ti + 1;
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

    num_time_intervals = size(graph.Nodes.PRESSURES, 2);
    pipes_solution = cell(num_pipes * num_time_intervals, num_variables);
    ti = 1;
    for p = 1:num_pipes
        for t = 1:size(graph.Edges.FLOWRATES, 2)
            pipes_solution{ti, 1} = pipes_name{p};
            pipes_solution{ti, 2} = graph.Edges.EndNodes(p, 1);
            pipes_solution{ti, 3} = graph.Edges.EndNodes(p, 2);
            pipes_solution{ti, 4} = t - 1;
            pipes_solution{ti, 5} = graph.Edges.FLOWRATES(p, t);
            ti = ti + 1;
        end
    end

    pipes_solution_tab = cell2table(pipes_solution, ...
        "VariableNames", pipes_tab_variables_name );

    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, "solution_pipe_flowrates", pipes_solution_tab);
    close(conn);
end

end