function [] = sql_populate_limits(db_path, graph)

if isfield(graph.Edges, 'COMP_ctrl_lim')
    return;
end

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    return;
end

pipes_name = map_pipes_name(graph);
pipes_type = map_pipes_type(graph);

%% Limits - Compressor
% Compressor (set limits)
% table compressor_limits

compressor_pipes = find((pipes_type == 1))';
num_compressor_pipes = size(compressor_pipes, 1);

if num_compressor_pipes > 0
    limits_comp = cell(num_compressor_pipes, 9);

    for r = 1:num_compressor_pipes
        p = compressor_pipes(r);
        limits_comp{r, 1} = pipes_name{p};
        limits_comp{r, 2} = graph.Edges.EndNodes(p, 1);
        limits_comp{r, 3} = graph.Edges.EndNodes(p, 2);
        limits_comp{r, 4} = graph.Edges.COMP_ctrl_lim.PWR_M(r);
        limits_comp{r, 5} = graph.Edges.COMP_ctrl_lim.p_out_M(r);
        limits_comp{r, 6} = graph.Edges.COMP_ctrl_lim.p_in_m(r);
        limits_comp{r, 7} = graph.Edges.COMP_ctrl_lim.beta_M(r);
        limits_comp{r, 8} = graph.Edges.COMP_ctrl_lim.beta_m(r);
        limits_comp{r, 9} = graph.Edges.COMP_ctrl_lim.Q_M(r);
    end

    limits_comp_tab_variables_name = ["p_name", "s_from",	"s_to",	"max_power", "max_outpress", "min_inpress",	"max_ratio", "min_ratio", "max_massflow"];
    limits_comp_tab = cell2table(limits_comp, ...
        "VariableNames", limits_comp_tab_variables_name );

    limits_comp_tab_name = "compressor_limits";
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, limits_comp_tab_name, limits_comp_tab);
    close(conn);
end


end