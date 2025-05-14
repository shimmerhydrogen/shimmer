function [] = sql_populate_limits(db_path, graph)

%% Limits - Nodes

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 1)
    return;
end

db = sqlite(db_path, 'readonly');
station_types_tab = sqlread(db, "station_types");
close(db);

stations_type = map_sql_stations_type(graph);

%% Limits - ReMi
% ReMi station w/o backflow (set pressure)
% Lmin = -300.0;
% Lmax = -10.0;
% Pmin = 60e5;
% Pmax = 80e5;

remi_stations = graph.Nodes.Nodes_ID((stations_type == 1));
num_remi_stations = size(remi_stations, 1);

if num_remi_stations > 0
    limits_remi = cell(num_remi_stations, 5);

    for r = 1:num_remi_stations
        s = remi_stations(r);
        limits_remi{r, 1} = graph.Nodes.Nodes_ID(s);
        limits_remi{r, 2} = -300.0;
        limits_remi{r, 3} = -10.0;
        limits_remi{r, 4} = 60.0e5;
        limits_remi{r, 5} = 80.0e5;
    end

    limits_remi_tab_variables_name = ["s_number", "lim_Lmin", "lim_Lmax", "lim_Pmin", "lim_Pmax"];
    limits_remi_tab = cell2table(limits_remi, ...
        "VariableNames", limits_remi_tab_variables_name );

    limits_remi_tab_name = station_types_tab.t_limits_table(1);
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, limits_remi_tab_name, limits_remi_tab);
    close(conn);
end

%% Limits - Inj
% Injection station w/ pressure control (set flow and pressure)
% Lmin = -300.0;
% Lmax = -10.0;
% Pmin = 60e5;
% Pmax = 80e5;
% f = 1.0;

inj_stations = graph.Nodes.Nodes_ID((stations_type == 2));
num_inj_stations = size(inj_stations, 1);

if num_inj_stations > 0
    limits_inj = cell(num_inj_stations, 6);

    for r = 1:num_inj_stations
        s = inj_stations(r);
        limits_inj{r, 1} = graph.Nodes.Nodes_ID(s);
        limits_inj{r, 2} = -300.0;
        limits_inj{r, 3} = -10.0;
        limits_inj{r, 4} = 60.0e5;
        limits_inj{r, 5} = 80.0e5;
        limits_inj{r, 6} = 1.0;
    end

    limits_inj_tab_variables_name = ["s_number", "lim_Lmin", "lim_Lmax", "lim_Pmin", "lim_Pmax", "parm_f"];
    limits_inj_tab = cell2table(limits_inj, ...
        "VariableNames", limits_inj_tab_variables_name );

    limits_inj_tab_name = station_types_tab.t_limits_table(2);
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, limits_inj_tab_name, limits_inj_tab);
    close(conn);
end

%% Limits - Cons
% Consumption point w/o pressure control (set flow)
% Lmin = -300.0;
% Lmax = -10.0;
% Pmin = 60e5;
% Pmax = 80e5;

cons_stations = graph.Nodes.Nodes_ID((stations_type == 3));
num_cons_stations = size(cons_stations, 1);

if num_cons_stations > 0
    limits_cons = cell(num_cons_stations, 5);

    for r = 1:num_cons_stations
        s = cons_stations(r);
        limits_cons{r, 1} = graph.Nodes.Nodes_ID(s);
        limits_cons{r, 2} = -300.0;
        limits_cons{r, 3} = -10.0;
        limits_cons{r, 4} = 60.0e5;
        limits_cons{r, 5} = 80.0e5;
    end

    limits_cons_tab_variables_name = ["s_number", "lim_Lmin", "lim_Lmax", "lim_Pmin", "lim_Pmax"];
    limits_cons_tab = cell2table(limits_cons, ...
        "VariableNames", limits_cons_tab_variables_name );

    limits_cons_tab_name = station_types_tab.t_limits_table(3);
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, limits_cons_tab_name, limits_cons_tab);
    close(conn);
end

%% Limits - Pipes

if ~isfield(graph.Edges, 'COMP_ctrl_lim')
    return;
end

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    return;
end

pipes_name = map_pipes_name(graph);
pipes_type = map_sql_pipes_type(graph);

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