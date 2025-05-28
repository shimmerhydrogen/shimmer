function [] = sql_populate_boundary_conditions(db_path, graph)

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 2)
    return;
end

num_pipes = size(graph.Edges.EndNodes, 1);

if (num_pipes < 1)
    return;
end


%% Nodes Boundary conditions

db = sqlite(db_path, 'readonly');
station_types_tab = sqlread(db, "station_types");
close(db);

stations_type = map_sql_stations_type(graph);

%% Boundary conditions - ReMi
% ReMi station w/o backflow (set pressure)
% MATLAB 1 trovi il valore in graph.Nodes.PRESSURES
% profiles table with time 0.0

remi_stations = graph.Nodes.Nodes_ID((stations_type == 1));
num_remi_stations = size(remi_stations, 1);

if num_remi_stations > 0
    profiles_remi = cell(num_remi_stations, 3);

    for r = 1:num_remi_stations
        s = remi_stations(r);
        profiles_remi{r, 1} = graph.Nodes.Nodes_ID(s);
        profiles_remi{r, 2} = 0.0;
        profiles_remi{r, 3} = graph.Nodes.PRESSURES(s);
    end

    profiles_remi_tab_variables_name = ["s_number", "prf_time",	"prf_Pset"];
    profiles_remi_tab = cell2table(profiles_remi, ...
        "VariableNames", profiles_remi_tab_variables_name );

    profiles_remi_tab_name = station_types_tab.t_profile_table(1);
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, profiles_remi_tab_name, profiles_remi_tab);
    close(conn);
end

%% Boundary conditions - Inj
% Injection station w/ pressure control (set flow and pressure)
% MATLAB 2/3 trovi il valore in graph.Nodes.G_EXE (G is flow exchanged)
% profiles table with time 0.0

inj_stations = graph.Nodes.Nodes_ID((stations_type == 2));
num_inj_stations = size(inj_stations, 1);

if num_inj_stations > 0
    profiles_inj = cell(num_inj_stations, 4);

    for r = 1:num_inj_stations
        s = inj_stations(r);
        profiles_inj{r, 1} = graph.Nodes.Nodes_ID(s);
        profiles_inj{r, 2} = 0.0;
        profiles_inj{r, 3} = graph.Nodes.PRESSURES(s);
        profiles_inj{r, 4} = graph.Nodes.G_EXE(s);
    end

    profiles_inj_tab_variables_name = ["s_number", "prf_time",	"prf_Pset", "prf_Lset"];
    profiles_inj_tab = cell2table(profiles_inj, ...
        "VariableNames", profiles_inj_tab_variables_name );

    profiles_inj_tab_name = station_types_tab.t_profile_table(2);
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, profiles_inj_tab_name, profiles_inj_tab);
    close(conn);
end

%% Boundary conditions - Cons
% Consumption point w/o pressure control (set flow)
% MATLAB 2/3 trovi il valore in graph.Nodes.G_EXE (G is flow exchanged)
% profiles table with time 0.0

cons_stations = graph.Nodes.Nodes_ID((stations_type == 3));
num_cons_stations = size(cons_stations, 1);

if num_cons_stations > 0
    profiles_cons = cell(num_cons_stations, 3);

    for r = 1:num_cons_stations
        s = cons_stations(r);
        profiles_cons{r, 1} = graph.Nodes.Nodes_ID(s);
        profiles_cons{r, 2} = 0.0;
        profiles_cons{r, 3} = graph.Nodes.G_EXE(s);
    end

    profiles_cons_tab_variables_name = ["s_number", "prf_time", "prf_Lset"];
    profiles_cons_tab = cell2table(profiles_cons, ...
        "VariableNames", profiles_cons_tab_variables_name );

    profiles_cons_tab_name = station_types_tab.t_profile_table(3);
    conn = sqlite(db_path, 'connect');
    sqlwrite(conn, profiles_cons_tab_name, profiles_cons_tab);
    close(conn);
end


%% Pipes Boundary conditions

if isfield(graph.Edges, 'COMP')
    pipes_name = map_pipes_name(graph);

    comp_pipes = find(graph.Edges.COMP == 1);
    num_comp_pipes = size(comp_pipes, 1);

    if num_comp_pipes > 0
        profiles_comp = cell(num_comp_pipes, 3);

        for r = 1:num_comp_pipes
            p = comp_pipes(r);
            profiles_comp{r, 1} = pipes_name{p};
            profiles_comp{r, 2} = graph.Edges.EndNodes(p, 1);
            profiles_comp{r, 3} = graph.Edges.EndNodes(p, 2);
            profiles_comp{r, 4} = 0.0;
            profiles_comp{r, 5} = graph.Edges.COMP_ctrl.RegType(r);
            profiles_comp{r, 6} = graph.Edges.COMP_ctrl.PWR(r);
            profiles_comp{r, 7} = graph.Edges.COMP_ctrl.p_out(r);
            profiles_comp{r, 8} = graph.Edges.COMP_ctrl.p_in(r);
            profiles_comp{r, 9} = graph.Edges.COMP_ctrl.beta(r);
            profiles_comp{r, 10} = graph.Edges.COMP_ctrl.Q(r);
        end

        profiles_comp_tab_variables_name = ["p_name", "s_from", "s_to", "prf_time",	"controlmode", "power", "outpress", "inpress", "ratio", "massflow"];
        profiles_comp_tab = cell2table(profiles_comp, ...
            "VariableNames", profiles_comp_tab_variables_name );

        profiles_comp_tab_name = "compressor_profile";
        conn = sqlite(db_path, 'connect');
        sqlwrite(conn, profiles_comp_tab_name, profiles_comp_tab);
        close(conn);
    end
end


end