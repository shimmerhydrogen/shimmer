function [] = sql_populate_boundary_conditions(db_path, graph)

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 2)
    return;
end

db = sqlite(db_path, 'readonly');
station_types_tab = sqlread(db, "station_types");
close(db);

stations_type = map_stations_type(graph);

%% Boundary conditions - ReMi
% ReMi station w/o backflow (set pressure)
% MATLAB 1 trovi il valore in graph.Nodes.PRESSURES
% MATLAB 2/3 trovi il valore in graph.Nodes.G_EXE (G is flow exchanged)
% Metti questo nelle tabelle profiles con time 0

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
% Injection station w/ pressure control (set flow)
% MATLAB 2/3 trovi il valore in graph.Nodes.G_EXE (G is flow exchanged)
% Metti questo nelle tabelle profiles con time 0
% ASK: why pressure, I set 0.0 on pressure?

inj_stations = graph.Nodes.Nodes_ID((stations_type == 2));
num_inj_stations = size(inj_stations, 1);

if num_inj_stations > 0
    profiles_inj = cell(num_inj_stations, 4);

    for r = 1:num_inj_stations
        s = inj_stations(r);
        profiles_inj{r, 1} = graph.Nodes.Nodes_ID(s);
        profiles_inj{r, 2} = 0.0;
        profiles_inj{r, 3} = 0.0; % ASK pressure
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
% MATLAB 1 trovi il valore in graph.Nodes.PRESSURES
% MATLAB 2/3 trovi il valore in graph.Nodes.G_EXE (G is flow exchanged)
% Metti questo nelle tabelle profiles con time 0

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

end