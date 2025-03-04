clear;

graph_path = fullfile(pwd, "graphs/test_case_1/trialNet.mat");
graph = load(graph_path);
graph = graph.gasNet_Res;
db_path = fullfile(pwd, "graphs/test_case_1/test_case_1.db");

clear_database = 1; % Set to zero to NOT clear database on start
reset_database = 1; % Set to zero to NOT clear database on start

%% create db from schema

if (reset_database)
    if exist(db_path, 'file')==2
      delete(db_path);
    end
end

[db_exists, is_empty] = sql_exists(db_path);

if is_empty
    clear_database = 0;
    db_schema = fullfile(pwd, "../../sqlite/shimmer.sql");

    sql_create(db_path, db_schema);
end


%% delete db

if (clear_database == 1)
    conn = sqlite(db_path, "connect");

    execute(conn, "delete from stations");
    execute(conn, "delete from limits_remi_wo");
    execute(conn, "delete from profiles_remi_wo");
    execute(conn, "delete from limits_injection_w");
    execute(conn, "delete from profiles_injection_w");
    execute(conn, "delete from limits_conspoint_wo");
    execute(conn, "delete from profiles_conspoint_wo");
    execute(conn, "delete from pipelines");
    execute(conn, "delete from pipe_parameters");

    close(conn);
end

%% Station Types

conn = sqlite(db_path, "connect");
station_types_tab = sqlread(conn, "station_types");
close(conn);

%% Stations
num_stations_converted = sql_populate_stations(db_path, graph);
assert(num_stations_converted == size(graph.Nodes.Nodes_ID, 1));
