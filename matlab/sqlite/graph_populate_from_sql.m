function [graph] = graph_populate_from_sql(db_path)

graph = [];

%% Stations
graph = graph_populate_stations(db_path, graph);

%% Pipelines
graph = graph_populate_pipes(db_path, graph);

end