clear all;
close all;

%% DB schema path
db_schema = fullfile(pwd, "../../sqlite/shimmer.sql");

%% Load NDF Matlab
graph_path = fullfile(pwd, "graph_example.mat");
graph = load(graph_path);
graph = graph.graph;

%% Load NDF SQlite
db_path = fullfile(pwd, "/graph_example.db");
if exist(db_path, 'file') == 2
    delete(db_path);
end

%% Fill NDF from graph Matlab to SQlite
sql_create(db_path, db_schema);
sql_populate_from_graph(db_path, graph);

%% Fill NDF from SQlite to Matlab
converted_graph = graph_populate_from_sql(db_path);

%% Test equality
assert(isequal(graph.Nodes.Nodes_ID, converted_graph.Nodes.Nodes_ID));
assert(isequal(graph.Nodes.Type, converted_graph.Nodes.Type));
assert(isequal(graph.Edges.EndNodes, converted_graph.Edges.EndNodes));
assert(isequal(graph.Edges.Length, converted_graph.Edges.Length));
assert(isequal(graph.Edges.Diameter, converted_graph.Edges.Diameter));
assert(isequal(graph.Edges.Epsi, converted_graph.Edges.Epsi));
assert(isequal(graph.Edges.FLOWRATES, converted_graph.Edges.FLOWRATES));