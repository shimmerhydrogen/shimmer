clear;

db_schema = fullfile(pwd, "../../sqlite/shimmer.sql");

%% test_case_1_13

graph_path = fullfile(pwd, "graphs/test_case_1_13/trialNet.mat");
graph = load(graph_path);
graph = graph.gasNet_Res;
db_path = fullfile(pwd, "graphs/test_case_1_13/test_case_1.db");

if exist(db_path, 'file') == 2
  delete(db_path);
end

sql_create(db_path, db_schema);
sql_populate_from_graph(db_path, graph);
graph_test = graph_populate_from_sql(db_path);
assert(isequal(graph.Nodes.Nodes_ID, graph_test.Nodes.Nodes_ID));
assert(isequal(graph.Nodes.Type, graph_test.Nodes.Type));
assert(isequal(graph.Nodes.PRESSURES(:, 1), graph_test.Nodes.PRESSURES));
assert(isequal(graph.Nodes.G_EXE(:, 1), graph_test.Nodes.G_EXE));
assert(isequal(graph.Edges.EndNodes, graph_test.Edges.EndNodes));
assert(isequal(graph.Edges.Length, graph_test.Edges.Length));
assert(isequal(graph.Edges.Diameter, graph_test.Edges.Diameter));
assert(isequal(graph.Edges.Epsi, graph_test.Edges.Epsi));
assert(isequal(graph.Edges.FLOWRATES(:, 1), graph_test.Edges.FLOWRATES));


%% test_case_1

graph_path = fullfile(pwd, "graphs/test_case_1/trialNet.mat");
graph = load(graph_path);
graph = graph.gasNet_Res;
db_path = fullfile(pwd, "graphs/test_case_1/test_case_1.db");

if exist(db_path, 'file') == 2
  delete(db_path);
end

sql_create(db_path, db_schema);
sql_populate_from_graph(db_path, graph);
graph_test = graph_populate_from_sql(db_path);
assert(isequal(graph.Nodes.Nodes_ID, graph_test.Nodes.Nodes_ID));
assert(isequal(graph.Nodes.Type, graph_test.Nodes.Type));
assert(isequal(graph.Nodes.PRESSURES(:, 1), graph_test.Nodes.PRESSURES));
assert(isequal(graph.Nodes.G_EXE(:, 1), graph_test.Nodes.G_EXE));
assert(isequal(graph.Edges.EndNodes, graph_test.Edges.EndNodes));
assert(isequal(graph.Edges.PIPE, graph_test.Edges.PIPE));
assert(isequal(graph.Edges.COMP, graph_test.Edges.COMP));
assert(isequal(graph.Edges.VALV, graph_test.Edges.VALV));
assert(isequal(graph.Edges.REDST, graph_test.Edges.REDST));
assert(isequal(graph.Edges.Length, graph_test.Edges.Length));
assert(isequal(graph.Edges.Diameter, graph_test.Edges.Diameter));
assert(isequal(graph.Edges.Epsi, graph_test.Edges.Epsi));
assert(isequal(graph.Edges.FLOWRATES(:, 1), graph_test.Edges.FLOWRATES));

%% test_gasco

graph_path = fullfile(pwd, "graphs/test_gasco/GasscoGrid.mat");
graph = load(graph_path);
graph = graph.gasNet_Res;
graph.Nodes.Nodes_ID =  graph.Nodes.Nodes_ID';
graph.Nodes.Type =  graph.Nodes.Type';
db_path = fullfile(pwd, "graphs/test_gasco/test_gasco.db");

if exist(db_path, 'file') == 2
  delete(db_path);
end

sql_create(db_path, db_schema);
sql_populate_from_graph(db_path, graph);
graph_test = graph_populate_from_sql(db_path);
assert(isequal(graph.Nodes.Nodes_ID, graph_test.Nodes.Nodes_ID));
assert(isequal(graph.Nodes.Type, graph_test.Nodes.Type));
assert(isequal(graph.Nodes.PRESSURES(:, 1), graph_test.Nodes.PRESSURES));
assert(isequal(graph.Nodes.G_EXE(:, 1), graph_test.Nodes.G_EXE));
assert(isequal(graph.Edges.EndNodes, graph_test.Edges.EndNodes));
assert(isequal(graph.Edges.Length, graph_test.Edges.Length));
assert(isequal(graph.Edges.Diameter, graph_test.Edges.Diameter));
assert(isequal(graph.Edges.Epsi, graph_test.Edges.Epsi));
assert(isequal(graph.Edges.FLOWRATES(:, 1), graph_test.Edges.FLOWRATES));

%% test_inrete

graph_path = fullfile(pwd, "graphs/test_inrete/InreteNet.mat");
graph = load(graph_path);
graph = graph.gasNet_Res;
db_path = fullfile(pwd, "graphs/test_inrete/test_inrete.db");

if exist(db_path, 'file') == 2
  delete(db_path);
end

sql_create(db_path, db_schema);
sql_populate_from_graph(db_path, graph);
graph_test = graph_populate_from_sql(db_path);
assert(isequal(graph.Nodes.Nodes_ID, graph_test.Nodes.Nodes_ID));
assert(isequal(graph.Nodes.Type, graph_test.Nodes.Type));
assert(isequal(graph.Nodes.altitude, graph_test.Nodes.altitude));
assert(isequal(graph.Nodes.coordinates_XY, graph_test.Nodes.coordinates_XY));
assert(isequal(graph.Nodes.PRESSURES(:, 1), graph_test.Nodes.PRESSURES));
assert(isequal(graph.Nodes.G_EXE(:, 1), graph_test.Nodes.G_EXE));
assert(isequal(graph.Edges.EndNodes, graph_test.Edges.EndNodes));
assert(isequal(graph.Edges.Length, graph_test.Edges.Length));
assert(isequal(graph.Edges.Diameter, graph_test.Edges.Diameter));
assert(isequal(graph.Edges.Epsi, graph_test.Edges.Epsi));
assert(isequal(graph.Edges.FLOWRATES(:, 1), graph_test.Edges.FLOWRATES));


%% test_sicilia

graph_path = fullfile(pwd, "graphs/test_sicilia/SicilyNetworkComplete.mat");
graph = load(graph_path);
graph = graph.SicilyNetworkComplete;
graph.Nodes.coordinates_XY = graph.Nodes.Coordinates_XY;
graph.Edges.PIPE(769) = 1;
graph.Edges.COMP(769) = 0;
graph.Edges.REDST(769) = 0;
graph.Edges.VALV(769) = 0;
db_path = fullfile(pwd, "graphs/test_sicilia/test_sicilia.db");

if exist(db_path, 'file') == 2
  delete(db_path);
end

sql_create(db_path, db_schema);
sql_populate_from_graph(db_path, graph);
graph_test = graph_populate_from_sql(db_path);
assert(isequal(graph.Nodes.Nodes_ID, graph_test.Nodes.Nodes_ID));
assert(isequal(graph.Nodes.Type, graph_test.Nodes.Type));
assert(isequal(graph.Nodes.coordinates_XY, graph_test.Nodes.coordinates_XY));
assert(isequal(graph.Nodes.PRESSURES(:, 1), graph_test.Nodes.PRESSURES));
assert(isequal(graph.Nodes.G_EXE(:, 1), graph_test.Nodes.G_EXE));
assert(isequal(graph.Edges.EndNodes, graph_test.Edges.EndNodes));
assert(isequal(graph.Edges.PIPE, graph_test.Edges.PIPE));
assert(isequal(graph.Edges.COMP, graph_test.Edges.COMP));
assert(isequal(graph.Edges.VALV, graph_test.Edges.VALV));
assert(isequal(graph.Edges.REDST, graph_test.Edges.REDST));
assert(isequal(graph.Edges.Length, graph_test.Edges.Length));
assert(isequal(graph.Edges.Diameter, graph_test.Edges.Diameter));
assert(isequal(graph.Edges.Epsi, graph_test.Edges.Epsi));
assert(isequal(graph.Edges.FLOWRATES(:, 1), graph_test.Edges.FLOWRATES));