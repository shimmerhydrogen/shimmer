clear all;
close all;

%% desidered DB

graph_path = fullfile(pwd, "graphs/test_gasco/gassco_res_small.mat");
graph = load(graph_path);
graph = graph.gasNet_Res;
graph.Nodes.Nodes_ID =  graph.Nodes.Nodes_ID';
graph.Nodes.coordinates_XY = [
    4.8419054048323495, 60.55916692208016
    5.534474571003653, 59.38768321539189
    2.4631164141539705, 57.77672651248003
    7.399360295413453, 53.67057855404977
    2.3139222514284628, 51.03760755736377
    ];

%% original DB

gas_net_path = "/home/geoscore/Dropbox/Polito/ExternalProjects/Shimmer/SHIMMER Project_DISMA/Compare_results/test_gasco/Rete5nodi_input";

gasNet = load(gas_net_path + "/" + "gasNet_input.mat").gasNet;
P_set = load(gas_net_path + "/" + "P_set.mat").P_set;
P_guess = load(gas_net_path + "/" + "Pguess.mat").Pguess;
G_cons = load(gas_net_path + "/" + "Gcons.mat").Gcons;
G_guess = load(gas_net_path + "/" + "Gguess.mat").Gguess;

no_remi_stations = graph.Nodes.Nodes_ID((gasNet.Nodes.Type ~= 1));

%% converted DB

gasNet_Res = gasNet;
gasNet_Res.Nodes.Nodes_ID =  gasNet_Res.Nodes.Nodes_ID';
gasNet_Res.Nodes.Type =  gasNet_Res.Nodes.Type';

gasNet_Res.Nodes.PRESSURES = P_set';
gasNet_Res.Nodes.PRESSURES(no_remi_stations) = P_guess(no_remi_stations);

gasNet_Res.Nodes.G_EXE = G_cons';

gasNet_Res.Edges.FLOWRATES = G_guess';

save(gas_net_path + "/gasNet_converted.mat", "gasNet_Res");
