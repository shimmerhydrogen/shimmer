clear all;
close all;

%% original DB

gas_net_path = "/home/geoscore/Dropbox/Polito/ExternalProjects/Shimmer/SHIMMER Project_DISMA/Compare_results/test_gasco/Rete10pipe_EACH";

gasNet = load(gas_net_path + "/" + "gasNet_input.mat").gasNet;
P_set = load(gas_net_path + "/" + "P_set.mat").P_set;
P_guess = load(gas_net_path + "/" + "Pguess.mat").Pguess;
G_cons = load(gas_net_path + "/" + "Gcons.mat").Gcons;
G_guess = load(gas_net_path + "/" + "Gguess.mat").Gguess;

no_remi_stations = gasNet.Nodes.Nodes_ID((gasNet.Nodes.Type ~= 1));

%% converted DB

gasNet_Res = gasNet;
gasNet_Res.Nodes.Nodes_ID =  gasNet_Res.Nodes.Nodes_ID';
gasNet_Res.Nodes.Type =  gasNet_Res.Nodes.Type';

gasNet_Res.Nodes.PRESSURES = P_set';
gasNet_Res.Nodes.PRESSURES(no_remi_stations) = P_guess(no_remi_stations);

gasNet_Res.Nodes.G_EXE = G_cons';

gasNet_Res.Edges.FLOWRATES = G_guess';

save(gas_net_path + "/gasNet_converted.mat", "gasNet_Res");
