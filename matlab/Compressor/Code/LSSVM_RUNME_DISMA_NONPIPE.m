%RUNME

close all
clear all
clc

% acquisizione dati da excel
nodei=xlsread('small_network_disma_pipes2.xlsx','B2:B18'); 
nodeo=xlsread('small_network_disma_pipes2.xlsx','C2:C18');
Diametri=xlsread('small_network_disma_pipes2.xlsx','E2:E18'); %m
L=xlsread('small_network_disma_pipes2.xlsx','D2:D18'); %m
epsi=xlsread('small_network_disma_pipes2.xlsx','G2:G18'); %m
S=pi*(Diametri.^2)/4;

% acquisizione non pipe elements
TYPE=xlsread('small_network_disma_branches_nonpipe.xlsx','E3:E19');
COMP=find(TYPE==1);
REDST=find(TYPE==2);
VALV=find(TYPE==3);
PIPE=find(TYPE==0);

% create graph
% gasNet = graph(nodei, nodeo);
gasNet.Edges.EndNodes(:,1)=nodei; 
gasNet.Edges.EndNodes(:,2)=nodeo; 
gasNet.Edges.Length = L;
gasNet.Edges.Diameter = Diametri;
gasNet.Edges.Epsi = epsi;
gasNet.Edges.COMP = zeros(size(TYPE));
gasNet.Edges.COMP(COMP)=1;
gasNet.Edges.REDST = zeros(size(TYPE));
gasNet.Edges.REDST(REDST)=1;
gasNet.Edges.VALV = zeros(size(TYPE));
gasNet.Edges.VALV(VALV)=1;
gasNet.Edges.PIPE = zeros(size(TYPE));
gasNet.Edges.PIPE(PIPE)=1;

if length(COMP)>0
ctrl_mode_COMP=readcell('small_network_disma_branches_nonpipe.xlsx',"Sheet", 'COMP','Range', 'A3:N3');
COMP_ctrl = table(ctrl_mode_COMP(:,1),ctrl_mode_COMP(:,3),ctrl_mode_COMP(:,4),ctrl_mode_COMP(:,5),ctrl_mode_COMP(:,6),ctrl_mode_COMP(:,7),ctrl_mode_COMP(:,8),...
                'VariableNames',["ID","RegType","PWR","p_out","p_in","beta","Q"]);
COMP_limits = table(ctrl_mode_COMP(:,1),ctrl_mode_COMP(:,3),ctrl_mode_COMP(:,9),ctrl_mode_COMP(:,10),ctrl_mode_COMP(:,11),ctrl_mode_COMP(:,12),ctrl_mode_COMP(:,13),ctrl_mode_COMP(:,14),...
                'VariableNames',["ID","RegType","PWR_M","p_out_M","p_in_m","beta_M","beta_m","Q_M"]);

gasNet.Edges.COMP_ctrl=COMP_ctrl;
gasNet.Edges.COMP_ctrl_lim=COMP_limits;
end


gasNet.Nodes.Nodes_ID = [1:15]';
% gasNet.Nodes.deg = [0; 0; 0; 0; 0; 1; 0; 0; 1; 0; 1; 0; 0];

grafo = graph(nodei, nodeo);

figure
plot(grafo)

% acquisizione dati consumo e pressione
Gcons = xlsread('small_network_disma_nodes_nonpipe.xlsx','Consumpt','M3:M17'); 
Pset = xlsread('small_network_disma_nodes_nonpipe.xlsx','Consumpt','F3:F17');  % bar
Gguess = xlsread('small_network_disma_pipes2.xlsx', 'F2:F18');
Pguess = xlsread('small_network_disma_nodes_nonpipe.xlsx','Consumpt', 'N3:N17');

% acquisizione profili di consumo
Profiles = xlsread('small_network_disma_nodes_nonpipe.xlsx','Profiles', 'B3:P27')';
% and the boundary conditions 
gasNet.Nodes.Type = xlsread('small_network_disma_nodes_nonpipe.xlsx', 'C3:C17');
gasNet.Nodes.Pset_bc = xlsread('small_network_disma_nodes_nonpipe.xlsx', 'F3:F17');
gasNet.Nodes.Lset_bc = xlsread('small_network_disma_nodes_nonpipe.xlsx', 'G3:G17');

[nn,TIME]=size(Profiles); %in hours
dt_min=60; %min
dt=dt_min*60;

Profiles_t=interp1([0:1:TIME-1],Profiles([1:end],:)',[0:dt/3600:TIME-1])';
tt=[0:dt:(TIME-1)*3600];
Gcons_t=Gcons.*Profiles_t;

Pset_ss=Pset.*Profiles_t(:,1);
Pset_ts = Pset.*Profiles_t(:,:);


%% FUNCTION CALLING
[gasNet_Res, pos] = PambourSS_DISMA_NONPIPE(gasNet,Gcons_t(:,1),Pset_ss,Gguess,Pguess);

save("trialNet", "gasNet_Res");

%% OUTPUT

RES_P(1,:) = gasNet_Res.Nodes.PRESSURES./10^5;%bar
RES_G(1,:) = gasNet_Res.Edges.FLOWRATES;%kg/s
RES_G_exe(1,:) = gasNet_Res.Nodes.G_EXE;%kg/s

%% FUNCTION CALLING TRANSIENT
[gasNet_Res, pos, BC_change_rec] = PambourTS_DISMA_NONPIPE_revMC(gasNet,Gcons_t(1:end,:),Pset_ts(:,:),RES_G(1,:),RES_P(1,:),dt);

%% OUTPUT

RES_P_T = gasNet_Res.Nodes.PRESSURES./10^5;%bar
RES_G_T = gasNet_Res.Edges.FLOWRATES;%kg/s
RES_G_exe_T = gasNet_Res.Nodes.G_EXE;%kg/s

RES_P_T(:,1) = RES_P(1,:)';%bar
RES_G_T(:,1) = RES_G(1,:)';%kg/s
RES_G_exe_T(:,1) = RES_G_exe(1,:)';%kg/s

% 
figure
bar(RES_P)
title('Nodal Pressures - bar')

figure
bar(RES_G)
title('Mass Flows - kg/s')

% Plot nodal pressures with dynamic legend
figure;

% Main plot
plot(0:size(RES_P_T, 2)-1, RES_P_T', LineWidth = 2); % Plot rows of RES_P_T against time steps
hold on;

% Dynamically generate legend labels
numNodes = size(RES_P_T, 1); % Number of rows (nodes) in RES_P_T
legendLabels = arrayfun(@(n) sprintf('Node %d', n), 1:numNodes, 'UniformOutput', false);

% Title, grid, and legend
title('Nodal Pressures - bar');
grid on;
legend(legendLabels, 'Location', 'northeastoutside');
% Plot data with custom legend entries
figure;

% Main plot
plot(RES_G_T', LineWidth = 2);

% Dynamically generate legend labels
numPipes = size(RES_G_T, 1); % Number of rows in RES_G_T
legendLabels = arrayfun(@(p) sprintf('Pipe %d', p), 1:numPipes, 'UniformOutput', false);

% Title, grid, and legend
title('Mass Flows - kg/s');
grid on;
legend(legendLabels, 'Location', 'northeastoutside');

% Plot data
figure;
hold on;

% Main plot
plot([0:size(RES_P_T,2)-1], RES_G_exe_T, LineWidth= 2);

% Dynamically generate legend labels
numNodes = size(RES_G_exe_T, 1); % Number of rows in RES_G_exe_T
legendLabels = arrayfun(@(n) sprintf('Node %d', n), 1:numNodes, 'UniformOutput', false);

% Apply legend
legend(legendLabels, 'Location', 'northeastoutside');

% Title and grid
title('Mass Flows IN-OUT - kg/s');
grid on;

% % 
% figure
% bar(RES_P)
% title('Nodal Pressures - bar')
% 
% figure
% bar(RES_G)
% title('Mass Flows - kg/s')
% 
% COLOR=jet(size(RES_G_T,1))
% figure
% plot([0:size(RES_P_T,2)-1],RES_P_T')%,'color',COLOR)
% set(gca, 'ColorOrder', COLOR)
% hold on
% legend()
% legend('Location','northeastoutside')
% plot([0:size(RES_P_T,2)-1],RES_P_T(1,:)','color','k')
% plot([0:size(RES_P_T,2)-1],RES_P_T(12,:)','color','m')
% title('Nodal Pressures - bar')
% grid on
% 
% figure
% plot(RES_G_T')
% set(gca, 'ColorOrder', COLOR)
% title('Mass Flows - kg/s')
% grid on
% legend()
% legend('Location','northeastoutside')
% 
% figure
% plot([0:size(RES_P_T,2)-1],RES_G_exe_T')
% set(gca, 'ColorOrder', COLOR)
% hold on
% legend()
% legend('Location','northeastoutside')
% plot([0:size(RES_P_T,2)-1],RES_G_exe_T(1,:)','color','k')
% plot([0:size(RES_P_T,2)-1],RES_G_exe_T(12,:)','color','m')
% title('Mass Flows IN-OUT - kg/s')
% grid on
