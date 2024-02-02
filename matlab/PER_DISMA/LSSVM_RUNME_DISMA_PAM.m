%RUNME

close all
clear all
clc

% acquisizione dati da excel
nodei=xlsread('small_network_disma_pipes.xlsx','B2:B16'); 
nodeo=xlsread('small_network_disma_pipes.xlsx','C2:C16');
Diametri=xlsread('small_network_disma_pipes.xlsx','E2:E16'); %m
L=xlsread('small_network_disma_pipes.xlsx','D2:D16'); %m
epsi=xlsread('small_network_disma_pipes.xlsx','G2:G16'); %m
S=pi*(Diametri.^2)/4;

% create graph
% gasNet = graph(nodei, nodeo);
gasNet.Edges.EndNodes(:,1)=nodei; 
gasNet.Edges.EndNodes(:,2)=nodeo; 
gasNet.Edges.Length = L;
gasNet.Edges.Diameter = Diametri;
gasNet.Edges.Epsi = epsi;

gasNet.Nodes.Nodes_ID = [1:13]';
gasNet.Nodes.deg = [0; 1; 1; 0; 0; 1; 1; 0; 1; 1; 1; 1; 1];

grafo = graph(nodei, nodeo);

figure
plot(grafo)

% acquisizione dati consumo e pressione
Gcons = xlsread('small_network_disma_nodes.xlsx','Consumpt','B3:B15'); 
Pset = 70; % bar
Gguess = xlsread('small_network_disma_pipes.xlsx', 'F2:F16');
Pguess = xlsread('small_network_disma_nodes.xlsx','Consumpt', 'C3:C15');

% acquisizione profili di consumo
Profiles = xlsread('small_network_disma_nodes.xlsx','Profiles', 'B3:N27')';
% and the boundary conditions 
gasNet.Nodes.isSource(1)=1;

[nn,TIME]=size(Profiles); %in hours
dt_min=60; %min
dt=dt_min*60;

Profiles_t=interp1([0:1:TIME-1],Profiles([1:end],:)',[0:dt/3600:TIME-1])';
tt=[0:dt:(TIME-1)*3600];
Gcons_t=Gcons.*Profiles_t;

Pset_t=Pset.*Profiles_t(1,:);


%% FUNCTION CALLING
[gasNet_Res, pos] = PambourSS_DISMA(gasNet,Gcons_t(:,1),Pset_t(1),Gguess,Pguess);

%% OUTPUT

RES_P(1,:) = gasNet_Res.Nodes.PRESSURES./10^5;%bar
RES_G(1,:) = gasNet_Res.Edges.FLOWRATES;%kg/s
RES_G_exe(1,:) = gasNet_Res.Nodes.G_EXE;%kg/s

%% FUNCTION CALLING TRANSIENT
[gasNet_Res, pos] = PambourTS_DISMA(gasNet,Gcons_t(1:end,:),Pset_t(1:end,:),RES_G(1,:),RES_P(1,:),dt);

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

figure
plot(RES_P_T')
title('Nodal Pressures - bar')
figure
plot(RES_G_T')
title('Mass Flows - kg/s')
figure
plot(RES_G_exe_T')
title('Mass Flows IN-OUT - kg/s')