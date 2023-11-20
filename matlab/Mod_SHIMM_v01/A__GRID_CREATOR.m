clear all
clc
close all
%
INPUT_PIPES = xlsread('INPUT_Pambour1.xlsx','PIPE');
INPUT_NODES = xlsread('INPUT_Pambour1.xlsx','PIPE_N','B3:Z5');

% TIME=0.25*3600;%s
TIME = 2*3600;
dt   = 1*60;
tt   = [0:dt:TIME];
dimt = length(tt);


Ainput_inner = INPUT_PIPES(:,1:3);% |--n Branch--|--NodeIN--|--NodeOUT--|
LL_inner  = INPUT_PIPES(:,4)*1e3; %m
DD_inner  = INPUT_PIPES(:,5); 	%m
epsi_inner= INPUT_PIPES(:,6)*1e-3;%m
COMP_inner= INPUT_PIPES(:,7);
REG_inner = INPUT_PIPES(:,8);
VAL_inner = INPUT_PIPES(:,9);
RES_inner = INPUT_PIPES(:,10);
nGridPoints = INPUT_PIPES(:,11);

HH_inner = INPUT_NODES(:,2);%m

%  dimb0=size(INPUT_PIPES,1);
%  dimn0=max(max(INPUT_PIPES(:,2:3)));
%
% BR0 = INPUT_PIPES(:,1);
% IN0 =  INPUT_PIPES(:,2);
% OUT0 = INPUT_PIPES(:,3);
% Asp_p0=sparse(IN0,BR0,ones(1,dimb0),dimn0,dimb0);
% Asp_m0=sparse(OUT0,BR0,-ones(1,dimb0),dimn0,dimb0);
% Asp0=Asp_p0+Asp_m0;
%
% delH_inner=-Asp0'*HH_inner;

IN = [];
OUT =[];
CORR = [];
nodi_notevoli = [];
nodi_notevoli_corr = [];
INPUT_b = [];
% INPUT_n = [];
COMP = [];
REG = [];
VAL = [];
RES = [];
NP  = [];
PIPE= [];
in = 0;

%WK:  numGridPoints? is the number of grid points in a pipe? number of nodes should
% be equal to numGridPoints ( but it seems it is numGridPoints +1).
% A wider example could help tio understand which are the right dimensions here.

for kk = 1:length(INPUT_PIPES(:,1))
	branch(kk).pos = zeros(1, nGridPoints(kk)+1);
	branch(kk).nodes = zeros(1,nGridPoints(kk)+1);

	branch(kk).pos  = [0:LL_inner(kk)./nGridPoints(kk):LL_inner(kk)];
	%nodi notevoli
	branch(kk).nodes(1) = INPUT_PIPES(kk,2); %in
	branch(kk).nodes(end) = INPUT_PIPES(kk,3); %out
	%nodi
	branch(kk).unkown = [in+1 : in+length(branch(kk).xx)];

	if kk>=2
		if any(branch(kk).nodes(1)==nodi_notevoli) %ingresso
			pos = find(branch(kk).nodes(1)==nodi_notevoli);
			branch(kk).xx(3,1) = nodi_notevoli_corr(pos(1));
			pos = [];
			flagIN = 1;
		end

		if any(branch(kk).nodes(end)==nodi_notevoli) %uscita
			pos = find(branch(kk).nodes(end)==nodi_notevoli);
			branch(kk).xx(3,end)=nodi_notevoli_corr(pos(end));
			pos = [];
			flagOUT = 1;
		end

		if flagIN==1 & flagOUT==1
			branch(kk).xx(3,2:end-1) = [in+1:in+size(branch(kk).xx,2)-2];
			in = branch(kk).xx(3,end-1);
		elseif flagIN==1 & flagOUT==0
			branch(kk).xx(3,2:end) = [in+1:in+size(branch(kk).xx,2)-1];
			in = branch(kk).xx(3,end);
		elseif flagIN==0 & flagOUT==1
			branch(kk).xx(3,1:end-1) = [in+1:in+size(branch(kk).xx,2)-1];
			in = branch(kk).xx(3,end-1);
		else
			branch(kk).xx(3,1:end-1) = [in+1:in+size(branch(kk).xx,2)-1];
			in = branch(kk).xx(3,end);
		end
	else
		in = branch(kk).xx(3,end);

	end



	nodi_notevoli = [nodi_notevoli; branch(kk).nodes(:)'];
	nodi_notevoli_corr = [nodi_notevoli_corr; branch(kk).xx(3,:)'];
	NN = [nodi_notevoli, nodi_notevoli_corr];

	flagIN  = 0;
	flagOUT = 0;

	% NODE_BR(kk).IN  = branch(kk).xx(3,1:end-1);
	% NODE_BR(kk).OUT = branch(kk).xx(3,2:end);
	IN  = [IN;  branch(kk).xx(3,1:end-1)'];
	OUT = [OUT; branch(kk).xx(3,2:end)'];


	lung(kk)   = length(branch(kk).xx);
	branch(kk).dx = LL_inner(kk)./nGridPoints(kk);
	branch(kk).diam = DD_inner(kk);
	branch(kk).epsi = epsi_inner(kk);

	branch(kk).xx(7,:) = linspace(HH_inner(INPUT_PIPES(kk,2)),HH_inner(INPUT_PIPES(kk,3)),nGridPoints(kk)+1);

	altitude(branch(kk).xx(3,:)) = branch(kk).xx(7,:);

	INPUT_b = [INPUT_b; branch(kk).xx(3,1:end-1)',branch(kk).xx(3,2:end)',branch(kk).xx(4:6,2:end)'];

	if COMP_inner(kk)==1
		COMP = [COMP;branch(kk).xx(3,1:end-1)'];
		NP   = [NP;  branch(kk).xx(3,1:end-1)'];
	elseif REG_inner(kk)==1
		REG  = [REG; branch(kk).xx(3,1:end-1)'];
		NP	 = [NP;  branch(kk).xx(3,1:end-1)'];
	elseif VAL_inner(kk)==1
		VAL  = [VAL; branch(kk).xx(3,1:end-1)'];
		NP   = [NP;  branch(kk).xx(3,1:end-1)'];
	elseif RES_inner(kk)==1
		RES  = [RES; branch(kk).xx(3,1:end-1)'];
		NP   = [NP;  branch(kk).xx(3,1:end-1)'];
	else
		PIPE = [PIPE; branch(kk).xx(3,1:end-1)'];
	end

	xx_old(kk).xx = branch(kk).pos(1,1:end-1);

end
dimb = size(INPUT_b,1);
BR = [1:1:dimb]';
INPUT_b = [BR,INPUT_b];
dimn = max(max(INPUT_b(:,2:3)));

INPUT_n = zeros(dimn,25);
corrispondenze = unique(NN(find(NN(:,1)~=0),:),'rows');
INPUT_n(:,1) = [1:1:dimn];
INPUT_n(corrispondenze(:,2),3:end) = INPUT_NODES(corrispondenze(:,1),3:end);
INPUT_n(:,2) = altitude;

Ainput=INPUT_b(:,1:3);

BR = INPUT_b(:,1);
IN =  INPUT_b(:,2);
OUT = INPUT_b(:,3);
Asp_p = sparse(IN,BR,ones(1,dimb),dimn,dimb);
Asp_m = sparse(OUT,BR,[-ones(1,dimb)],dimn,dimb);
Asp = Asp_p + Asp_m;

%figure(100)
%GRAPH = digraph(IN,OUT);
%plot(GRAPH,'layout','force');

%GRAPH1 = graph(IN,OUT);
%gradi  = degree(GRAPH1);
%INNER_H = find(gradi>2);


save('DATA_INPUT_Trial2023_3.mat')



