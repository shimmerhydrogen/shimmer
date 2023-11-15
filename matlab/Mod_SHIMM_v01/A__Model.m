clear all
close all
clc

check_constr=0;
A__GRID_CREATOR
toll = 1e-4;  %tolerance on Residue calculation of linearized fluid-dynamic problem
% load('DATA_INPUT_Trial2023.mat')

Gas_Prof0=xlsread('INPUT_Pambour1.xlsx','Profili_M2','C3:RO5')';


TIME=length(Gas_Prof0(:,1))*180-1;
dt=180; %s
tt=[0:dt:TIME];
dimt=length(tt);

%this command for possible interpolation
%Gas_Prof(:,1)=interp1([0:3600:TIME],Gas_Prof0(:,1),tt)';
Gas_Prof=Gas_Prof0;

time0=cputime;

A  = Asp;
LL = INPUT_b(:,4);            %m  %length of the pipelines
dx = LL;
dxe = dx;                         % vectoor size: #pipelines
HH =INPUT_n(:,2);             %m  % altitude of each node - size: #node
delH = -A'*HH;                %m  % h_out-h_in / delta height %m - size: #pipeline
DD = INPUT_b(:,5);            %m  % Diameters of the pipelines
epsi = INPUT_b(:,6);          %m  % absolute roughness
AA = pi.*DD.^2/4;             %m2 % pipeline cross section

eta = 1e-5;                   %kg/m/s % dynamic viscosity
RD1 = 0.6;                    %[-]    % Relative Density
rho_air_std = 1.225;          %kg/Sm3 @p=1 bar T=15 C - 288K
%rho_air_std = 1.292;         %kg/Nm3 @p=1 bar T=0 C - 273K

rho_g_std = RD1*rho_air_std;  %kg/Sm3  % density of natural gas @standard cond

G_ext = INPUT_n(:,3).*0.77;       %kg/s   % vector size: #nodes - external exchange of gas mass flows (-) inwards / (+) outwards.
P_in  = INPUT_n(:,4).*1e6+101325; %Pa     % vector size: #nodes - nodal pressures -it is basically a tentative value...
%%
% P_in(1)=120.*1e5+101325; %Pa         %... except at nodes with fixed pressures where it is set as a boudary condition
                                       % NB!! MODIFICATION REQUIRED
                                       % we should provide the list of
                                       % nodes where the pressure is fixed.
%%
INLET  = find(G_ext<0);  %nodi di ingresso
OUTLET = find(G_ext>0);  %nodi di uscita
INNER  = find(G_ext==0); %nodi interni - non scambiano con l'esterno

%% Input concentrazioni
MOLEFRAC_INPUT = INPUT_n(INLET,5:25);  % udm: %;  initial composition at nodes that are INLETS.
nComponents = size(MOLEFRAC_INPUT,2);  % size= #columns of MOLEFRAC_INPUT = #chemical species

MASSFRAC_INPUT = MOL2MASS_CONV(MOLEFRAC_INPUT);  % udm: %
                                               % conversion form molar fraction to mass fraction
                                               % by means of function "MOL2MASS_CONV

UGC = 8.3144598*1e3;                  %J/kmol/K  % Universal Gas Constant
MM = MolarMassGERG1(reshape(MOLEFRAC_INPUT,[size(INLET,1),nComponents])); %kg/kmol
                                                                        % calculation of the molar mass of gas(es) at the inlet.
RR = UGC./MM'; %J/kg/K  % gas constant of the gas at the inlet


%% profiles

%%% GAS FLOWS PROFILES (IN ENERGY TERMS - beacuse of the mutli-component feature)

G_ext_t1 = zeros(length(G_ext),dimt);         % vector size: #nodes x #timesteps
                                              % initialization of gas flow time-profiles for each node of the network
G_ext_t1(corrispondenze(:,2),:) = Gas_Prof'.*0.77;  % kg/s (kWh/h)

%%% PRESSURES PROFILES

profiloP = 1.*(tt<=TIME);      % vector size: #timesteps
                               % "per unit" profile of the pressure -
                               % alway constant and equal to the one fixed
                               % above (line 47)

%%% COMPOSITION
% molar profiles
profiloCONC = zeros(dimn,nComponents,length(tt));
profiloCONC(INLET,:,:) = MOLEFRAC_INPUT(1,:).*ones(length(INLET),nComponents,length(tt));
% we are constructing a 3D matrix #inlets x #components x #timesteps
% where theoretically we could shape a profile of varying molar share of
% each component of the INLET gas at every timesteps
if length(INLET)>1
	profiloCONC(INLET(2),:,1:end) = MOLEFRAC_INPUT(2,:).*ones(length(INLET(2)),nComponents,length(tt(1:end)));
end

% mass profiles  %same as above but with mass fractions
profiloCONC_M = zeros(dimn,nComponents,length(tt));
profiloCONC_M(INLET,:,:) = MASSFRAC_INPUT(1,:).*ones(length(INLET),nComponents,length(tt));
if length(INLET)>1
	profiloCONC_M(INLET(2),:,1:end) = MASSFRAC_INPUT(2,:).*ones(length(INLET(2)),nComponents,length(tt(1:end)));
end

% figure(3)
% ComponentsM1= reshape(profiloCONC_M(INLET(1),:,:),nComponents,length(tt))';
% plot(tt/3600,ComponentsM1,'linewidth',3);

%% BC - BOUNDARY CONDITIONS
% fixed G pattern at the outlets

%it is actually the vector/matrix G_ext_t1
figure
plot(tt/3600,Gas_Prof(:,[2,3]))
grid on
xlabel('time [h]')
ylabel('mass flow rate [kg/s]')
title('outlet mass flow rate')
legend('G_e_x_t _2', 'G_e_x_t _3')
% fixed P pattern at inlet where the pressure is set
p_in_t = P_in(INLET(1)).*profiloP;
figure
plot(tt/3600,p_in_t/10^6)
grid on
xlabel('time [h]')
ylabel('pressure [MPa]')
title('Inlet node fixed pressure')
% concentration pattern at inlet set already in the section before


%% STD Conditions
p_std = 1.01325*1e5; %Pa
T_std = 273.15 + 15; %K
% isothermal conditions
Tn = (273.15 + 20)*ones(dimn,1);% K nodal temperature
Tb = (273.15 + 20)*ones(dimb,1);% K pipes temperature


%% IC - Initial Conditions
% @ t=0 hyp: CONSTANT MASSFLOW
p_0 = zeros(dimn,dimt);
p_0(1,:) = p_in_t; %nodal pressure status of all nodes (tentative/guessed value)

%concentration
CC_0(INLET,:) = profiloCONC(INLET,:,1);
if length(INNER)~=0
	CC_0(INNER,:) = padarray(profiloCONC(INLET(1),:,1),length(INNER)-1,'symmetric','pre');
end
CC_0(OUTLET,:) = padarray(profiloCONC(INLET(1),:,1),length(OUTLET)-1,'symmetric','pre');
% we are making the assumption that at the Initial Condition the gas flow
% concentration is the same as the one at the INLET (1)
% udm: % - out of 100

% here we are preparing the 3D matrix that will store the results in terms
% of concentration of all the networks node
CC_gas = zeros(dimn,nComponents,length(tt));
CC_gas(:,:,1) = CC_0./100;                            % udm: [-] - out of 1
CC_gas(INLET,:,:) = profiloCONC(INLET,:,:)./100;
% CC_gas(INNER,:,1) = padarray(profiloCONC(INLET(1),:,1),length(INNER)-1,'symmetric','pre')./100; %istante iniziale
CC_gas(OUTLET,:,1) = padarray(profiloCONC(INLET(1),:,1),length(OUTLET)-1,'symmetric','pre')./100; %istante iniziale

% here we are preparing the 3D matrix that will store the results in terms
% of concentration specifically for the nodes that has an exchange with the
% outside of the network
CC_gas_ext = zeros(dimn,nComponents,length(tt));
CC_gas_ext([INLET;OUTLET],:,1) = CC_gas([INLET;OUTLET],:,1);
%these are useful to quantify the mass flow rates as we are working with
%fixed thermal energy request (kWh/h) at outlet nodes, not directly mass or
%volume flow rates


%here we are creating the same variables but with mass based concentration.
CC_0M = MOL2MASS_CONV(CC_0);  % !!!NB out of 100!!
CC_gasM = MOL2MASS_CONV3(CC_gas);
CC_gasM_ext = MOL2MASS_CONV3(CC_gas_ext);

%condizione iniziale concentrazione specie 1 in tutti i nodi e in tuttii
%rami
for kk=1:INPUT_b_inner(end,1)
	vel_b(kk).ccb = CC_gasM(branch(kk).xx(3,2:end)-1,:,1); %imputo la concentrazione del nodo precedente al pezzo di ramo (batch) successiva)
	vel_n(kk).C_N(:,:,1) = CC_gasM(branch(kk).xx(3,:),:,1);
end


MM = MolarMassGERG1(reshape(CC_gas(:,:,1),dimn,21)*100); %kg/kmol  % size: 1 x #nodes
% calculation of the molar mass of gas for all the nodes of the network
% at the initial timestep.

RR = UGC./MM'; %J/kg/K   % size: #nodes x 1
% gas constant of the gas at the inlet for all the nodes of the network
% at the initial timestep.

MHV(:,1) = MASS_HV(reshape(CC_gasM(:,:,1),[dimn 21])/100)'; % MJ/kg
% mass based higher heating value of the natural gas calculated according
% to the mass concentration at all the nodes at the first timestep

H_ext_t = G_ext_t1.*(MHV(:,1)*1000); %kWh/h  % THERMAL ENERGY DEMANDS/SUPPLY
% H_ext_t=G_ext_t1; %kWh/h  % THERMAL ENERGY DEMANDS/SUPPLY
% G_ext_t1=H_ext_t./(MHV(:,1)*1000); %kW * kg/MJ = kg/s
%calculation of GAS MASS FLOW RATES based on the composition at the first
%timestep detailed for all the nodes.

s = 0;  % assumption: we neglect the gravitational terms in this phase in which
% we want to calculate a first good approximation of the 1st timestep
% fluid-dynamic equilibrium of the network
Aplus  = (A==1).*1;
Aminus = (A==-1).*1;
Aminus_s0 = Aminus %Aminus.*repmat(exp(s/2),1,dimn)'; %WK: these dimensions are not consistent
ADP_0 = (-Aminus_s0 + Aplus)';
% %
% function that uses another type gas network model (steady state) used to
% calculate a good guessed value for the solution of nodal pressure and pipeline
% mass flow rates

[p_0 G_0] = SteadyState(A,Aplus,Aminus,P_in, Tb, G_ext_t1(:,1),AA, DD, dx, HH, epsi, RR, MM, CC_gas, NP, PIPE,	dimn, dimb, dimt);

% % % [p_0 G_0] = Q_linearization_fun_08(A,Aplus,Aminus,p_0(1,1),Tb,G_ext_t1(:,1),AA,DD,dxe,epsi,eta,RR,dt,p_0(:,1),G_0(:,1),dimn,Aplus'*reshape(CC_gas(:,:,1),dimn,21));
save('Trial2023_3_P0','p_0');
save('Trial2023_3_G0','G_0');

%   % I am saving this so that for further variations on the same initial
%   % condition we can simply load the result with the command commented
%   % below

%   % load('Greenstream_BASEmA_SUM_p0_120_120','p_0');
%   % load('Greenstream_BASEmA_SUM_G0_120_120','G_0');

%   % G_0(:,:)=40*ones(dimn,dimt);
%   % p_0(:,1)=P_in(1)*ones(dimn,1);


%%

G_n = G_0(:,1);            % matrix (bxt) of pipeline mass flow rates for each timestep
p_n = p_0(:,1);            % matrix (nxt) of nodal pressures for each timestep
L_0(:,1) = G_ext_t1(:,1);     % matrix (nxt) of nodal inlet(-)/outlet(+) mass flow rates for each timestep
                              % quasi-fixed term: it is mass flow rated fixed
                              % by the THERMAL REQUEST thus depends on the
                              % gas quality at each node


ii = 1;       % time loop counter
k = 1;        % linearization loop counter

%initialization of the linearization loop (k-referred)
p_k(:,k) = p_0(:,ii);
Pin  = Aplus' *p_k(:,k);
Pout = Aminus'*p_k(:,k);
pm(:,1) = (2.0/3.0) * averagePressure(Pin, Pout);
MM = MolarMassGERG1(reshape(CC_gas(:,:,1),dimn,21)*100);%kg/kmol
RR = UGC./MM'; %J/kg/K
MolMass(:,1) = MM;

%initialization of the main parameters related to the
%equation of state used in this case
%EoS: GERG-2008

%NODAL BASED CALCULATION
iFlag = 0;
%[Den, ierr, herr] = DensityGERG(iFlag, Tn, p_0(:,1)/1e3,reshape(CC_gas(:,:,1),dimn,21)); not of use
x_ni = reshape(CC_gas(:,:,1),dimn,21);
gerg_ni = UtilitiesGERG(x_ni, dimn );

% Equation of State: it gives us the value of Zm
% (compressibility factor) and the density related to the
% pressure given (in this case the nodal pressure)
% NB the pressure should be given in kPa
[Zm, Den] = PropertiesGERG(iFlag, p_0(:,1)/1e3, Tn, x_ni, dimn, gerg_ni);

%NB: this section can be simplified by choosing another (and
%simpler Equation of State or can be generalized implementing
%the choice among different equation of state


ZZn  = Zm;            % [-] compressibility factor (node based)
cc2n = ZZn.*RR.*Tn;   % [m2/s2] squared speed of sound (node based)
                    % it is a quantity directly available from the EoS
%     rho_n(:,1)=p_0(:,1)./cc2n;
rho_n(:,1) = Den.*MM'; % [kg/m3] actual density of the gas (node based)
                    % taken directly from the EoS solver. Otherwise
                    % it could be calculated with the expression at
                    % the previous line actual means that it is the density of
					% the gas at the pressure and the temperature at the specific
					% node or pipe recall of the EoS to calculate the nodal based
					% quantities at STANDARD CONDITION of T and p so to be able
					% to calculate the quantities listed few lines below
iFlag = 0;
p_is = p_std*ones(size(p_0(:,ii)))/1e3;
T_is = T_std*ones(size(Tn));
[Zm0(:,1), Den1] = PropertiesGERG(iFlag, p_is, T_is, x_ni, dimn, gerg_ni);


ZZN(:,1) = ZZn;
RHO_S(:,1)= rho_n(:,1).*(p_std./p_0(:,1)).*(Tn./T_std).*(ZZN(:,1)./Zm0(:,1));
RD(:,1)  = RHO_S(:,1)./rho_air_std;
HHV(:,1) = MHV(:,1).*RHO_S(:,1);
WI(:,1)  = HHV(:,1)./sqrt(RD(:,1));

% PIPELINE BASED CALCULATION
% the same EoS function is fed by pipeline based quantities such
% as:pm (pressure averaged over the pipeline) or
% Aplus'*reshape(CC_gas(:,:,1),dimn,21) -this is a way to
% transfer the outlet nodal quantities to the related pipe
iFlag = 0;
x_bi = Aplus'*reshape(CC_gas(:,:,1), dimn, 21);
gerg_bi = UtilitiesGERG(x_bi, dimb);
[Zm, Den]  = PropertiesGERG(iFlag, pm(:,1)/1e3, Tb, x_bi, dimb, gerg_bi);

ZZb = Zm;			 % [-] compressibility factor (node based)% [-] compressibility factor (pipe based)
RRb = Aplus'*RR;     % [%J/kg/K]  gas constant (pipe based)-the nodal value of the outlet node has been trasferred to the previous pipe
cc2b = ZZb.*RRb.*Tb; % [m2/s2] squared speed of sound (pipe based)

LP(:,ii) = dx.*pm(:,1).*AA./cc2b; %[kg] Linepack - amount of gas accumulated in dx

% alfa_G = 0.8;
alfa_P = 0;
alphaC = 1.5;

%% Time Cycle / for the transient model
for ii = 2:dimt

    ii %time step counter

    k2  = 1;
    res = 1;  % residue of the linearized fluid-dynamic problem
    %res2=1;  % another different calculation method for convergence gas quality
    res3= 1;  % residue/error for the convergence loop on gas quality

    CC_gasM_k(:,:,k2) = CC_gasM(:,:,ii-1);  %( tra 0-100)
    CC_gas_k(:,:,k2)  = CC_gas(:,:,ii-1);   %( tra 0-1)

    CC_gasM_ext_k(:,:,k2) = CC_gasM_ext(:,:,ii-1);
	% CC_gasM_k(:,:,k2) = CC_gasM(:,:,1);
	% CC_gas_k(:,:,k2) = CC_gas(:,:,1);

    iter_max2=1;
    %% CYCLE FOR THE QUALITY TRACKING CONVERGENCE
    while res3>0.5*1e-2 && iter_max2<500
		iter_max2 = iter_max2+1;

		k = 1;  % linearization loop counter

        CC_gasM_ext_k(:,:,k2) = CC_gasM_ext(:,:,ii-1);

		p_n = p_0(:,ii-1);
		G_n = G_0(:,ii-1);

		p_k(:,k) = p_0(:,ii-1);
		G_k(:,k) = G_0(:,ii-1);
		L_k(:,k) = L_0(:,ii-1);

		x_n = reshape(CC_gas_k(:,:,k2),dimn,21);
		x_b = Aplus' * x_n;

		gerg_b = UtilitiesGERG(x_b, dimb);
		gerg_n = UtilitiesGERG(x_n, dimn);

		mass_frac = reshape(CC_gasM_k(:,:,k2), dimn, nComponents);
		MM  = MolarMassGERG2(mass_frac);%kg/kmol %% mappa concentrazioni di tutti i nodi all'istante iniziale
		RR  = UGC./MM';
		RRb = Aplus'*RR;

		% MHV_ext(:,k2)=MASS_HV(reshape(CC_gas_k_ext(:,:,k2),dimn,21)/100)';
		MHV_ext(:,k2) = MASS_HV(x_n/100)';
		MHV_ext(find(isnan(MHV_ext(:,k2))),k2) = 1;
		G_ext_t(:,ii) = H_ext_t(:,ii)./(MHV_ext(:,k2)*1e3);

		% xx_old_k=xx_old;
		% % vel_b_k.ccb=vel_b.ccb;
		% vel_b_k=vel_b;
		% iter_max=1;

		% xx_old_k=xx_old;
		% for xxx=1:INPUT_b_inner(end,1)
		% 	vel_b_k(xxx).ccb=vel_b(xxx).ccb;
		% end
		iter_max = 1;

		%% LINEARIZED FLUID-DYNAMIC PROBLEM CYCLE
		%for any timestep, a solution for the linearized problem should be found.
		while res>toll && iter_max<500

			%this are underrelaxation coefficients to help the convergence, if
			%needed
			alfa_G = 0;
			alfa_P = 0;

			% if iter_max>10
			%    alfa_G=1;
			%    alfa_P=0.8;
			%    iter_max
			% else
			%    alfa_G=0;
			% end

			iter_max = iter_max+1; %fluid-dynamic iteration counter

			%% MOMENTUM EQUATION: each branch section CV - dimb

			%average pressure update
			p_in_k  = Aplus' * p_k(:,k);
			p_out_k = Aminus'* p_k(:,k);
			pm(:,k) = (2.0/3.0) * averagePressure(p_in_k, p_out_k);

			%update of all the PIPELINE BASED properties with Equation of State
            [Zm, Den] = PropertiesGERG(iFlag, pm(:,k)/1e3, Tb, x_n, dimb, gerg_b); %WK: why x_n instead of x_b?

			ZZb  = Zm;
			cc2b = ZZb.*RRb.*Tb;
			%rho(:,ii) = pm(:,k)./cc2b;
			rho(:,ii)  = Den.*(Aplus'*MM');       % [kg/m3] actual density of the gas (pipeline based)
			vel(:,ii)  = G_k(:,k)./AA./rho(:,ii); % [m/s] velocity of the gas within pipes.

			s = 2 * 9.81 * delH./cc2b;
			Aminus_s  = Aminus.*repmat(exp(s),1,dimn)';
			Aminus_s1 = Aminus.*repmat(exp(s/2),1,dimn)';
			ADP = (-Aminus_s1 + Aplus)';

			Ri = zeros(dimb,1);                 % Inertia Resistance - initialization
			Ri = 2*dxe.*pm(:,k)./(AA*dt)./(abs(ADP)*p_k(:,k)); % Inertia Resistance

			MOLEFRAC = x_b;
			[lambda Re viscosity] =  friction(Tb,epsi,G_k(:,k),DD,MOLEFRAC);
			Rf = 16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p_k(:,k)); % Fluid-dynamikc Resistance
			rr_k = (2*Rf.*(abs(G_k(:,k)))+Ri); % composite resistance linearized problem (R)
			R_k = sparse(diag(rr_k));          % transformed into sparse diagonal matrix

			%% CONTINUITY EQUATION: each node CV - dimn
			%update of all the NODE BASED properties with Equation of State
            [Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n, dimn, gerg_n);

			ZZn  = Zm;
			cc2n = ZZn.*RR.*Tn;
			Vol  = pi/8*abs(A)*(DD.^2.*dx);
			phi = Vol./(cc2n.*dt);
			PHI = sparse(diag(phi));

			Ap = A;

			II = sparse(eye(size(PHI)));

			% %   % NON PIPELINES ELEMENTS CONTROL MODE - p_out fixed
			% %   ADP(NP,in_np)=diag(zeros(size(NP))); %p_in
			% %   ADP(NP,out_np)=diag(ones(size(NP))); %p_out
			% %   R_k(NP,:)=zeros(size(R_k(NP,:)));
			% % %     R_k(54,54)=+1e7;
			% if p_0(in_np,ii-1)/p_0(out_np,ii-1)>=1.0%max(p_0(out_np:end,ii-1))<=min(p_0(1:in_np,ii-1))
			%     % NON PIPELINES ELEMENTS CONTROL MODE - bypass
			%     ADP(NP,in_np)=-diag(ones(size(NP))); %p_in
			%     ADP(NP,out_np)=diag(ones(size(NP))); %p_out
			%     R_k(NP,NP)=zeros(size(R_k(NP,NP)));
			%     D=zeros(size(NP));
			% else
			%     % NON PIPELINES ELEMENTS CONTROL MODE - OFF
			%
			%     ADP(NP,in_np)=-zeros(ones(size(NP))); %p_in
			%     ADP(NP,out_np)=zeros(ones(size(NP))); %p_out
			%     R_k(NP,NP)=-ones(size(R_k(NP,NP)));
			%     D=zeros(size(NP));
			% end

			%% MATRIX PROBLEM CONSTRUCTION:

			MATRIX_P_k=[ADP,-R_k,zeros(dimb,dimn)];  %momentum equation

			MATRIX_k=[PHI,Ap,II;                      %cont. eq
						MATRIX_P_k;                   %mom. eq
						zeros(dimn,dimn+dimb),II];    %boundary condition (pressure and gas flow)

			% KNOWN TERM
			% TN_P=PHI*p_n-II./2*L_0(:,ii-1);
			TN_P = PHI*p_n;
			TN_M_k = (-Rf.*abs(G_k(:,k)).*G_k(:,k) - Ri.*G_n);
			% TN_M_k(NP)=P_in(out_np);%(p_in_t(ii)-p_0(54,ii-1));%-G_0(54,ii-1);
			% TN_M_k(NP)=D;
			TN_L=G_ext_t(:,ii);

			%  %no backflow al gate
			%  if p_in_t(ii)>p_k(2,k)
			if vel(1,ii)>0
				MATRIX_k(dimn+dimb+1,1)=1;
				MATRIX_k(dimn+dimb+1,dimn+dimb+1)=0;
				TN_L(1)=p_in_t(ii);
			else
				MATRIX_k(dimn+dimb+1,1)=0;
				MATRIX_k(dimn+dimb+1,dimn+dimb+1)=1;
				TN_L(1)=0;
			end

			TN_k=[TN_P; TN_M_k; TN_L];        % full vector of KNOWN TERMs composition

			%% LINEAR SOLUTION OF THE LINEARIZED FLUID-DYNAMIC PROBLEM

			XXX_k = MATRIX_k\TN_k;

			%%
			p_k(:,k+1) = XXX_k(1:dimn);           % extraction results for nodal pressures
			G_k(:,k+1) = XXX_k(dimn+1:dimn+dimb); % extraction results for pipeline mass flows
			L_k(:,k+1) = XXX_k(dimn+dimb+1:end);  % extraction results for nodal mass flows exchanged with outside.

			%check! to avoid infinite resistances, if a pipe has 0 as mass flow it
			%is apporximated to 10^8
			if any(G_k(PIPE,k+1)==0)==1
				G_k(find(G_k(PIPE,k+1)==0),k+1)=1e-8;
			end

			k = k+1; %linearization index update

			%% update of all the quantities and residual calculation
		    % MOMENTUM EQUATION
		    % average pressure update
		    pm(:,ii) = (2.0/3.0)*averagePressure(Aplus'*p_k(:,k), Aminus'*p_k(:,k));
		    % update of all the PIPELINE BASED properties with Equation of State
			[Zm, Den] = PropertiesGERG(iFlag, pm(:,ii)/1e3, Tb, x_b, dimb, gerg_b);

			ZZb  = Zm;
			cc2b = ZZb.*RRb.*Tb;
			%rho(:,ii)= pm(:,ii)./cc2b;
			rho(:,ii) = Den.*(Aplus'*MM');
			vel(:,ii) = G_k(:,k)./AA./rho(:,ii); %m/s

			s = 2*9.81*delH./cc2b;
			Aminus_s  = Aminus.*repmat(exp(s),1,dimn)';
			Aminus_s1 = Aminus.*repmat(exp(s/2),1,dimn)';
			ADP = (-Aminus_s1+Aplus)';

			Ri = zeros(dimb,1);
			Ri = 2*dxe.*pm(:,ii)./(AA*dt)./(abs(ADP)*p_k(:,k));

			[lambda Re viscosity] = friction(Tb,epsi,G_k(:,k),DD,x_b);
			Rf = 16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p_k(:,k));
			rr_k = (2*Rf.*(abs(G_k(:,k)))+Ri);
			R_k  = sparse(diag(rr_k));

			RES_P(k) = norm(ADP(:,:)*p_k(:,k) - (Rf(:).*(abs(G_k(:,k)).*G_k(:,k))+Ri(:).*(G_k(:,k)-G_n(:))));
			%  norm(ADP(PIPE,:)*p_k(:,k) - (Rf(PIPE).*(abs(G_k(PIPE,k)).*G_k(PIPE,k))+Ri(PIPE).*(G_k(PIPE,k)-G_n(PIPE))));

			% CONTINUITY EQUATION
			% update of all the NODE BASED properties with Equation of State
            [Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n, dimn, gerg_n);

			ZZn = Zm;
			cc2n = ZZn.*RR.*Tn;
			Vol = pi/8*abs(A)*(DD.^2.*dx);
			phi = Vol./(cc2n.*dt);
			PHI = sparse(diag(phi));

			II=sparse(eye(size(PHI)));

			%RES_M(k) = norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II./2*(L_0(:,ii-1)+L_k(:,k)));
			 RES_M(k) = norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II*L_k(:,k));
			%RES_M(k) = norm(PHI*p_k(:,k) + A(:,PIPE)*G_k(PIPE,k) - PHI*p_n + II*L_k(:,k));

			RES_C(k) = 0;

			res = max([RES_P(k),RES_M(k),RES_C(k)])

			p_kf = p_k(:,k);
			G_kf = G_k(:,k);
			p_k(:,k) = alfa_P*p_k(:,k-1)+(1-alfa_P)*p_k(:,k);
			G_k(:,k) = alfa_G*G_k(:,k-1)+(1-alfa_G)*G_k(:,k);

			Residui(ii).Fluid(k2).RES_P(k) = RES_P(k);
			Residui(ii).Fluid(k2).RES_M(k) = RES_M(k);

			Res_P1(:,k) = ((p_k(:,k)-p_k(:,k-1))./p_k(:,k-1))*100;
			Res_G1(:,k) = ((G_k(:,k)-G_k(:,k-1))./G_k(:,k-1))*100;

		end

		if iter_max>=500
			disp('Warning: No Convergence FLUID')
			NO_FL_CONV=[NO_FL_CONV,ii];
			%break
		end

        p_0(:,ii) = p_kf;
        G_0(:,ii) = G_kf;
        L_0(:,ii) = L_k(:,k);
        LP(:,ii)  = dx.*pm(:,ii).*AA./cc2b; %kg

		% figure
		% plot(G_0)
		% figure
		% plot(p_0)

		CC_gasM_k(:,:,k2+1) = CC_gasM_k(:,:,k2);
		CC_gasM_ext_k(:,:,k2+1) = CC_gasM_ext_k(:,:,k2);

		E_conc = (CC_gasM_k(:,:,k2+1)-CC_gasM_k(:,:,k2))./(CC_gasM_k(:,:,k2))*100;
		% E_conc(find(isnan(E_conc)))=zeros(size(find(isnan(E_conc))));
		% E_conc(find(isinf(E_conc)))=zeros(size(find(isinf(E_conc))));

		ERR_C(k2) = max(max(abs(E_conc)));
		res3 = ERR_C(k2)/100

		% Residui(ii).RES_C(1,k2)=RES_C(k2);
		Residui(ii).RES_C(2,k2) = res3;
		%
		k2 = k2 + 1;
		mass_frac = reshape(CC_gasM_k(:,:,k2),dimn,nComponents);
		MM = MolarMassGERG2(mass_frac);%kg/kmol %% mappa concentrazioni di tutti i nodi all'istante iniziale
		CC_gas_k(:,:,k2) = MASS2MOL_CONV3(CC_gasM_k(:,:,k2),MM)./100;
		CC_gasM(:,:,ii)  = CC_gasM_k(:,:,k2);
		CC_gasM(:,:,ii)  = CC_gasM_k(:,:,k2);
		RR = UGC./MM';

	end

	if iter_max2>=500
		disp('No Convergence COMP')
		NO_FL_CONC = [NO_FL_CONC,ii];
		%break
	end
	%
	% xx_old = xx_old_k;
	%
	% for xxx = 1:INPUT_b_inner(end,1)
	% 	vel_b(xxx).ccb = vel_b_k(xxx).ccb;
	% end

	CC_gasM(:,:,ii)     = CC_gasM_k(:,:,k2);     %( tra 0-100)
	CC_gasM_ext(:,:,ii) = CC_gasM_ext_k(:,:,k2); %( tra 0-100)

	MM = MolarMassGERG2(reshape(CC_gasM(:,:,ii),dimn,nComponents));%kg/kmol %% mappa concentrazioni di tutti i nodi all'istante iniziale
	RR = UGC./MM';
	MolMass(:,ii) = MM;

	CC_gas(:,:,ii)=CC_gas_k(:,:,k2);%MASS2MOL_CONV3(CC_gasM(:,:,ii),MM)./100;
	MHV(:,ii) = MASS_HV(reshape(CC_gasM(:,:,ii),[dimn 21])/100)';

    iFlag = 0;
	gerg_nn = UtilitiesGERG(x_n, dimn);
	[Zm, Den, gamma]=PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n,dimn,gerg_nn);

	ZZn  = Zm;
	cc2n = ZZn.*RR.*Tn;
	rho_n(:,ii) = Den.*MM';
	[Zm0(:,ii), Den1] = PropertiesGERG(iFlag, p_std*ones(size(p_0(:,ii)))/1e3, T_std*ones(size(Tn)), x_n,dimn, gerg_nn);

	ZZN(:,ii)  = ZZn;
	RHO_S(:,ii)= rho_n(:,ii).*(p_std./p_0(:,ii)).*(Tn./T_std).*(ZZN(:,ii)./Zm0(:,ii));%[kg/Sm3]
	RD(:,ii)  = RHO_S(:,ii)./rho_air_std;
	HHV(:,ii) = MHV(:,ii).*RHO_S(:,ii);%[MJ/Sm3]
	WI(:,ii)  = HHV(:,ii)./sqrt(RD(:,ii));

	CONC_H2(:,ii) = CC_gas(:,15,ii);
end

N_Courandt_ADV = vel(:,end).*dt./dx
N_Courandt_FLU = (sqrt(cc2b)).*dt./dx
time = cputime - time0

save('Trial2023_3.mat')


