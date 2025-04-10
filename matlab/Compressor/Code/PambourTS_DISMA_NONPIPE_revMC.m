function [gasNet_Res, pos, BC_change_rec] = PambourTS_DISMA_NONPIPE_revMC(gasNet,Gsnam,P_set,Gguess,Pguess,dt)
disp("===================================================");
disp('                  TIME SOLVER                     ')
disp("===================================================");
time0=cputime;
toll=1e-4; % 1e-6 iniziale
%dt=1;
dimt=length(Gsnam(1,:));
MAX_ITER=1500;
%% sparse incidence matrix construction
nodei=gasNet.Edges.EndNodes(:,1); 
nodeo=gasNet.Edges.EndNodes(:,2);

dimb=length(nodei);
dimn=max(max([nodei,nodeo]));
BR=[1:1:dimb]';

Asp_p=sparse(nodei,BR,ones(1,dimb),dimn,dimb);
Asp_m=sparse(nodeo,BR,[-ones(1,dimb)],dimn,dimb);
Asp=Asp_p+Asp_m;

A=Asp;
Aplus=(A==1).*1;
Aminus=(A==-1).*1;

% f1=figure;
% p1=plot(gasNet,'Xdata',gasNet.Nodes.coordinates_XY(:,1),'Ydata',gasNet.Nodes.coordinates_XY(:,2));

%% STD Conditions
p_std = 1.01325*1e5; %Pa
T_std = 273.15 + 15;%K
% isothermal conditions  
T=20; %°C
Tn=(273.15 + T)*ones(dimn,1);% K nodal temperature
Tb=(273.15 + T)*ones(dimb,1);% K pipes temperature

UGC=8.3144598*1e3;                  %J/kmol/K  % Universal Gas Constant
% MM=MolarMassGERG1(reshape(MOLEFRAC_INPUT,[size(INLET,1),nComponents])); %kg/kmol 
MM=16.*ones(dimn,1);% kg/kmol
% calculation of the molar mass of gas(es) at the inlet.
RR=UGC./MM; %J/kg/K  % gas constant of the gas at the inlet

%% network technological parameters

LL=gasNet.Edges.Length;               %m  %length of the pipelines 
dx=LL;
dxe=dx;                               % vectoor size: #pipelines
HH=zeros(dimn,1);                     %m  % altitude of each node - size: #node
delH=-A'*HH;                          %m  % h_out-h_in / delta height %m - size: #pipeline  
DD=gasNet.Edges.Diameter;             %m  % Diameters of the pipelines
epsi=gasNet.Edges.Epsi;               %m  % absolute roughness
AA=pi.*DD.^2/4;                       %m2 % pipeline cross section

eta=1e-5;                             %kg/m/s % dynamic viscosity
RD1=0.6;                              %[-]    % Relative Density
rho_air_std=1.225;                    %kg/Sm3 @p=1 bar T=15 °C - 288K
%rho_air_std=1.292;                   %kg/Nm3 @p=1 bar T=0 °C - 273K
rho_g_std=RD1*rho_air_std;            %kg/Sm3  % density of natural gas @standard cond

%% boundary conditions of mass flow rates and pressures

% INLET=find(G_ext<0);  %nodi di ingresso
% OUTLET=find(G_ext>0); %nodi di uscita
% INNER=find(G_ext==0); %nodi interni - non scambiano con l'esterno

G_ext=zeros(size(gasNet.Nodes.Nodes_ID)); %kg/s
G_ext=zeros(size(Gsnam)); %kg/s
G_ext=Gsnam;
% OUTLETS=find(gasNet.Nodes.deg==1);
% % for ii=1:
% G_ext(OUTLETS,:)=Gsnam(OUTLETS,:);

P_set=P_set*10^5;
Pguess = Pguess*10^5;

% KNOWN=find(gasNet.Nodes.isSource==1);

users=find(G_ext>0);

G_n(:)=Gguess;%G_0(:,1);    % matrix (bxt) of pipeline mass flow rates for each timestep
p_n(:)=Pguess;%p_0(:,1);    % matrix (nxt) of nodal pressures for each timestep
L_0(:,1)=Gsnam(:,1);        % matrix (nxt) of nodal inlet(-)/outlet(+) mass flow rates for each timestep
                            % quasi-fixed term: it is mass flow rated fixed
                            % by the THERMAL REQUEST thus depends on the
                            % gas quality at each node

%% LINEARIZED FLUID-DYNAMIC PROBLEM CYCLE
%for any timestep, a solution for the linearized problem should be found.
ii=1;
res(1)=1;
iter_max=0;
k=1;
%initialization of the linearization loop (k-referred)
p_k(:,k)=Pguess;
G_k(:,k)=Gguess;
L_k(:,k)=G_ext(:,ii+1);
system_size = 2 * size(Pguess,2) +  size(Gguess,2); 
XXX_time = zeros(system_size, dimt); 
XXX_time(:,1) = [Pguess'; Gguess'; G_ext(:,ii+1)]; 
% G_n=zeros(size(G_k));
% p_n=zeros(size(p_k));
G_n=G_k(:,k);
p_n=p_k(:,k);
RRb=Aplus'*RR;

pm(:,k)=2/3*(Aplus'*p_k(:,k).^2+Aminus'*p_k(:,k).^2+(Aplus'*p_k(:,k)).*(Aminus'*p_k(:,k)))./((Aplus'*p_k(:,k))+(Aminus'*p_k(:,k)));

ZZb=Papay(pm(:,k),Tb);
cc2b=ZZb.*RRb.*Tb;

LP(:,ii)=dx.*pm(:,1).*AA./cc2b; %[kg] Linepack - amount of gas accumulated in dx

II=sparse(eye(dimn));
III=II; % Matrix ext mass flow set BC
OObn=zeros(dimb,dimn);
OOnn=zeros(dimn,dimn);% Matrix pressure set BC
OOnb=zeros(dimn,dimb);

% boundary conditions (time-varying)
entry_p=find(gasNet.Nodes.Type==1);
entry_l=find(gasNet.Nodes.Type==2);
exit=find(gasNet.Nodes.Type==3);

%boudary condition equations
%(ENTRY) pressure set points
OOnn(entry_p,entry_p)=eye(length(entry_p));
III(entry_p,entry_p)=zeros(length(entry_p));

%boundary condition values:
TN_L=G_ext; %default condition (there is already also the one of injection)
TN_L(entry_p,:)=P_set(entry_p,:);


%==================================================================
% START TIME
%==================================================================
for ii=2:dimt
    disp("========================================================");
    disp(strcat(" Solving at time ..........",int2str(ii-1)));
    disp("========================================================");
    res(1)=1;
    iter_max=0;
    k=1;
    FLAG_BC=1;

    %==================================================================
    % START ITER HARD CONSTRAINTS 
    %==================================================================

    while FLAG_BC~=0 


        disp("****************************************************************");
        disp(strcat(" Iteration CONSTRAINTS it ...............",int2str(iter_max)," ... at time ", int2str(ii-1)));
        disp("****************************************************************");


        res(1)=1;
        iter_max=0;
        k=1;
        XXX_iter = [];
        TN_iter = [];
        %==================================================================
        % START FLUID SOLVER
        %==================================================================
    
        while res(k)>toll && iter_max<MAX_ITER
        
            %underrelaxation coefficients to help the convergence, if needed
            alfa_G=0;
            alfa_P=0;
            
            
            iter_max = iter_max+1; %fluid-dynamic iteration counter
            
            %% MOMENTUM EQUATION: each branch section CV - dimb
            
            %average pressure update
            pm(:,k)=2/3*(Aplus'*p_k(:,k).^2+Aminus'*p_k(:,k).^2+(Aplus'*p_k(:,k)).*(Aminus'*p_k(:,k)))./((Aplus'*p_k(:,k))+(Aminus'*p_k(:,k)));
            % %update of all the PIPELINE BASED properties with Equation of State  
            
            % ZZb=Zm;
            ZZb=Papay(pm(:,k),Tb);
            cc2b=ZZb.*RRb.*Tb;
            rho(:,ii)=pm(:,k)./cc2b;
            % rho(:,ii)=Den.*(Aplus'*MM');       % [kg/m3] actual density of the gas (pipeline based)
            vel(:,ii)=G_k(:,k)./AA./rho(:,ii); % [m/s] velocity of the gas within pipes.
            
            s=2*9.81*delH./cc2b;
            Aminus_s=Aminus.*repmat(exp(s),1,dimn)';
            Aminus_s1=Aminus.*repmat(exp(s/2),1,dimn)';
            ADP=(-Aminus_s1+Aplus)';
            
            %     Ri=zeros(dimb,1);                 % Inertia Resistance - initialization
            %     ZERO IN CASE OF STEADY STATE 
            Ri=2*dxe.*pm(:,k)./(AA*dt)./(abs(ADP)*p_k(:,k)); % Inertia Resistance
            
            % [lambda Reyn viscosity]=Frictionfactoraverage(Tb,epsi,G_k(:,k),DD,Aplus'*reshape(CC_gas_k(:,:,k2),dimn,21));
            lambda=Frictionfactoraverage(Tb,epsi,G_k(:,k),DD);
            Rf=16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p_k(:,k)); % Fluid-dynamikc Resistance
            rr_k=(2*Rf.*(abs(G_k(:,k)))+Ri); % composite resistance linearized problem (R)
            R_k=sparse(diag(rr_k));          % transformed into sparse diagonal matrix
            
            %% CONTINUITY EQUATION: each node CV - dimn

            % ZZn=Zm;
            ZZn=Papay(p_k(:,k),Tn);
            cc2n=ZZn.*RR.*Tn;
            Vol=pi/8*abs(A)*(DD.^2.*dx);
            % Vol=zeros(size(cc2n));
            phi=Vol./(cc2n.*dt);
            PHI=sparse(diag(phi));
            
            Ap=A;
            
            
            
            %% MATRIX PROBLEM CONSTRUCTION:
            
            % KNOWN TERM
            %TN_P=PHI*p_n-II./2*L_0(:,ii-1);
            TN_P=PHI*p_n;
            TN_M_k=(-Rf.*abs(G_k(:,k)).*G_k(:,k) -Ri.*G_n);
            
            %NON PIPELINE
            if length(find(gasNet.Edges.COMP>0))==0
            else
                COMP=find(gasNet.Edges.COMP>0);
                for ccc=1:length(COMP)
                
                if gasNet.Edges.COMP_ctrl.RegType{ccc}==1
                    %PowerControl
                    eta_comp=0.9;
                    kappa=1.4;
                    ck=(kappa-1)/kappa;
                    Ki=ZZb(COMP(ccc))*Tb(COMP(ccc))*RRb(COMP(ccc))/eta_comp;
                    pin=p_k(gasNet.Edges.EndNodes(COMP(ccc),1),k);
                    pout=p_k(gasNet.Edges.EndNodes(COMP(ccc),2),k);
                    beta=pout/pin;
                    % beta=(gasNet.Edges.COMP_ctrl.PWR/G_k(COMP(ccc),k)*ck/Ki+1)^(1/ck)
                    % beta=1.1;
                    c1=-Ki*G_k(COMP(ccc),k)/pin*beta^ck;
                    c2=Ki*G_k(COMP(ccc),k)/pout*beta^ck;
                    c3=-Ki/ck*(beta^ck-1);
                    d=gasNet.Edges.COMP_ctrl.PWR{1}; %kW
                elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==2
                    %Pout
                    c1=0; c2=1; c3=0; d=gasNet.Edges.COMP_ctrl.p_out{1};
                elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==3
                    %Pin
                    c1=1; c2=0; c3=0; d=gasNet.Edges.COMP_ctrl.p_in{1};
                elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==4
                    %beta
                    c1=-gasNet.Edges.COMP_ctrl.beta{1}; c2=1; c3=0; d=0;
                elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==5
                    %Q
                    c1=0; c2=0; c3=-1; d=gasNet.Edges.COMP_ctrl.Q{1};
                else
                    disp('Error COMP control mode')
                end
                    ADP(COMP(ccc),gasNet.Edges.EndNodes(COMP(ccc),1))=c1;
                    ADP(COMP(ccc),gasNet.Edges.EndNodes(COMP(ccc),2))=c2;
                    R_k(COMP(ccc),COMP(ccc))=c3;
                    TN_M_k(COMP(ccc),1)=d;
                end
            end
               
            MATRIX_P_k = [ADP,-R_k,OObn];   %momentum equation 
            MATRIX_k = [PHI,Ap,II;          %cont. eq
            MATRIX_P_k;                     %mom. eq
            OOnn,OOnb,III];       %boundary condition (pressure and gas flow)
            
            
            % %boundary condition values            
            TN_k=[TN_P; TN_M_k; TN_L(:,ii)];        % full vector of KNOWN TERMs composition  
            
            %% LINEAR SOLUTION OF THE LINEARIZED FLUID-DYNAMIC PROBLEM
            
            XXX_k = MATRIX_k\TN_k;
            XXX_iter(:,iter_max) = XXX_k;
            TN_iter(:,iter_max) = TN_k;
            %%     
            p_k(:,k+1)=XXX_k(1:dimn);           % extraction results for nodal pressures
            G_k(:,k+1)=XXX_k(dimn+1:dimn+dimb); % extraction results for pipeline mass flows
            L_k(:,k+1)=XXX_k(dimn+dimb+1:end);  % extraction results for nodal mass flows exchanged with outside.

            k=k+1; %linearization index update
            
            %% update of all the quantities and residual calculation
            % MOMENTUM EQUATION
            %average pressure update
            pm(:,k)=2/3*(Aplus'*p_k(:,k).^2+Aminus'*p_k(:,k).^2+(Aplus'*p_k(:,k)).*(Aminus'*p_k(:,k)))./((Aplus'*p_k(:,k))+(Aminus'*p_k(:,k)));
            % %update of all the PIPELINE BASED properties with Equation of State  
            
            % ZZb=Zm;
            ZZb=Papay(pm(:,k),Tb);
            cc2b=ZZb.*RRb.*Tb;
            rho(:,ii)=pm(:,k)./cc2b;
            % rho(:,ii)=Den.*(Aplus'*MM');       % [kg/m3] actual density of the gas (pipeline based)
            vel(:,ii)=G_k(:,k)./AA./rho(:,ii); % [m/s] velocity of the gas within pipes.
            
            s=2*9.81*delH./cc2b;
            Aminus_s=Aminus.*repmat(exp(s),1,dimn)';
            Aminus_s1=Aminus.*repmat(exp(s/2),1,dimn)';
            ADP=(-Aminus_s1+Aplus)';
            
            %     Ri=zeros(dimb,1);
            %     ZERO IN CASE OF STEADY STATE 
            Ri=2*dxe.*pm(:,k)./(AA*dt)./(abs(ADP)*p_k(:,k)); % Inertia Resistance
            
            % [lambda Reyn viscosity]=Frictionfactoraverage(Tb,epsi,G_k(:,k),DD,Aplus'*reshape(CC_gas_k(:,:,k2),dimn,21));
            lambda=Frictionfactoraverage(Tb,epsi,G_k(:,k),DD);
            Rf=16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p_k(:,k));
            rr_k=(2*Rf.*(abs(G_k(:,k)))+Ri);
            R_k=sparse(diag(rr_k));
            
            TN_M_k=(-Rf.*abs(G_k(:,k)).*G_k(:,k) -Ri.*G_n);
            
            %NON PIPELINE
            if length(find(gasNet.Edges.COMP>0))==0
            else
                COMP=find(gasNet.Edges.COMP>0);
                for ccc=1:length(COMP)
                
                    if gasNet.Edges.COMP_ctrl.RegType{ccc}==1
                        %PowerControl
                        eta_comp = 0.9;
                        kappa    = 1.4;
                        ck  = (kappa-1)/kappa;
                        Ki  = ZZb(COMP(ccc))*Tb(COMP(ccc))*RRb(COMP(ccc))/eta_comp;
                        pin = p_k(gasNet.Edges.EndNodes(COMP(ccc),1),k);
                        pout= p_k(gasNet.Edges.EndNodes(COMP(ccc),2),k);
                        beta= pout/pin;
                        % beta=(gasNet.Edges.COMP_ctrl.PWR/G_k(COMP(ccc),k)*ck/Ki+1)^(1/ck);
                        c1 = -Ki*G_k(COMP(ccc),k)/pin*beta^ck; 
                        c2 = Ki*G_k(COMP(ccc),k)/pout*beta^ck; 
                        c3 = -Ki/ck*(beta^ck-1); 
                        d  = gasNet.Edges.COMP_ctrl.PWR{1}; %W
                        
                        control_name = "pwd";
                    elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==2
                        %Pout
                        c1=0; c2=1; c3=0; d=gasNet.Edges.COMP_ctrl.p_out{1};
                        control_name = "pout";
                    elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==3
                        %Pin
                        c1=1; c2=0; c3=0; d=gasNet.Edges.COMP_ctrl.p_in{1};
                        control_name = "pin";
                    elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==4
                        %beta
                        c1= -gasNet.Edges.COMP_ctrl.beta{1}; c2=1; c3=0; d=0;
                        control_name = "beta";
                    elseif gasNet.Edges.COMP_ctrl.RegType{ccc}==5
                        %Q
                        c1=0; c2=0; c3= -1; d=gasNet.Edges.COMP_ctrl.Q{1};
                        control_name = "flux";
                    else
            
                    end
                    ADP(COMP(ccc),gasNet.Edges.EndNodes(COMP(ccc),1))=c1;
                    ADP(COMP(ccc),gasNet.Edges.EndNodes(COMP(ccc),2))=c2;
                    R_k(COMP(ccc),COMP(ccc))=c3;
                    TN_M_k(COMP(ccc),1)=d;
                end
            end
            
            RES_P0(k) = norm([ADP,-R_k,zeros(dimb,dimn)]*XXX_k-TN_M_k);

            % CONTINUITY EQUATION

            
            % ZZn=Zm;
            ZZn=Papay(p_k(:,k),Tn);
            Vol=pi/8*abs(A)*(DD.^2.*dx);
            %Vol=zeros(size(cc2n));
            phi=Vol./(cc2n.*dt);
            PHI=sparse(diag(phi));  
            

            %RES_M(k)= norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II./2*(L_0(:,ii-1)+L_k(:,k)));
            RES_M(k)= norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II*L_k(:,k));
            %RES_M(k)= norm(PHI*p_k(:,k) + A(:,PIPE)*G_k(PIPE,k) - PHI*p_n + II*L_k(:,k));
            
            RES_P(k)=norm(p_k(:,k)-p_k(:,k-1));
            RES_M(k)=norm(G_k(:,k)-G_k(:,k-1));
            RES_C(k)=0;   
            
            res(k)=max([RES_P(k),RES_M(k),RES_C(k)]);
            
            p_kf=p_k(:,k);
            G_kf=G_k(:,k);
            p_k(:,k)=alfa_P*p_k(:,k-1)+(1-alfa_P)*p_k(:,k);
            G_k(:,k)=alfa_G*G_k(:,k-1)+(1-alfa_G)*G_k(:,k);
            
            Res_P1(:,k)=((p_k(:,k)-p_k(:,k-1))./p_k(:,k-1))*100;
            Res_G1(:,k)=((G_k(:,k)-G_k(:,k-1))./G_k(:,k-1))*100;
            
        end %% CHIUSURA CICLO LINEARIZATION

        %==================================================================
        % END FLUID SOLVER
        %==================================================================

        if iter_max>=MAX_ITER
            disp('No Convergence FLUID')
            % NO_FL_CONV=[NO_FL_CONV,ii];
            %break
            pos=1;
            WARN=1;
        else
            pos=0;
            WARN=0;
        end
        
        % check on boundary conditions
        if any(L_k(entry_p,k) > 0 )
            disp("*******************************************");
            disp("WARNING: HARD CONSTRAINT violated for REMI");
            disp("*******************************************");
            FLAG_P(ii) = 1;
            entry_p_change_l = entry_p(find(L_k(entry_p,k)> 0));
            OOnn(entry_p_change_l,entry_p_change_l)= zeros(length(entry_p_change_l));
            III(entry_p_change_l,entry_p_change_l) = eye(length(entry_p_change_l));
            TN_L(entry_p_change_l,ii:end) = zeros(length(entry_p_change_l),dimt-ii+1);
            % gasNet.Nodes.Type(aa)==2;
            % ii = ii-1; %% WARN!!! da cambiare 
            BC_change_rec(ii).P = entry_p_change_l;
            % entry_p_change_l=zeros(size(entry_p_change_l));
            entry_p_change_l = [];
        else
            disp("*******************************************");
            disp(" HARD CONSTRAINT in  REMI............ [PASS]");
            disp("*******************************************");

            FLAG_P(ii) =0;
        end
        
        
        if any(L_k(entry_l,k) > 0) || any(p_kf(entry_l) > 75*1e5) % I assumed 65 as maximum pressure
            disp("*******************************************");
            disp("WARNING: HARD CONSTRAINT violated for INJ");
            disp("*******************************************");
            % p_k(aa,k+1) = gasNet.Nodes.Pset_bc(aa);
            % gasNet.Nodes.Type(aa)==1;
            % ii = ii-1;
            FLAG_L(ii) =1;
            entry_l_change_p=entry_l(find(p_kf(entry_l) > 75*1e5));
            OOnn(entry_l_change_p,entry_l_change_p)=eye(length(entry_l_change_p));
            III(entry_l_change_p,entry_l_change_p)=zeros(length(entry_l_change_p));
            TN_L(entry_l_change_p,ii:end)=75*1e5;%P_set(entry_l_change_p,dimt-ii);%zeros(length(entry_p_change_l),dimt-ii);
            % gasNet.Nodes.Type(aa)==2;
            % ii = ii-1; %% WARN!!! da cambiare 
            BC_change_rec(ii).L=entry_l_change_p;
            % entry_l_change_p=zeros(size(entry_l_change_p));
            entry_l_change_p=[];
        else
            disp("*******************************************");
            disp(" HARD CONSTRAINT in  INJ............. [PASS]");
            disp("*******************************************");

            FLAG_L(ii)=0;
        end
        
        FLAG_BC=any([FLAG_P(ii),FLAG_L(ii)]>0);
        %% per saltare il check delle condizioni (commentare le righe sopra e
        %% decommentare queste due sotto)
        % FLAG_BC=0;
        % BC_change_rec(ii)=0;
        
        for kkk=1:length(exit)
            aa=exit(kkk);
            if L_k(aa,k) < 0
            disp("L < 0 in exit station, check input data")
            end
        end
        
        
        
        
    end %% CHIUSURA CICLO FLAG
   
    %==================================================================
    % END ITER HARD CONSTRAINTS 
    %==================================================================

    
    p_0(:,ii)=p_kf;
    G_0(:,ii)= G_kf;
    L_0(:,ii)= L_k(:,k);
    LP(:,ii)=dx.*pm(:,k).*AA./cc2b; %kg
    
    G_n=G_0(:,ii);
    p_n=p_0(:,ii);


    XXX_time(:,ii) = XXX_k; 
end  %% CHIUSURA CICLO TIME

%==================================================================
% END TIME
%==================================================================
if WARN==1
    gasNet_Res=gasNet;
    gasNet_Res.Edges.FLOWRATES= G_kf.*NaN ; %kg/s
    gasNet_Res.Nodes.PRESSURES=p_kf.*NaN; %Pa
    gasNet_Res.Nodes.G_EXE= L_k(:,k).*NaN; %kg/s
    % gasNet_Res.Nodes.MolComposition= MOLFRAC(:,:,end).*NaN; %p.u
else
    gasNet_Res=gasNet;
    gasNet_Res.Edges.FLOWRATES= G_0 ; %kg/s
    gasNet_Res.Nodes.PRESSURES=p_0; %Pa
    gasNet_Res.Nodes.G_EXE= L_0; %kg/s
    % gasNet_Res.Nodes.MolComposition= MOLFRAC(:,:,end); %p.u.
end


% ABS_ERR_PRESS=abserrorp(end);%Pa
% rel_error_P_perc=max(abs(Pguess_rec(:,end-1)-Pguess_rec(:,end))./Pguess_rec(:,end))*100;%
% 
% ABS_ERR_Gpipe=abserrorg(end);%kg/s
% rel_error_Gpipe_perc=max(abs(Gguess_rec(:,end-1)-Gguess_rec(:,end))./Gguess_rec(:,end))*100;%
% 

save(strcat("../Comparison/var_time_DENERG_", control_name, ".dat"), "XXX_time", "-ascii")

time=cputime-time0;
end


% figure
% plot(p_k(:,end))
% hold on
% plot(p_k(:,end-1))
% title('Pressure')
% 
% p_k(146,end)/10^5
% 
% figure
% plot(G_k(:,end))
% hold on
% plot(G_k(:,end-1))
% title('G branch')
% figure
% plot(L_k(:,end))
% hold on
% plot(L_k(:,end-1))
% title('G ext')
