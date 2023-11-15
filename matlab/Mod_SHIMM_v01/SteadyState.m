function [p_0 G_0] = SteadyState(A, Aplus, Aminus, P_in, Tb, G_ext_t, AA, DD, dxe, delH,epsi,RR,MM,CC_gas_k,NP,PIPE,dimn,dimb,dimt)
	k  = 1;
	k2 = 1;
	p_k(:,k)= P_in(1)*ones(dimn,1);
	p_n(:)  = P_in(1)*ones(dimn,1);
	G_n(:)  = ones(dimn,1);
	G_k(:,k)= G_n(:);		%KCommit: mass flow \dot[m]
	res = 1;
	iter_max = 1;
	toll = 1e-4;
	ii = 1;
	Tn = (273.15 + 20)*ones(dimn,1)
	while res>toll && iter_max<500

		alfa_G = 0;
		alfa_P = 0;

		% if iter_max>10
		%    alfa_G=1;
		%    alfa_P=0.8;
		%    iter_max
		% else
		%    alfa_G=0;
		% end
		%
		% if iter_max>1500
		% alfa_G=0.1;
		% alfa_P=0.99;
		% iter_max
		% else
		%  alfa_G=0;
		%  alfa_P=0;
		% end

		iter_max = iter_max + 1; %fluid-dynamic iteration counter

		x_1= reshape(CC_gas_k(:,:,1),dimn,21);
		x_k2 = reshape(CC_gas_k(:,:,k2),dimn,21) ;

		%% A. MOMENTUM EQUATION: each branch section CV - dimb

		% A.1. UPDATE VARIABLES

		% KCommit: Eq pag 28 pm = (p_in^2 + p_in*p_out + p_out^2) /(p_in+p_out)
		p_k_in  =  Aplus' * p_k(:,k);
		p_k_out =  Aminus'* p_k(:,k);
		pm(:,k) = (2.0/3.0)*averagePressure(p_k_in, p_k_out);

		% update of all the PIPELINE BASED properties with Equation of State
		iFlag = 0;
		x_b = Aplus' * x_1;
		gerg_b = UtilitiesGREG(x_b, dimb);
		[Zm, Den] = PropertiesGERG(iFlag, pm(:,k)/1e3, Tb, Aplus'*x_k2, dimb,gerg_b); %WK why not x_b? instead of  Aplus'*x_k2
		%-----------------------------------------------------------------------
		% A.2  ADP Computation
		RRb = Aplus'*RR;
		ZZb = Zm;
		cc2b = ZZb.*RRb.*Tb;
		% rho(:,ii) = pm(:,k)./cc2b;
		rho(:,ii) = Den.*(Aplus'*MM');       % [kg/m3] actual density of the gas (pipeline based)
		vel(:,ii) = G_k(:,k)./AA./rho(:,ii); % [m/s] velocity of the gas within pipes.

		s = 2*9.81*delH./cc2b;
		Aminus_s  = Aminus.*repmat(exp(s),1,dimn)';
		Aminus_s1 = Aminus.*repmat(exp(s/2),1,dimn)';
		ADP = (-Aminus_s1+Aplus)';			% KCommit:  Ag't pag 27-28
		%-----------------------------------------------------------------------
		% A.3. COMPUTATION OF Rf/ Ri (Ri = 0 due to steady state assumption)
		Ri = zeros(dimb,1);                 % Inertia Resistance - initialization
		% Ri=2*dxe.*pm(:,k)./(AA*dt)./(abs(ADP)*p_k(:,k)); % Inertia Resistance

		[lambda Re viscosity] = friction(Tb,epsi,G_k(:,k),DD,Aplus'*x_k2);
		Rf = 16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p_k(:,k)); % Fluid-dynamic Resistance
		rr_k = (2*Rf.*(abs(G_k(:,k)))+Ri); % composite resistance linearized problem (R)
		R_k  = sparse(diag(rr_k));         % transformed into sparse diagonal matrix
		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		%% B. CONTINUITY EQUATION: each node CV - dimn
		% B.1 update of all the NODE BASED properties with Equation of State
		gerg = UtilitiesGREG(x, dimn);
		[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_k2, dimn, gerg);  %WK why not x? instead of  x_k2? ...check what x is

		ZZn = Zm;
		cc2n = ZZn.*RR.*Tn;					%KCommit: Sound speed squared
		Vol = pi/8*abs(A)*(DD.^2.*dxe);
		% phi = Vol./(cc2n.*dt);
		phi = zeros(size(Vol));
		PHI = sparse(diag(phi));

		Ap = A;

		II = sparse(eye(size(PHI)));

		% %     % NON PIPELINES ELEMENTS CONTROL MODE - p_out fixed
		% %     ADP(NP,in_np)=diag(zeros(size(NP))); %p_in
		% %     ADP(NP,out_np)=diag(ones(size(NP))); %p_out
		% %     R_k(NP,:)=zeros(size(R_k(NP,:)));
		% % %     R_k(54,54)=+1e7;
		% if p_0(in_np,ii-1)/p_0(out_np,ii-1)>=1.0%max(p_0(out_np:end,ii-1))<=min(p_0(1:in_np,ii-1))
		%     % NON PIPELINES ELEMENTS CONTROL MODE - bypass
		%     ADP(NP,in_np)  =-diag(ones(size(NP))); %p_in
		%     ADP(NP,out_np) = diag(ones(size(NP))); %p_out
		%     R_k(NP,NP)=zeros(size(R_k(NP,NP)));
		%     D = zeros(size(NP));
		% else
		%     % NON PIPELINES ELEMENTS CONTROL MODE - OFF
		%
		%     ADP(NP,in_np)  =-zeros(ones(size(NP))); %p_in
		%     ADP(NP,out_np) = zeros(ones(size(NP))); %p_out
		%     R_k(NP,NP) = -ones(size(R_k(NP,NP)));
		%     D = zeros(size(NP));
		% end

		%% MATRIX PROBLEM CONSTRUCTION:

		MATRIX_P_k = [ADP,-R_k, zeros(dimb,dimn)];    %momentum equation
		MATRIX_k   = [PHI,Ap,II;                      %cont. eq
					  MATRIX_P_k;                     %mom. eq
					  zeros(dimn,dimn+dimb),II];      %BCond (pressure and gas flow)

		% KNOWN TERM
		% TN_P = PHI*p_n-II./2*L_0(:,ii-1);
		TN_P = PHI*p_n';
		TN_M_k = (-Rf.*abs(G_k(:,k)).*G_k(:,k) -Ri.*G_n');
		% TN_M_k(NP) = P_in(out_np);%(p_in_t(ii)-p_0(54,ii-1));%-G_0(54,ii-1);
		% TN_M_k(NP) = D;
		TN_L = G_ext_t;

		% %no backflow al gate
		%  if p_in_t(ii)>p_k(2,k)
		if vel(1,ii)>0
			MATRIX_k(dimn+dimb+1,1) = 1;
			MATRIX_k(dimn+dimb+1,dimn+dimb+1) = 0;
			TN_L(1) = P_in(1);
		else
			MATRIX_k(dimn+dimb+1,1) = 0;
			MATRIX_k(dimn+dimb+1,dimn+dimb+1) = 1;
			TN_L(1)=0;
		end

		TN_k = [TN_P; TN_M_k; TN_L];        % full vector of KNOWN TERMs composition

		%% LINEAR SOLUTION OF THE LINEARIZED FLUID-DYNAMIC PROBLEM

		XXX_k = MATRIX_k\TN_k;

		p_k(:,k+1) = XXX_k(1:dimn);           % extraction results for nodal pressures
		G_k(:,k+1) = XXX_k(dimn+1:dimn+dimb); % extraction results for pipeline mass flows
		L_k(:,k+1) = XXX_k(dimn+dimb+1:end);  % extraction results for nodal mass flows exchanged with outside.

		%check! to avoid infinite resistances, if a pipe has 0 as mass flow it
		%is apporximated to 10^8
		if any(G_k(PIPE,k+1)==0)==1
		   G_k(find(G_k(PIPE,k+1)==0),k+1) = 1e-8;
		end

		k = k+1; %linearization index update

		%% update of all the quantities and residual calculation
		% MOMENTUM EQUATION
		%---------------------------------------------------------------------------
		% 1. Find ADP
		% average pressure update
		p_in  =  Aplus' * p_k(:,k);
		p_out =  Aminus'* p_k(:,k);
		pm(:,ii) = (2.0/3.0)*averagePressure(p_in, p_out);
		% update of all the PIPELINE BASED properties with Equation of State
		[Zm, Den] = PropertiesGERG(iFlag, pm(:,ii)/1e3, Tb, Aplus'*x_k2,dimb, gerg_b);

		ZZb = Zm;
		cc2b = ZZb.*RRb.*Tb;
		%rho(:,ii) = pm(:,ii)./cc2b;
		rho(:,ii) = Den.*(Aplus'*MM');
		vel(:,ii) = G_k(:,k)./AA./rho(:,ii); %m/s

		s = 2*9.81*delH./cc2b;
		Aminus_s  = Aminus.*repmat(exp(s),1,dimn)';
		Aminus_s1 = Aminus.*repmat(exp(s/2),1,dimn)';
		ADP = (-Aminus_s1+Aplus)';
		%---------------------------------------------------------------------------
		% 2. COMPUTATION OF Rf/ Ri (Ri = 0 due to steady state assumption)
		Ri = zeros(dimb,1);
		% Ri =2*dxe.*pm(:,ii)./(AA*dt)./(abs(ADP)*p_k(:,k));

		[lambda Re viscosity] = friction(Tb,epsi,G_k(:,k),DD,Aplus'*x_k2);
		Rf = 16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p_k(:,k));
		rr_k = (2*Rf.*(abs(G_k(:,k)))+Ri);
		R_k  = sparse(diag(rr_k));

		%---------------------------------------------------------------------------
		% 3. RESIDUAL
		RES_P(k) = norm(ADP(:,:)*p_k(:,k) - (Rf(:).*(abs(G_k(:,k)).*G_k(:,k))+Ri(:).*(G_k(:,k)-G_n(:))));
		% norm(ADP(PIPE,:)*p_k(:,k) - (Rf(PIPE).*(abs(G_k(PIPE,k)).*G_k(PIPE,k))+Ri(PIPE).*(G_k(PIPE,k)-G_n(PIPE))));

		% CONTINUITY EQUATION
		% update of all the NODE BASED properties with Equation of State
		[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_k2, dimn,gerg);

		ZZn = Zm;
		cc2n = ZZn.*RR.*Tn;
		Vol = pi/8*abs(A)*(DD.^2.*dxe);
		% phi = Vol./(cc2n.*dt);
		PHI = sparse(diag(phi));

		II = sparse(eye(size(PHI)));

		%RES_M(k) = norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II./2*(L_0(:,ii-1)+L_k(:,k)));
		RES_M(k) = norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n' + II*L_k(:,k));
		%RES_M(k) = norm(PHI*p_k(:,k) + A(:,PIPE)*G_k(PIPE,k) - PHI*p_n + II*L_k(:,k));

		RES_C(k)=0;

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
	G_0=G_kf;
	p_0=p_kf;

end

% WK questions:
% What is k2 intended for? is only defined at the beginning

