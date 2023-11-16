function [p_0 G_0] = SteadyState(A, Aplus, Aminus, P_in, Tb, G_ext_t, AA, DD, dxe, delH,epsi,RR,MM,CC_gas_k,NP,PIPE,dimn,dimb,dimt)
	k  = 1;
	k2 = 1;
	p_k(:,k)= P_in(1)*ones(dimn,1);
	p_n = P_in(1)*ones(dimn,1);
	G_n = ones(dimn,1);
	G_k(:,k)= G_n(:);		%KCommit: mass flow \dot[m]
	res = 1;
	toll = 1e-4;
	ii = 1;
	Tn = (273.15 + 20)*ones(dimn,1)
	iFlag = 0;
	x_n = reshape(CC_gas_k(:,:,k2),dimn,21) ;
	x_b = Aplus' * x_n;

	gerg_b = UtilitiesGERG(x_b, dimb);
	gerg_n = UtilitiesGERG(x_n, dimb);
	RRb  = Aplus'*RR;

	dt = 1.0; % WK: Do it by me. To impose Ri = 0; Since it is used in computation of resistance and it is in the denominator. This should be done in a better way
	steady = true;

	%% LINEARIZED FLUID-DYNAMIC PROBLEM CYCLE
	%for any timestep, a solution for the linearized problem should be found.

	[rho(:,ii), vel(:,ii), pm, p_k, p_kf, G_k, G_kf, L_k, ResiduiS, res] = ...
									SIMPLE(toll, k, k2, ii, iFlag, ...
									A, AA, Aplus, Aminus, ...
									p_k, p_n, P_in(1), G_k, G_n, G_ext_t(:,ii),...
									x_n, x_b, gerg_b, gerg_n, RRb, Tb, Tn, MM, RR,...
									dimb, dimn, dxe, dxe, dt, delH, ...
									epsi, DD, PIPE, res, steady);

	%{
	iter_max = 1;
	while res>toll && iter_max<500

		alfa_G = 0;
		alfa_P = 0;
		iter_max = iter_max + 1; %fluid-dynamic iteration counter

		%% A. MOMENTUM EQUATION: each branch section CV - dimb
		% A.1. UPDATE: pm and all PIPELINE-BASED properties with Equation of State
		p_k_in  =  Aplus' * p_k(:,k);
		p_k_out =  Aminus'* p_k(:,k);
		pm(:,k) = (2.0/3.0)*averagePressure(p_k_in, p_k_out);

		[Zm, Den] = PropertiesGERG(iFlag, pm(:,k)/1e3, Tb, x_b, dimb, gerg_b); %WK why not x_b1?
		cc2b = speedSound(Zm,RRb,Tb);

		% rho(:,ii) = pm(:,k)./cc2b;
		rho(:,ii) = Den.*(Aplus'*MM');       % [kg/m3] actual density of the gas (pipeline based)
		vel(:,ii) = G_k(:,k)./AA./rho(:,ii); % [m/s] velocity of the gas within pipes.

		%-----------------------------------------------------------------------
		% A.2  ADP Computation
		ADP = matrixADP(Aminus, Aplus, cc2b, delH, dimn);
		%-----------------------------------------------------------------------
		% A.3. COMPUTATION OF Rf/ Ri (Ri = 0 due to steady state assumption)
		% pm =0 to induce Ri = 0
		Omega = Resistance(dxe, 1.0, zeros(dimb,1), p_k(:,k), G_k(:,k), AA, ADP, x_b, Tb, epsi, DD, cc2b);

		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		%% B. CONTINUITY EQUATION: each node CV - dimn
		% B.1 update of all the NODE BASED properties with Equation of State
		% gerg_1 = UtilitiesGERG(x_1, dimn); %WK is there any difference using  reshape(CC_gas_k(:,:,k2),dimn,21)  or  reshape(CC_gas_k(:,:,1),dimn,21). It is mixed in the whole steady state,but k2 =1 and never changes.
											 %WK: what is the need to recompute here again gerg_1 (which is equal to gerg_n)
		[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n, dimn, gerg_n); %WK why not x_1? instead of  x_n?
		cc2n = speedSound(Zm, RR, Tn);	%KCommit: Sound speed squared. Why not RRb?

		%WK:  DIFF
		phi = zeros(size(DD));
		PHI = sparse(diag(phi));
		II = sparse(eye(size(PHI)));

		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		%% 			C. MATRIX PROBLEM CONSTRUCTION:
		%-----------------------------------------------------------------------
		%% C.3 SOLVE: LINEAR SOLUTION OF THE LINEARIZED FLUID-DYNAMIC PROBLEM
		[p_k(:,k+1), G_k(:,k+1), L_k(:,k+1)] = matrix_system(dimn, dimb,...
											A, ADP, Omega, II, PHI, vel(1,ii),...
											p_n, P_in(1), G_k(:,k), G_n, G_ext_t);
		%WK: DIFF => P_in(1) change by p_int_t
		%-------------------------------------------------------------------

		%check! to avoid infinite resistances, if a pipe has 0 as mass flow it
		%is apporximated to 10^8
		if any(G_k(PIPE,k+1)==0)==1
		   G_k(find(G_k(PIPE,k+1)==0),k+1) = 1e-8;
		end

		%--kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
		k = k+1; %linearization index update
		%--kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

		%% update of all the quantities and residual calculation
		%---------------------------------------------------------------------------
		%---------------------------------------------------------------------------
		% MOMENTUM EQUATION
		%---------------------------------------------------------------------------
		% 1.UPDATE: pm and PIPELINE BASED properties with Equation of State
		p_in  =  Aplus' * p_k(:,k);
		p_out =  Aminus'* p_k(:,k);
		pm(:,ii) = (2.0/3.0)*averagePressure(p_in, p_out);

		[Zm, Den] = PropertiesGERG(iFlag, pm(:,ii)/1e3, Tb, x_b, dimb, gerg_b);
		cc2b = speedSound(Zm, RRb, Tb);

		rho(:,ii) = Den.*(Aplus'*MM');
		vel(:,ii) = G_k(:,k)./AA./rho(:,ii); %m/s
		%---------------------------------------------------------------------------
		ADP = matrixADP(Aminus, Aplus, cc2b, delH, dimn);
		%---------------------------------------------------------------------------
		% 2. RESISTANCES: %WK: DIFF =>  Ri = 0 due to steady state assumption)
		% pm = 0, dt =1.0  to induce Ri = 0
		Omega = Resistance(dxe, 1.0, zeros(dimb,1), p_k(:,k), G_k(:,k), AA, ADP, x_b, Tb, epsi, DD, cc2b);
		%---------------------------------------------------------------------------
		%---------------------------------------------------------------------------
		% 					CONTINUITY EQUATION
		%---------------------------------------------------------------------------
		% UPDATE: NODE-BASED properties with Equation of State
		[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n, dimn, gerg_n);
		cc2n = speedSound(Zm, RR,Tn);
		%WK: DIFF phi not updated
		PHI = sparse(diag(phi));	%WK:  here phi is not updated. So it is just zeros?
		II = sparse(eye(size(PHI)));
		%---------------------------------------------------------------------------
		%---------------------------------------------------------------------------
		% 						RESIDUALS
		%---------------------------------------------------------------------------
		RES_P(k) = norm(ADP(:,:)*p_k(:,k) ...
					- (Omega.Rf(:).*(abs(G_k(:,k)).*G_k(:,k))
						+	Omega.Ri(:).*(G_k(:,k) - G_n(:))));
		RES_M(k) = norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II*L_k(:,k)); %WK: inconsistent dimension PHI*p_n' chnaged to => PHI*p_n
		RES_C(k)=0;

		Residui(ii).Fluid(k2).RES_P(k) = RES_P(k);
		Residui(ii).Fluid(k2).RES_M(k) = RES_M(k);

		res = max([RES_P(k),RES_M(k),RES_C(k)])
		%---------------------------------------------------------------------------
		% Update p* and G*
		p_kf = p_k(:,k);
		G_kf = G_k(:,k);
		p_k(:,k) = alfa_P*p_k(:,k-1)+(1-alfa_P)*p_k(:,k);
		G_k(:,k) = alfa_G*G_k(:,k-1)+(1-alfa_G)*G_k(:,k);
		%---------------------------------------------------------------------------
		% WK: Res_P1/Res_G1 are not call anywhere else. Then I commented them
		Res_P1(:,k) = ((p_k(:,k)-p_k(:,k-1))./p_k(:,k-1))*100;
		Res_G1(:,k) = ((G_k(:,k)-G_k(:,k-1))./G_k(:,k-1))*100;

	end
%}
	G_0 = G_kf;
	p_0 = p_kf;

end

% WK questions:
% What is k2 intended for? is only defined at the beginning
% What are Residui, Res_P1, Res_G1 intended for here? not used and not set as outputs

