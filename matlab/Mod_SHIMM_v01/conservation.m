function [ADP, Omega, PHI, II, pm, rho, vel] = conservation(iFlag,...
									A, AA, Aplus, Aminus, ...
									p_k,  G_k, ...
									x_n, x_b, gerg_b, gerg_n, RR, Tb, Tn, MM, ...
									dimb, dimn, dxe, dx, dt, delH, ...
									epsi, DD, steady)
	%-------------------------------------------------------------------
	%-------------------------------------------------------------------
	%% 			A. MOMENTUM EQUATION: each branch section CV - dimb
	%-------------------------------------------------------------------
	% A.1 UPDATE: PIPELINE-BASED properties with Equation of State
	pm = (2.0/3.0) * averagePressure(Aplus' * p_k, Aminus'* p_k);
    [Zm, Den] = PropertiesGERG(iFlag, pm/1e3, Tb, x_b, dimb, gerg_b);
	cc2b = speedSound(Zm,Aplus'*RR, Tb);
	rho  = Den.*(Aplus'*MM');   % [kg/m3] actual density of the gas (pipeline based)
	vel  = G_k./AA./rho; 		% [m/s] velocity of the gas within pipes.
	%-------------------------------------------------------------------
	ADP = matrixADP(Aminus, Aplus, cc2b, delH, dimn);
	%-------------------------------------------------------------------
	if(steady)
		% pm =0 to induce Ri = 0. dt = 1.0
		Omega = Resistance(dxe, dt, zeros(dimb,1), p_k,G_k, AA, ADP, x_b, Tb, epsi, DD, cc2b);
	else
		Omega = Resistance(dxe, dt, pm, p_k,G_k, AA, ADP, x_b, Tb, epsi, DD, cc2b);
	end


	%-------------------------------------------------------------------
	%-------------------------------------------------------------------
	%% 			B. CONTINUITY EQUATION: each node CV - dimn
	%-------------------------------------------------------------------

	%-------------------------------------------------------------------
	if(~steady)
		% B.1 UPDATE: NODE-BASED properties with Equation of State
		[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k/1e3, Tn, x_n, dimn, gerg_n);
		cc2n = speedSound(Zm,RR, Tn);
		[PHI, II] = matrixPHI(dx, dt, cc2n, A, DD);
	else
		phi = zeros(size(DD));
		PHI = sparse(diag(phi));
		II = sparse(eye(size(PHI)));
	end
end
