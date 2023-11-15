function [p, G, L] = matrix_system(dimn, dimb, A, ADP, Omega, II, PHI, vel, p_n, p_in, G_k, G_n, G_ext)
	Ap = A;

	MATRIX_P_k = [ADP,- Omega.R, zeros(dimb,dimn)];  %momentum equation

	MATRIX_k = [PHI, Ap, II;                  %cont. eq
				MATRIX_P_k;                   %mom. eq
				zeros(dimn,dimn+dimb),II];    %boundary condition (pressure and gas flow)
	%  %no backflow al gate
	%  if p_in_t(ii)>p_k(2,k)
	if vel > 0
		MATRIX_k(dimn+dimb+1,1)=1;
		MATRIX_k(dimn+dimb+1,dimn+dimb+1)=0;
	else
		MATRIX_k(dimn+dimb+1,1)=0;
		MATRIX_k(dimn+dimb+1,dimn+dimb+1)=1;
	end
	%-------------------------------------------------------------------
	% C.2 Load vector
	TN_P = PHI * p_n;
	TN_M_k = (- Omega.Rf .* abs(G_k) .* G_k - Omega.Ri .* G_n);
	TN_L = G_ext;

	%  %no backflow al gate
	%  if p_in_t(ii)>p_k(2,k)
	if vel > 0
		TN_L(1) = p_in;
	else
		TN_L(1) = 0;
	end
	TN_k = [TN_P; TN_M_k; TN_L];        % full vector of KNOWN TERMs composition

	%-------------------------------------------------------------------
	%% C.3 SOLVE: LINEAR SOLUTION OF THE LINEARIZED FLUID-DYNAMIC PROBLEM
	sol = MATRIX_k\TN_k;

	p = sol(1:dimn);           % extraction results for nodal pressures
	G = sol(dimn+1:dimn+dimb); % extraction results for pipeline mass flows
	L = sol(dimn+dimb+1:end);  % extraction results for nodal mass flows exchanged with outside.
end
