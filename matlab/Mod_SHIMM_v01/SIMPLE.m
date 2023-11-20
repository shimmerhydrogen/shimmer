function  [rho, vel, pm,  p_k, p_kf, G_k, G_kf, L_k, Residui, res] = SIMPLE(tol, k, k2, ii, iFlag,...
									A, AA, Aplus, Aminus, ...
									p_k, p_n, p_in_t, G_k, G_n, G_ext_t,...
									x_n, x_b, gerg_b, gerg_n, RRb, Tb, Tn, MM, RR,...
									dimb, dimn, dxe, dx, dt, delH, ...
									epsi, DD, PIPE, res, steady)

	iter_max = 1;

	while res>tol && iter_max<500

		%Underrelaxation coefficients to help the convergence
		alfa_G = 0;
		alfa_P = 0;
		iter_max = iter_max+1; %fluid-dynamic iteration counter
		%-------------------------------------------------------------------
		[ADP, Omega, PHI, II, pm(:,k), rho, vel] = conservation(iFlag,...
								A, AA, Aplus, Aminus, ...
								p_k(:,k),  G_k(:,k), ...
								x_n, x_b, gerg_b, gerg_n, RR, Tb, Tn, MM, ...
								dimb, dimn, dxe, dx, dt, delH, ...
								epsi, DD, steady);

		%-------------------------------------------------------------------
		[p_k(:,k+1), G_k(:,k+1), L_k(:,k+1)] = matrix_system(dimn, dimb, ...
									A, ADP, Omega, II, PHI, vel(1),...
									p_n, p_in_t, G_k(:,k), G_n, G_ext_t);

		%check! to avoid infinite resistances, if a pipe has 0 as mass flow it
		%is apporximated to 10^8
		if any(G_k(PIPE,k+1)==0)==1
			G_k(find(G_k(PIPE,k+1)==0),k+1)=1e-8;
		end
		%-------------------------------------------------------------------
		%% update of all the quantities and residual calculation
		%-------------------------------------------------------------------
		[ADP, Omega, PHI, II, pm(:,ii), rho, vel] = conservation(iFlag,...
								A, AA, Aplus, Aminus, ...
								p_k(:,k+1),  G_k(:,k+1), ...
								x_n, x_b, gerg_b, gerg_n, RR, Tb, Tn, MM, ...
								dimb, dimn, dxe, dx, dt, delH, ...
								epsi, DD, steady)

		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% 						Find Residuals
		%-------------------------------------------------------------------
		RES_P = norm(ADP(:,:)* p_k(:,k+1)...
				- (Omega.Rf(:).*(abs(G_k(:,k+1)) .* G_k(:,k+1))...
						+ Omega.Ri(:).*(G_k(:,k+1) - G_n(:))));
		RES_M = norm(PHI*p_k(:,k+1) + A*G_k(:,k+1) - PHI*p_n + II*L_k(:,k+1));
		RES_C = 0;

		res = max([RES_P,RES_M,RES_C])

		Residui.Fluid(k2).RES_P(k+1) = RES_P;
		Residui.Fluid(k2).RES_M(k+1) = RES_M;

		%-------------------------------------------------------------------
		% Update p* and G*
		p_kf = p_k(:,k+1);
		G_kf = G_k(:,k+1);
		p_k(:,k+1) = alfa_P * p_k(:,k) + (1.0 - alfa_P) * p_k(:,k+1);
		G_k(:,k+1) = alfa_G * G_k(:,k) + (1.0 - alfa_G) * G_k(:,k+1);
		%-------------------------------------------------------------------
		% WK: Res_P1/Res_G1 are not call anywhere else. Then I commented them
		%Res_P1(:,k+1) = ((p_k(:,k+1) - p_k(:,k))./p_k(:,k))*100;
		%Res_G1(:,k+1) = ((G_k(:,k+1) - G_k(:,k))./G_k(:,k))*100;

		k = k + 1;

	end

	if iter_max>=500
		disp('WARNING: No convergence FLUID')
		NO_FL_CONV=[NO_FL_CONV,ii];
		%break
	end

end
