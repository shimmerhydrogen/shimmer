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
			%-------------------------------------------------------------------
			%% 			A. MOMENTUM EQUATION: each branch section CV - dimb
			%-------------------------------------------------------------------
			% A.1 UPDATE:PIPELINE-BASED properties with Equation of State
			pm(:,k) = (2.0/3.0) * averagePressure(Aplus' * p_k(:,k), Aminus'* p_k(:,k));
            [Zm, Den] = PropertiesGERG(iFlag, pm(:,k)/1e3, Tb, x_n, dimb, gerg_b); %WK: why x_n instead of x_b?
			cc2b = speedSound(Zm,RRb, Tb);
			rho  = Den.*(Aplus'*MM');       % [kg/m3] actual density of the gas (pipeline based)
			vel  = G_k(:,k)./AA./rho; % [m/s] velocity of the gas within pipes.
			%-------------------------------------------------------------------
			ADP = matrixADP(Aminus, Aplus, cc2b, delH, dimn);
			%-------------------------------------------------------------------
			if(steady)
				% pm =0 to induce Ri = 0. dt = 1.0
				Omega = Resistance(dxe, dt, zeros(dimb,1), p_k(:,k),G_k(:,k), AA, ADP, x_b, Tb, epsi, DD, cc2b);
			else
				Omega = Resistance(dxe, dt, pm(:,k), p_k(:,k),G_k(:,k), AA, ADP, x_b, Tb, epsi, DD, cc2b);
			end
			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
			%% 			B. CONTINUITY EQUATION: each node CV - dimn
			%-------------------------------------------------------------------

			%-------------------------------------------------------------------
			if(~steady)
				% B.1 UPDATE: NODE-BASED properties with Equation of State
				[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n, dimn, gerg_n);
				cc2n = speedSound(Zm,RR, Tn);
				[PHI, II] = matrixPHI(dx, dt, cc2n, A, DD);
			else
				phi = zeros(size(DD));
				PHI = sparse(diag(phi));
				II = sparse(eye(size(PHI)));
			end
			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
			%% 			C. FLUID-DYNAMICS MATRIX PROBLEM:
			%-------------------------------------------------------------------
			[p_k(:,k+1), G_k(:,k+1), L_k(:,k+1)] = matrix_system(dimn, dimb, ...
										A, ADP, Omega, II, PHI, vel(1),...
										p_n, p_in_t, G_k(:,k), G_n, G_ext_t);
			%-------------------------------------------------------------------

			%check! to avoid infinite resistances, if a pipe has 0 as mass flow it
			%is apporximated to 10^8
			if any(G_k(PIPE,k+1)==0)==1
				G_k(find(G_k(PIPE,k+1)==0),k+1)=1e-8;
			end

			%--kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
			k = k+1; %linearization index update
			%--kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

			%% update of all the quantities and residual calculation
			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
		    % 						B. MOMENTUM EQUATION
			%-------------------------------------------------------------------
		    % UPDATE: PIPELINE-BASED properties with Equation of State
		    pm(:,ii) = (2.0/3.0)*averagePressure(Aplus'*p_k(:,k), Aminus'*p_k(:,k)); %WK:  why pm here is using ii instead of k. previously was save with k. Then, it is being overwritten every iteration. If I dont put pm as an output it will be out of bounds. I dont understand how this is organized.
			[Zm, Den] = PropertiesGERG(iFlag, pm(:,ii)/1e3, Tb, x_b, dimb, gerg_b);
			cc2b = speedSound(Zm, RRb, Tb);
			%rho(:,ii)= pm(:,ii)./cc2b;
			rho = Den.*(Aplus'*MM');
			vel = G_k(:,k)./AA./rho; %m/s
			%-------------------------------------------------------------------
			ADP = matrixADP(Aminus, Aplus, cc2b, delH, dimn);
			%-------------------------------------------------------------------
			if(steady)
				%WKcommit: pm = 0 and dt = 1.0 to induce Ri = 0.
				Omega = Resistance(dxe, dt, zeros(dimb,1), p_k(:,k), G_k(:,k), AA, ADP, x_b, Tb, epsi, DD, cc2b)
			else
				Omega = Resistance(dxe, dt, pm(:,ii), p_k(:,k), G_k(:,k), AA, ADP, x_b, Tb, epsi, DD, cc2b)
			end

			%-------------------------------------------------------------------
			% 					CONTINUITY EQUATION
			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
			if(~steady)
				% UPDATE: NODE-BASED properties with Equation of State
				[Zm, Den, gamma] = PropertiesGERG(iFlag, p_k(:,k)/1e3, Tn, x_n, dimn, gerg_n);
				cc2n = speedSound(Zm, RR, Tn);
				[PHI, II] = matrixPHI(dx, dt, cc2n, A, DD);
			else
				PHI = sparse(diag(phi));	%WK:  here phi is not updated. So it is just zeros?
				II = sparse(eye(size(PHI)));
			end
			%-------------------------------------------------------------------
			%-------------------------------------------------------------------
			% 						Find Residuals
			%-------------------------------------------------------------------
			RES_P = norm(ADP(:,:)*p_k(:,k)...
					- (Omega.Rf(:).*(abs(G_k(:,k)) .* G_k(:,k))...
							+ Omega.Ri(:).*(G_k(:,k) - G_n(:))));
			RES_M = norm(PHI*p_k(:,k) + A*G_k(:,k) - PHI*p_n + II*L_k(:,k));
			RES_C = 0;

			res = max([RES_P,RES_M,RES_C])

			Residui.Fluid(k2).RES_P(k) = RES_P;
			Residui.Fluid(k2).RES_M(k) = RES_M;

			%-------------------------------------------------------------------
			% Update p* and G*
			p_kf = p_k(:,k);
			G_kf = G_k(:,k);
			p_k(:,k) = alfa_P * p_k(:,k-1) + (1.0 - alfa_P) * p_k(:,k);
			G_k(:,k) = alfa_G * G_k(:,k-1) + (1.0 - alfa_G) * G_k(:,k);
			%-------------------------------------------------------------------
			% WK: Res_P1/Res_G1 are not call anywhere else. Then I commented them
			%Res_P1(:,k) = ((p_k(:,k) - p_k(:,k-1))./p_k(:,k-1))*100;
			%Res_G1(:,k) = ((G_k(:,k) - G_k(:,k-1))./G_k(:,k-1))*100;
		end

		if iter_max>=500
			disp('WARNING: No convergence FLUID')
			NO_FL_CONV=[NO_FL_CONV,ii];
			%break
		end

end
