function [PHI, II] = matrixPHI(dx, dt, cc2n, A, DD)
	Vol  = pi/8*abs(A)*(DD.^2.*dx);
	phi = Vol./(cc2n.*dt);
	PHI = sparse(diag(phi));
	II = sparse(eye(size(PHI)));
end
