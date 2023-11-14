function [P1, Z, D, gamma]=GERG(iFlag,P, T, x, dimn, gerg)%(T, D, x)
%Sub PropertiesGERG(ByVal T As Double, ByVal D As Double, x() As Double, ByRef P As Double, ByRef Z As Double, ByRef dPdD As Double, ByRef d2PdD2 As Double, ByRef d2PdTD As Double, ByRef dPdT As Double, ByRef U As Double, ByRef H As Double, ByRef S As Double, ByRef Cv As Double, ByRef Cp As Double, ByRef W As Double, ByRef G As Double, ByRef JT As Double, ByRef Kappa As Double, Optional ByRef A As Double)
%Sub PropertiesGERG(T, D, x, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, Optional A)

%Calculate thermodynamic properties as a function of temperature and density.
%If the density is not known, call subroutine DensityGERG first with the known values of pressure and temperature.
%Many of the formulas below do not appear in Part 2 of AGA 8, but rather in Part 1, which uses a dimensional Helmholtz equation with more direct formulas for quick calculation.
%
%Inputs:
%     T - temperature (K)
%     D - density (mol/l)
%   x() - composition (mole fraction)
%
%Outputs:
%     P - pressure (kPa)
%     Z - compressibility factor
%  dPdD - first derivative of pressure with respect to density [kPa/(mol/l)]
%d2PdD2 - second derivative of pressure with respect to density [kPa/(mol/l)^2]
%d2PdTD - second derivative of pressure with respect to temperature and density [kPa/(mol/l)/K]
%  dPdT - first derivative of pressure with respect to temperature (kPa/K)
%     U - internal energy (J/mol)
%     H - enthalpy (J/mol)
%     S - entropy (J/mol-K)
%    Cv - isochoric heat capacity (J/mol-K)
%    Cp - isobaric heat capacity (J/mol-K)
%     W - speed of sound (m/s)
%     G - Gibbs energy (J/mol)
%    JT - Joule-Thomson coefficient (K/kPa)
% Kappa - Isentropic Exponent
%     A - Helmholtz energy (J/mol)

%   MaxFlds=length(x);
% global xold;
%   for i = 1 : MaxFlds
%     xold(i) = 0;
%   end

%Dim a0(2) As Double, ar(3, 3) As Double, Mm As Double, R As Double, RT As Double
RGERG = 8.314472;
%Calculate molar mass
Mm = MolarMassGERG(x);%Call MolarMassGERG(x, Mm)

[D, ierr, herr, alpha_r]=DensityGERG(iFlag, T, P, x, dimn,gerg);

%Calculate the ideal gas Helmholtz energy, and its first and second derivatives with respect to temperature.
a0= Alpha0GERG(T, D, x);%Call Alpha0GERG(T, D, x, a0)

%Calculate the real gas Helmholtz energy, and its derivatives with respect to temperature and/or density.
%   alpha_r = AlpharGERG(1,T,D,x);%Call AlpharGERG(1, T, D, x, ar)

R = RGERG;
RT = R * T;
%%
Z = 1 + alpha_r(:,1, 2);
P1 = D .* RT .* Z;
dPdD = RT .* (1 + 2 * alpha_r(:,1, 2) + alpha_r(:,1, 3));
dPdT = D .* R .* (1 + alpha_r(:,1, 2) - alpha_r(:,2, 2));
%   d2PdTD = R .* (1 + 2 * alpha_r(:,1, 2) + alpha_r(:,1, 3) - 2 * alpha_r(:,2, 2) - alpha_r(:,2, 3));
%   A = RT .* (a0(:,1) + alpha_r(:,1, 1));
%   G = RT .* (1 + alpha_r(:,1, 2) + a0(:,1) + alpha_r(:,1, 1));
%   U = RT .* (a0(:,2) + alpha_r(:,2, 1));
%   H = RT .* (1 + alpha_r(:,1, 2) + a0(:,2) + alpha_r(:,2, 1));
%   S = R .* (a0(:,2) + alpha_r(:,2, 1) - a0(:,1) - alpha_r(:,1, 1));
Cv = -R .* (a0(:,3) + alpha_r(:,3, 1));
if D > 0
	p = zeros(size(D));
	%  d2PdD2=zeros(size(D));
	%  JT=zeros(size(D));
	%
    Cp(find(D>0)) = Cv(find(D>0)) + T(find(D>0)) .* (dPdT(find(D>0)) ./ D(find(D>0))) .^ 2 ./ dPdD(find(D>0));
	% d2PdD2(find(D>0)) = RT(find(D>0)) .* (2 * alpha_r((find(D>0)),1, 2) + 4 * alpha_r((find(D>0)),1, 3) + alpha_r((find(D>0)),1, 4)) ./ D(find(D>0));
	% JT(find(D>0)) = (T(find(D>0)) ./ D(find(D>0)) .* dPdT(find(D>0)) ./ dPdD(find(D>0)) - 1) ./ Cp(find(D>0)) ./ D(find(D>0));  %=(dB/dT*T-B)/Cp for an ideal gas, but dB/dT is not known
else
	Cp(find(D<=0)) = Cv(find(D<=0)) + R;
	% d2PdD2(find(D<=0)) = 0;
	% JT(find(D<=0)) = 1E+20;
end
gamma = Cp./Cv;

%   W = 1000 .* Cp ./ Cv .* dPdD ./ Mm;
%   if any(W < 0)
%       W(find(W<0)) = 0;
%   end
%   W = sqrt(W);
%   Kappa = W .^ 2 .* Mm ./ (RT .* 1000 .* Z);
end

