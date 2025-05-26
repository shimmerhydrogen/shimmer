function[P,Z,dPdDsave,alpha_r] = PressureGERG(T,D,x,dimn,Tr,Dr)
%Sub PressureGERG(T, D, x, P, Z)

%Calculate pressure as a function of temperature and density.  The derivative d(P)/d(D) is also calculated
%for use in the iterative DensityGERG subroutine (and is only returned as a common variable).

%Inputs:
%    T - temperature (K)
%    D - density (mol/l)
%    x() - composition (mole fraction)
%          Do not send mole percents or mass fractions in the x() array, otherwise the output will be incorrect.
%          The sum of the compositions in the x() array must be equal to one.
% 
%Outputs:
%      P - pressure (kPa)
%      Z - compressibility factor
%  dPdDsave - d(P)/d(D) [kPa/(mol/l)]
RGERG = 8.314472;
  
  alpha_r = AlpharGERG(0, T, D, x,dimn,Tr,Dr); %dim(ar)=(3,3)
  Z = 1 + alpha_r(:,1, 2);
  P = D .* RGERG .* T .* Z;
  dPdDsave = RGERG .* T .* (1 + 2 .* alpha_r(:,1, 2) + alpha_r(:,1, 3));
  
end 