function[alpha_r] = AlpharGERG(iprop,T,D,x,dimn,Tr,Dr)
 
% Private Sub AlpharGERG(iprop, T, D, x, ar)
% 
% Calculate dimensionless residual Helmholtz energy and its derivatives with respect to tau and delta.
% 
% Inputs:
%  iprop - set to 1 to return all derivatives or 0 to return only pressure related properties [ar(0,1) and ar(0,2)]
%      T - temperature (K)
%      D - density (mol/l)
%    x() - composition (mole fraction)
% 
% Outputs:

%  ar(0,0) - residual Helmholtz energy (dimensionless, =a/RT)
%  ar(0,1) -     del*partial  (ar)/partial(del)
%  ar(0,2) -   del^2*partial^2(ar)/partial(del)^2
%  ar(0,3) -   del^3*partial^3(ar)/partial(del)^3
%  ar(1,0) -     tau*partial  (ar)/partial(tau)
%  ar(1,1) - tau*del*partial^2(ar)/partial(tau)/partial(del)
%  ar(2,0) -   tau^2*partial^2(ar)/partial(tau)^2

% %   Dim i As Integer, j As Integer, k As Integer, mn As Integer
% %   Dim Tr As Double, Dr As Double, del As Double, tau As Double
% %   Dim lntau As Double, ex As Double, ex2 As Double, ex3 As Double, cij0 As Double, eij0 As Double
% %   Dim delp(7) As Double, Expd(7) As Double, ndt As Double, ndtd As Double, ndtt As Double, xijf As Double

  %%%For i = 0 To 3: For j = 0 To 3: ar(i, j) = 0: Next: Next
  alpha_r = zeros(dimn,4,4);

%Set up del, tau, log(tau), and the first 7 calculations for del^i
%   [Tr,Dr] = ReducingParametersGERG(x);
  nn=length(D);
  del = D ./ Dr;
  tau = Tr ./ T;
  lntau = log(tau);
  delp = zeros(nn,7);
  Expd = zeros(nn,7);
  delp(:,1) = del;
  Expd(:,1) = exp(-delp(:,1));
  for i = 2 : 7
    delp(:,i) = delp(:,i - 1) .* del;
    Expd(:,i) = exp(-delp(:,i));
  end

%   Told = 0;
%   Trold2 = 0;
  
% global Told
% global Trold2
  Told = zeros(size(x,1),1);
  Trold2 = zeros(size(x,1),1);
global taup 
global taupijk

%If temperature has changed, calculate temperature dependent parts
  if norm(T - Told) > 0.00000001 || norm(Tr - Trold2) > 0.00000001
     [taup0, taup, taupijk] = tTermsGERG(tau, lntau, x); %Call tTermsGERG(tau, lntau, x)
  end
  Told = T;
  Trold2 = Tr;
  
  %taup0 = zeros(size(taup0));
  %taup = zeros(size(taup));
  %taupijk = zeros(size(taupijk));

%Calculate pure fluid contributions
NcGERG = size(x,2);
load('doik.mat','doik');
load('toik.mat','toik');
load('kpol.mat','kpol');
load('coik.mat','coik');
load('mNumb.mat','mNumb');
load('kexp.mat','kexp');
load('fij.mat','fij');
load('kpolij.mat','kpolij');
load('dijk.mat','dijk');
load('kexpij.mat','kexpij');
load('cijk.mat','cijk');
load('eijk.mat','eijk');
load('nijk.mat','nijk');
load('gijk.mat','gijk');
load('tijk.mat','tijk');

for nn = 1 : size(x,1)
  for i = 1 : NcGERG
    if x(nn,i) > 0
      for k = 1 : kpol(i)
        ndt(nn) = x(nn,i) * delp(nn,doik(i, k)) * taup(nn,i, k);
        ndtd(nn) = ndt(nn) * doik(i, k);
        alpha_r(nn,1, 2) = alpha_r(nn,1, 2) + ndtd(nn);
        alpha_r(nn,1, 3) = alpha_r(nn,1, 3) + ndtd(nn)* (doik(i, k) - 1);
        if iprop > 0
          ndtt(nn) = ndt(nn) * toik(i, k);
          alpha_r(nn,1, 1) = alpha_r(nn,1, 1) + ndt(nn);
          alpha_r(nn,2, 1) = alpha_r(nn,2, 1) + ndtt(nn);
          alpha_r(nn,3, 1) = alpha_r(nn,3, 1) + ndtt(nn) * (toik(i, k) - 1);
          alpha_r(nn,2, 2) = alpha_r(nn,2, 2) + ndtt(nn) * doik(i, k);
          alpha_r(nn,2, 3) = alpha_r(nn,2, 3) + ndtt(nn) * doik(i, k) * (doik(i, k) - 1);
          alpha_r(nn,1, 4) = alpha_r(nn,1, 4) + ndtd(nn) * (doik(i, k) - 1) * (doik(i, k) - 2);
        end
      end
      for k = (1 + kpol(i)):(kpol(i) + kexp(i))
        ndt(nn) = x(nn,i) * delp(nn,doik(i, k)) * taup(nn,i, k) * Expd(nn,coik(i, k));
        ex(nn) = coik(i, k) * delp(nn,coik(i, k));
        ex2(nn)= doik(i, k) - ex(nn);
        ex3(nn) = ex2(nn) * (ex2(nn) - 1);
        alpha_r(nn,1, 2) = alpha_r(nn,1, 2) + ndt(nn) * ex2(nn);
        alpha_r(nn,1, 3) = alpha_r(nn,1, 3) + ndt(nn) * (ex3(nn) - coik(i, k) * ex(nn));
        if iprop > 0 
          ndtt(nn) = ndt(nn) * toik(i, k);
          alpha_r(nn,1, 1) = alpha_r(nn,1, 1) + ndt(nn);
          alpha_r(nn,2, 1) = alpha_r(nn,2, 1) + ndtt(nn);
          alpha_r(nn,3, 1) = alpha_r(nn,3, 1) + ndtt(nn) * (toik(i, k) - 1);
          alpha_r(nn,2, 2) = alpha_r(nn,2, 2) + ndtt(nn) * ex2(nn);
          alpha_r(nn,2, 3) = alpha_r(nn,2, 3) + ndtt(nn) * (ex3(nn) - coik(i, k) * ex(nn));
          alpha_r(nn,1, 4) = alpha_r(nn,1, 4) + ndt(nn) * (ex3(nn) * (ex2(nn) - 2) - ex(nn) * (3 * ex2(nn) - 3 + coik(i, k)) * coik(i, k));
        end
      end
    end
  end
end
  
 %Calculate mixture contributions
 for nn = 1 : size(x,1)
  for i = 1 :(NcGERG - 1)
    if x(nn,i) > 0
      for j = (i + 1) : NcGERG
        if x(nn,j) > 0
          mn = mNumb(i, j);
          if mn >= 0
            xijf = x(nn,i) * x(nn,j) * fij(i, j);
            for k = 1 : kpolij(mn)
              ndt(nn) = xijf * delp(nn,dijk(mn, k)) * taupijk(nn,mn, k);
              ndtd(nn) = ndt(nn) * dijk(mn, k);
              alpha_r(nn,1, 2) = alpha_r(nn,1, 2) + ndtd(nn);
              alpha_r(nn,1, 3) = alpha_r(nn,1, 3) + ndtd(nn) * (dijk(mn, k) - 1);
              if iprop > 0
                ndtt(nn) = ndt(nn) * tijk(mn, k);
                alpha_r(nn,1, 1) = alpha_r(nn,1, 1) + ndt(nn);
                alpha_r(nn,2, 1) = alpha_r(nn,2, 1) + ndtt(nn);
                alpha_r(nn,3, 1) = alpha_r(nn,3, 1) + ndtt(nn) * (tijk(mn, k) - 1);
                alpha_r(nn,2, 2) = alpha_r(nn,2, 2) + ndtt(nn) * dijk(mn, k);
                alpha_r(nn,2, 3) = alpha_r(nn,2, 3) + ndtt(nn) * dijk(mn, k) * (dijk(mn, k) - 1);
                alpha_r(nn,1, 4) = alpha_r(nn,1, 4) + ndtd(nn)* (dijk(mn, k) - 1) * (dijk(mn, k) - 2);
              end
            end
            for k = (1 + kpolij(mn)):(kpolij(mn) + kexpij(mn))
              cij0 = cijk(mn, k) * delp(nn,2);
              eij0 = eijk(mn, k) * del(nn);
              ndt(nn) = xijf * nijk(mn, k) * delp(nn,dijk(mn, k)) * exp(cij0 + eij0 + gijk(mn, k) + tijk(mn, k) * lntau(nn));
              ex(nn) = dijk(mn, k) + 2 * cij0 + eij0;
              ex2(nn) = (ex(nn) * ex(nn) - dijk(mn, k) + 2 * cij0);
              alpha_r(nn,1, 2) = alpha_r(nn,1, 2) + ndt(nn) * ex(nn);
              alpha_r(nn,1, 3) = alpha_r(nn,1, 3) + ndt(nn) * ex2(nn);
              if iprop > 0
                ndtt(nn) = ndt(nn) * tijk(mn, k);
                alpha_r(nn,1, 1) = alpha_r(nn,1, 1) + ndt(nn);
                alpha_r(nn,2, 1) = alpha_r(nn,2, 1) + ndtt(nn);
                alpha_r(nn,3, 1) = alpha_r(nn,3, 1) + ndtt(nn) * (tijk(mn, k) - 1);
                alpha_r(nn,2, 2) = alpha_r(nn,2, 2) + ndtt(nn) * ex(nn);
                alpha_r(nn,2, 3) = alpha_r(nn,2, 3) + ndtt(nn) * ex2(nn);
                alpha_r(nn,1, 4) = alpha_r(nn,1, 4) + ndt(nn) * (ex(nn) * (ex2(nn) - 2 * (dijk(mn, k) - 2 * cij0)) + 2 * dijk(mn, k));
              end
            end
          end 
        end 
      end
    end
  end
 end
  
end 