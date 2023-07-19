function [taup0, taup, taupijk] = tTermsGERG(tau, lntau, x)
%Private Sub tTermsGERG(ByVal tau As Double, ByVal lntau As Double, x() As Double)
%Private Sub tTermsGERG(tau, lntau, x)

%Calculate temperature dependent parts of the GERG-2008 equation of state

  %Dim i As Integer, j As Integer, k As Integer, mn As Integer, taup0(12) As Double
NcGERG = size(x,2);

load('toik.mat','toik');
load('kpol.mat','kpol');
load('kexp.mat','kexp');
load('noik.mat','noik');

  i = 5;  %Use propane to get exponents for short form of EOS
  for k = 1 : (kpol(i) + kexp(i))
    taup0(:,k) = exp(toik(i, k) .* lntau);
  end
% taup0=ones(1,12);
for nn = 1 : size(x,1)
  for i = 1 : NcGERG
    if x(nn,i) > 0 
      if i > 4 && i ~= 15 && i ~= 18 && i ~= 20
        for k = 1 : (kpol(i) + kexp(i))
          taup(nn,i, k) = noik(i, k) * taup0(nn,k);
        end
      else
        for k = 1 : (kpol(i) + kexp(i))
          taup(nn,i, k) = noik(i, k) * exp(toik(i, k) .* lntau(nn));
        end
      end
    end
  end
end
  
  load('mNumb.mat','mNumb');
  load('kpolij.mat','kpolij');
  load('nijk.mat','nijk');
  load('tijk.mat','tijk');
  taupijk=[];
  for nn = 1 : size(x,1)
  for i = 1 : (NcGERG - 1)
    if x(nn,i) > 0
      for j = (i + 1) : NcGERG
        if x(nn,j) > 0 
          mn = mNumb(i, j);
          if mn >= 0 
            for k = 1 : kpolij(mn)
              taupijk(nn,mn, k) = nijk(mn, k) * exp(tijk(mn, k) * lntau(nn));
            end
          end
        end
      end
    end
  end
  end
if max(size(taupijk))==0
    taupijk=0;
end
end 