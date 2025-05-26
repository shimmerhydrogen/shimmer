function [a0] = Alpha0GERG(T, D, x)
%Private Sub Alpha0GERG(ByVal T As Double, ByVal D As Double, x() As Double, a0() As Double)
% Private Sub Alpha0GERG(T, D, x, a0)
% 
% Calculate the ideal gas Helmholtz energy and its derivatives with respect to tau and delta.
% This routine is not needed when only P (or Z) is calculated.
% 
% Inputs:
%      T - temperature (K)
%      D - density (mol/l)
%    x() - composition (mole fraction)
% 
% Outputs:
%  a0(0) - ideal gas Helmholtz energy (dimensionless [i.e., divided by RT])
%  a0(1) - tau*partial(a0)/partial(tau)
%  a0(2) - tau^2*partial^2(a0)/partial(tau)^2
  
%   Dim i As Integer, j As Integer
%   Dim LogT As Double, LogD As Double, LogHyp As Double, th0T As Double, LogxD As Double
%   Dim SumHyp0 As Double, SumHyp1 As Double, SumHyp2 As Double
%   Dim em As Double, ep As Double, hcn As Double, hsn As Double
%   a0(0) = 0: a0(1) = 0: a0(2) = 0

  NcGERG = size(x,2);
  load('th0i.mat','th0i');
  load('n0i.mat','n0i');
  
  nn=size(x,1);
  a0 = zeros(nn,3);
  %if D > 0 
  LogD = zeros(size(D));
     LogD(find(D>0)) = log(D(find(D>0)));
%   else
      LogD(find(D<=0)) = log(1E-20);
  %end
  LogT = log(T);
  
  for nn=1 : size(x,1)
  for i = 1 : NcGERG
    if x(nn,i) > 0
      LogxD = LogD(nn) + log(x(nn,i));
      SumHyp0 = 0;
      SumHyp1 = 0;
      SumHyp2 = 0;
      for j = 4 : 7
        if th0i(i, j) > 0 
          th0T = th0i(i, j)./ T(nn);
          ep = exp(th0T);
          em = 1 ./ ep;
          hsn = (ep - em) / 2;
          hcn = (ep + em) / 2;
          if j == 4 || j == 6
            LogHyp = log(abs(hsn));
            SumHyp0 = SumHyp0 + n0i(i, j) * LogHyp;
            SumHyp1 = SumHyp1 + n0i(i, j) * th0T * hcn / hsn;
            SumHyp2 = SumHyp2 + n0i(i, j) * (th0T / hsn) ^ 2;
          else
            LogHyp = log(abs(hcn));
            SumHyp0 = SumHyp0 - n0i(i, j) * LogHyp;
            SumHyp1 = SumHyp1 - n0i(i, j) * th0T * hsn / hcn;
            SumHyp2 = SumHyp2 + n0i(i, j) * (th0T / hcn) ^ 2;
          end
        end
      end
      a0(nn,1) = a0(nn,1) + x(nn,i) * (LogxD + n0i(i, 1) + n0i(i, 2) / T(nn) - n0i(i, 3) * LogT(nn) + SumHyp0);
      a0(nn,2) = a0(nn,2) + x(nn,i) * (n0i(i, 3) + n0i(i, 2) / T(nn) + SumHyp1);
      a0(nn,3) = a0(nn,3) - x(nn,i) * (n0i(i, 3) + SumHyp2);
    end
  end
  end
end