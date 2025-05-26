function [Tr,Dr] = ReducingParametersGERG(x)
%The following routines are low-level routines that should not be called outside of this code.

%Private Sub ReducingParametersGERG(x, Tr, Dr)
% 
% Calculate reducing variables.  Only need to call this if the composition has changed.
% 
% Inputs:
%    x() - composition (mole fraction)
% 
% Outputs:
%     Tr - reducing temperature (K)
%     Dr - reducing density (mol/l)
global Drold
global Trold

 NcGERG = size(x,2);
  %Check to see if a component fraction has changed.  If x is the same as the previous call, then exit.
  icheck = 1;
  
  global MaxFlds;
MaxFlds=size(x,2);
global xold;
%   for i = 1 : MaxFlds
%     xold(i) = 0;
%   end
xold=zeros(size(x));

  global xold
  %size(xold)
%   for i = 1 : NcGERG
    if norm(abs(x - xold)) > 0.00000001
        icheck = 1;
    end
    xold = x;
%   end
 
  if icheck == 0 
    Dr = Drold;
    Tr = Trold;
%     Exit Sub
return
  end
  
  global Told Trold2;
      
  Told = zeros(size(x,1),1);
  Trold2 = zeros(size(x,1),1);

%Calculate reducing variables for T and D
  Dr = zeros(size(x,1),1);
  Vr = zeros(size(x,1),1);
  Tr = zeros(size(x,1),1);
  
  
  load('gvij.mat','gvij');
  load('bvij.mat','bvij');
  load('btij.mat','btij');  
  load('gtij.mat','gtij'); 
  
  for nn=1:size(x,1)
  for i = 1:NcGERG
    if x(nn,i) > 0 
      F = 1;
      for j = i : NcGERG
        if x(nn,j) > 0 
          xij = F * (x(nn,i) * x(nn,j)) * (x(nn,i) + x(nn,j));
          Vr(nn) = Vr(nn) + xij * gvij(i, j) / (bvij(i, j) * x(nn,i) + x(nn,j));
          Tr(nn) = Tr(nn) + xij * gtij(i, j) / (btij(i, j) * x(nn,i) + x(nn,j));
          F = 2;
        end 
      end
    end
  end
  end
  
  if any(Vr > 0)
  Dr = 1 ./ Vr;
  end
  Drold = Dr;
  Trold = Tr;
end 