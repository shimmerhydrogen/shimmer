function [D, ierr, herr, alpha_r] = DensityGERG(iFlag, T, P, x, dimn,Tr,Dr,Tcx,Dcx,Vcx)
%Sub DensityGERG(iFlag As Integer, ByVal T As Double, ByVal P As Double, x() As Double, ByRef D As Double, ByRef ierr As Integer, ByRef herr As String)

% Sub DensityGERG(iFlag, T, P, x, D, ierr, herr)
% 
% Calculate density as a function of temperature and pressure.  This is an iterative routine that calls PressureGERG
% to find the correct state point.  Generally only 6 iterations at most are required.
% If the iteration fails to converge, the ideal gas density and an error message are returned.
% No checks are made to determine the phase boundary, which would have guaranteed that the output is in the gas phase (or liquid phase when iFlag=2).
% It is up to the user to locate the phase boundary, and thus identify the phase of the T and P inputs.
% If the state point is 2-phase, the output density will represent a metastable state.
% 
% Inputs:
%  iFlag - set to 0 for strict pressure solver in gas phase without checks (fastest mode, but output state may not be stable single phase)
%          set to 1 to make checks for possible 2-phase state (result may still not be stable single phase, but many unstable states will be identified)
%          set to 2 to search for liquid phase (and make the same checks when iFlag=1)
%      T - temperature (K)
%      P - pressure (kPa)
%    x() - composition (mole fraction)
% (An initial guess for the density can be sent in D as the negative of the guess for roots that are in the liquid phase instead of using iFlag=2)
% 
% Outputs:
%      D - density (mol/l)
%   ierr - Error number (0 indicates no error)
%   herr - Error message if ierr is not equal to zero

% %   Dim it As Integer, nFail As Integer, iFail As Integer
% %   Dim plog As Double, vlog As Double, P2 As Double, Z As Double, dpdlv As Double, vdiff As Double, tolr As Double, vinc As Double
% %   Dim Tcx As Double, Dcx As Double
% %   Dim dPdD As Double, d2PdD2 As Double, d2PdTD As Double, dPdT As Double, U As Double, H As Double, S As Double, A As Double
% %   Dim Cv As Double, Cp As Double, W As Double, G As Double, JT As Double, Kappa As Double, PP As Double

% global MaxFlds;
% MaxFlds=size(x,2);
% global xold;
% %   for i = 1 : MaxFlds
% %     xold(i) = 0;
% %   end
% xold=zeros(size(x));


  ierr = 0;
  herr = '';
  nFail = 0;
  iFail = 0;
  if any(P == 0)
  D(find(P==0)) = 0 ;
  return
  end
  
  D=zeros(size(P));
  RGERG = 8.314472;
  index=find(sum(x)>0);
  
  tolr = 1e-9;
%   [Tcx,Dcx]=PseudoCriticalPointGERG(x,dimn,index); %Call PseudoCriticalPointGERG(x, Tcx, Dcx)
  
  if any(D >= 0)
    D(find(D>=0)) = P(find(D>=0))./RGERG./T(find(D>=0));               %Ideal gas estimate for vapor phase
    if iFlag == 2
        D(find(D>=0)) = Dcx(find(D>=0)) * 3;   %Initial estimate for liquid phase
    end
  else
        D(find(D<0)) = abs(D(find(D<0))) ;                     %If D<0, then use as initial estimate
  end
  plog = log(P);
  vlog = -log(D);
  
    
  for it = 1 : 50
    if any(vlog < -7) ||any(vlog > 100) || it == 20 || it == 30 || it == 40 || iFail == 1
      %Current state is bad or iteration is taking too long.  Restart with completely different initial state
      iFail = 0;
      if nFail > 2
          DError(P,RGERG,T);
          return
      end
      nFail = nFail + 1;
      if nFail == 1 
         D = Dcx * 3 ;   %If vapor phase search fails, look for root in liquid region
      elseif nFail == 2
        D = Dcx * 2.5;  %If liquid phase search fails, look for root between liquid and critical regions
     elseif nFail == 3 
        D = Dcx * 2;    %If search fails, look for root in critical region
      end
      vlog = -log(D);
    end
    D = exp(-vlog);
    [P2,Z,dPdDsave,alpha_r] = PressureGERG(T,D,x,dimn,Tr,Dr); %Call PressureGERG(T, D, x, P2, Z)
    if any(dPdDsave) < 0 || any(P2) <= 0
      %Current state is 2-phase, try locating a different state that is single phase
      vinc = 0.1;
      if D > Dcx
         vinc = -0.1;
      end
      if it > 5 
         vinc = vinc/2;
      end
      if it > 10 && it < 20
         vinc = vinc/5;
      end
      vlog = vlog + vinc;
    else
      %Find the next density with a first order Newton's type iterative scheme, with
      %log(P) as the known variable and log(v) as the unknown property.
      %See AGA 8 publication for further information.
      dpdlv = -D .* dPdDsave;     %d(p)/d[log(v)]
      vdiff = (log(P2) - plog) .* P2 ./ dpdlv;
      vlog = vlog - vdiff;
      if max(abs(vdiff)) < tolr
        %Check to see if state is possibly 2-phase, and if so restart
        if any(dPdDsave) < 0
          iFail = 1;
        else
          D = exp(-vlog);
          Converged(iFlag,P,RGERG,T);
%           [ierr, herr, D]=Converged(iFlag,P,RGERG,T);
          return
        end
      end
    end
  end
  herr = '';
  ierr = 0;

%Iteration failed (above loop did find a solution or checks made below indicate possible 2-phase state)
function [ierr, herr, D] = DError(P,RGERG,T)
  ierr = 1;
  herr = "Calculation failed to converge in GERG method, ideal gas density returned.";
  D = P./ RGERG./ T;
  return
end

%Iteration converged
function []=Converged(iFlag,P,RGERG,T)
%     iFlag=0;
%     herr = '';
%     ierr = 0;
%  %NumbItGERG = it    'Used to pass back the number of iterations for information only

  %If requested, check to see if point is possibly 2-phase
  if iFlag > 0
%     [PP, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A] = PropertiesGERG(T, D, x);
    if any(PP <= 0) || any(dPdD <= 0) || any(d2PdTD <= 0)
        [ierr, herr, D]=DError(P,RGERG,T);
        ierr=0;
        return
    end
    if any(Cv <= 0) || any(Cp <= 0) || any(W <= 0)
        [ierr, herr, D]=DError(P,RGERG,T);
        ierr=0;
        return
    end
  end 
end
  end
   