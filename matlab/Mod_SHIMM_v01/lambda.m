%%% HOFER 1973 explicit approximations for Coolebrook-White Equation %%%
function [lambda]=lambda(Reyn,D,epsi)

lambda=(2.*log10(4.518./Reyn.*log10(Reyn./7)+epsi./(3.71.*D))).^(-2);

end


