function Y=MASS2MOL_CONV(X,MM)

%MM=MolarMassGERG1(Y);
X=X/100;
% [a b c]=size(x);
load('MMiGERG.mat','MMiGERG');
MMiGERG=MMiGERG';
%Y=permute(Y,[2,1,3]);

%X=Y.*MMiGERG./MM'*100;
Y=X./MMiGERG.*MM'*100;
