function X=MOL2MASS_CONV(Y)

MM=MolarMassGERG1(Y);
Y=Y/100;
% [a b c]=size(x);
load('MMiGERG.mat','MMiGERG');
% global MMiGERG
MMiGERG=MMiGERG';

%Y=permute(Y,[2,1,3]);

X=Y.*MMiGERG./MM'*100;