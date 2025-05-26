function [MM]=MolarMassGERG1(y)
%x = MOLEFRAC as a row vector %% concentrazioni MOLARI
%MOLEFRAC normalized to 1
y=y/100;
% [a b c]=size(x);
load('MMiGERG.mat','MMiGERG');
MMiGERG=MMiGERG';
% MMiGERG_vect=padarray(MMiGERG,c-1,'symmetric','pre')


%MM=sum(x'.*MMiGERG)';
y=permute(y,[2,1,3]);
for ii=1:size(y,3)
    
MM(ii,:)=MMiGERG*y;

% %Molar masses (g/mol)
%   MMiGERG(1) = 16.04246;    %Methane
%   MMiGERG(2) = 28.0134;     %Nitrogen
%   MMiGERG(3) = 44.0095;     %Carbon dioxide
%   MMiGERG(4) = 30.06904;    %Ethane
%   MMiGERG(5) = 44.09562;    %Propane
%   MMiGERG(6) = 58.1222;     %Isobutane
%   MMiGERG(7) = 58.1222;     %n-Butane
%   MMiGERG(8) = 72.14878;    %Isopentane
%   MMiGERG(9) = 72.14878;    %n-Pentane
%   MMiGERG(10) = 86.17536;   %Hexane
%   MMiGERG(11) = 100.20194;  %Heptane
%   MMiGERG(12) = 114.22852;  %Octane
%   MMiGERG(13) = 128.2551;   %Nonane
%   MMiGERG(14) = 142.28168;  %Decane
%   MMiGERG(15) = 2.01588;    %Hydrogen
%   MMiGERG(16) = 31.9988;    %Oxygen
%   MMiGERG(17) = 28.0101;    %Carbon monoxide
%   MMiGERG(18) = 18.01528;   %Water
%   MMiGERG(19) = 34.08088;   %Hydrogen sulfide
%   MMiGERG(20) = 4.002602;   %Helium
%   MMiGERG(21) = 39.948;     %Argon
    
end