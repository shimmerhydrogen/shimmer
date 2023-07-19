function [MHV]=MASS_HV(MASSFRAC)

%MOLEFRAC_INPUT=[xCH4 , xN2 , xCO2 , xC2H6 , xC3H8 , xiC4H10 , xnC4H10 , xiC5H12 , xnC5H12, xC6H14, xC7H16, xC8H18, xC9H20, xC10H22, xH2, xO2, xCO, xH2O , xH2S , xHe, xAr]
%MASSFRAC_INPUT=[xCH4 , xN2 , xCO2 , xC2H6 , xC3H8 , xiC4H10 , xnC4H10 , xiC5H12 , xnC5H12, xC6H14, xC7H16, xC8H18, xC9H20, xC10H22, xH2, xO2, xCO, xH2O , xH2S , xHe, xAr]
load('MMiGERG.mat','MMiGERG');
MM_COMPONENTS=MMiGERG;
%MM_COMPONENTS_HC=[16.04; 30.07; 44.095; 58.12; 58.12; 72.1488; 72.1488; 86.1754];%kg/kmol % 28.15, 44.01, 2, 28.01];

% CALCULATE THE MOLAR HEATING VALUE OF THE EQUIVALENT HYDROCARBON GAS
%Superior (gross) calorific value (kJ/mol)
%891.56  Methane
%1562.14 Ethane
%2221.10 Propane
%2879.76 n-butane
%2869.0 Iso-butane;   source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C75285&Units=SI&Mask=1#Thermo-Gas
%3535.4 n-Pentane;    source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C109660&Units=SI&Mask=1#Thermo-Gas
%3528.4 Iso-Pentane;  source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C78784&Units=SI&Mask=1#Thermo-Gas
%4163.2 Hexane;       source: http://pubchem.ncbi.nlm.nih.gov/compound/8058#section=Decomposition

            %[ xCH4 , xN2 , xCO2  xC2H6,  xC3H8,   xiC4H10,  xnC4H10, xiC5H12, xnC5H12, xC6,]
HHV_MOL_HC = [891.56; 0   ; 0   ; 1562.14; 2221.10; 2869.0;  2879.76;  3528.4; 3535.4; 4163.2];%MJ/kmol

HHV_MASS_HC = HHV_MOL_HC./MM_COMPONENTS(1:10);

%Superior (gross) calorific value (MJ/kg)
% 141.80  Hydrogen
% 10.096  CO
                           %xC7H16, xC8H18, xC9H20, xC10H22, xH2,   xO2, xCO,  xH2O , xH2S , xHe, xAr]
HHV_MASS_ALL = [HHV_MASS_HC;   0  ;   0   ;   0   ;   0    ; 141.80; 0; 10.096; 0    ; 0   ;  0 ; 0  ] ;

MHV=(MASSFRAC*HHV_MASS_ALL)';

% MOLARHV=891.56*xCH4+1562.14*xC2H6+2221.10*xC3H8+2879.76*xC4H10;
% 
% MASSHV=MOLARHV./MCH;
