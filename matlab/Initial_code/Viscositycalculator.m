function out=Viscositycalculator(MOLEFRAC,T)
%% Inputs 
% X - vector of volume fractions for next componemt of gas mixture:

X(:,1)=MOLEFRAC(:,2);       %   X(1)  - volume fraction of N2   [m^3.m^-3]
X(:,2)=MOLEFRAC(:,16);      %   X(2)  - volume fraction of O2   [m^3.m^-3]
X(:,3)=MOLEFRAC(:,3);       %   X(3)  - volume fraction of CO2  [m^3.m^-3]
X(:,4)=MOLEFRAC(:,18);      %   X(4)  - volume fraction of H2O  [m^3.m^-3]
X(:,5)=zeros(size(T));      %   X(5)  - volume fraction of SO2  [m^3.m^-3]
X(:,6)=MOLEFRAC(:,17);       %   X(6)  - volume fraction of CO   [m^3.m^-3]
X(:,7)=MOLEFRAC(:,15);       %   X(7)  - volume fraction of H2   [m^3.m^-3]
X(:,8)=MOLEFRAC(:,1);       %   X(8)  - volume fraction of CH4  [m^3.m^-3]
X(:,9)=zeros(size(T));      %   X(9)  - volume fraction of C2H4 [m^3.m^-3]
X(:,10)=MOLEFRAC(:,4);      %   X(10) - volume fraction of C2H6 [m^3.m^-3]
X(:,11)=MOLEFRAC(:,5);      %   X(11) - volume fraction of C3H8 [m^3.m^-3]
X(:,12)=MOLEFRAC(:,7);      %   X(12) - volume fraction of n-C4H10[m^3.m^-3]
X(:,13)=MOLEFRAC(:,19);     %   X(13) - volume fraction of H2S  [m^3.m^-3]
X(:,14)=MOLEFRAC(:,6);     %   X(14) - volume fraction of i-C4H10  [m^3.m^-3]
X(:,15)=MOLEFRAC(:,9);     %   X(15) - volume fraction of n-C5H12  [m^3.m^-3]
X(:,16)=MOLEFRAC(:,8);     %   X(16) - volume fraction of i-C5H12  [m^3.m^-3]
X(:,17)=sum(MOLEFRAC(:,[10:14]),2);     %   X(17) - volume fraction of C6  [m^3.m^-3]

%     T - Temperature of mixture                                   [K]


%% Description
%
% Function for dynamic viscozity calculation 
% of gas with composition X and temperature T.
%
%       mu1 . (T1 + 1,47.Tc)       T^1,5
% mu = ---------------------- . ---------- . 10^(-6) . X
%           T + 1,47.Tc            T1^1,5
%

%% Used coeficients and functions
%
% TK - constant for conversion Celsius degrees to Kelvins      [K]
% mu1(1).. mu1(13) - dynamic viscosity of gas component 
%                    for temperature T1(1)...T(13)             [Pa.s]
% Tc(1) .. Tc(13)  - boiling temperature of gas component      [K]
%
%% Output
%
% mu - Dynamic viscosity of gas mixture                        [Pa.s]
%      defined for temperature T =<273,15;2500>.
%
%% Copyright (C) 2011
% Authors     : Kukurugya Jan, Terpak Jan
% Organization: Technical University of Kosice
% e-mail      : jan.kukurugya@tuke.sk, jan.terpak@tuke.sk
% Revision    : 28.12.2011
%
%     1      2      3      4     5      6     7      8      9      10    11    12    13
%     N2     O2     CO2    H2O   SO2    CO    H2     CH4    C2H4   C2H6  C3H8  C4H10 H2S   
%mu1=[17.499 20.194 14.568 14.2  12.626 16.6  8.78   10.745 10.152 9.181 8.181 7.359 12.674];
%T1 =[20.0   20.0   20.0   20.0  20.0   0.0   20.0   20.0   20.0   20.0  28.0  20.0  33.0];
%Tc =[-196.0 -183.0 -78.5  99.63 -10.0 -191.5 -253.0 -258.7 -103.7 -88.0 -42.1 -0.5  -60.0];
%     1          2         3        4         5        6        7         8          9        10        11         12         13        14       15         16        17
%     N2         O2       CO2      H2O       SO2      CO        H2        CH4       C2H4      C2H6      C3H8     n-C4H10      H2S     i-C4H10  n-C5H12    i-C5H12    C6
mu1=[17.580    20.194    14.689   13.022    12.626   17.426    8.8129   11.024     10.152    9.2071    8.0115     7.359      12.674    7.3744    7.2967    7.3340    7.6545 ];
T1 =[20.0       20.0      20.0    120.0      20.0     20.0     20.0      20.0      20.0      20.0      20.0       20.0        33.0      20.0      40.0      28.0      70.0  ];
Tc =[-195.795  -182.962   -78.4    99.9743   -10.0   -191.51  -252.754  -161.483   -103.7    -88.598   -42.09     -0.49       -60.0    -11.749    36.06     27.82      68.71 ];
TK=273.15;

  for i=1:1:17
    T11(1:size(X(:,1),1),i)=T1(i)+TK;    
    Tc1(1:size(X(:,1),1),i)=Tc(i)+TK;    
  end

  T(find(T<TK))=TK;
  T(find(T>2500))=2500;

  for i=1:1:17
    muT(1:size(X(:,1),1),i)=((mu1(i)*10^(-6)).*X(:,i).*(T11(:,i) + 1.47*Tc1(:,i)).*T.^1.5)./((T + 1.47*Tc1(:,i)).*T11(:,i).^1.5);
  end
  out=sum(muT')';
end