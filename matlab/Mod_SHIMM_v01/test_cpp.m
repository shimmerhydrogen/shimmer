clear;
clc;

%% Test SHMGERG


x = [0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
P = 5.101325000000000e+06;
T = 2.931500000000000e+02;
iFlag = 0;
dimn = 1;

[Tr,Dr] = ReducingParametersGERG(x);
[Tcx,Dcx,Vcx] = PseudoCriticalPointGERG(x,dimn);

[P1, Z, D, gamma] = PropertiesGERG(iFlag,P,T,x,dimn,Tr,Dr,Tcx,Dcx,Vcx);


%% Test Speed of sound

x = [1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
P = 5.101325000000000e+06;
T = 2.931500000000000e+02;
iFlag = 0;
dimn = 1;

[Tr,Dr] = ReducingParametersGERG(x);
[Tcx,Dcx,Vcx] = PseudoCriticalPointGERG(x,dimn);

[P1, Z, D, gamma] = PropertiesGERG(iFlag,P,T,x,dimn,Tr,Dr,Tcx,Dcx,Vcx);

