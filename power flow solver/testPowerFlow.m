% Script for testing the power flow module
clc; close all; clear;

% Construct grid
Z = [1 2 0.1 0.1 1;
     2 3 0.1 0.1 1;
     2 4 0.1 0.1 1;
     3 4 0.1 0.1 1;
     3 5 0.1 0.1 1;
     5 6 0.1 0.1 1;
     6 7 0.1 0.1 1;
     6 8 0.1 0.1 1;
     7 9 0.1 0.1 1;];
nBus = 9;
% Setup load flow object
param.type = [0 ones(1,nBus-1)];
param.vBase = 400*ones(1,nBus);

% Construct power flow object
pf = powerFlow(param);
 
% Get admittance matrix
Y = pf.lDataToY(Z);
Z = pf.zbuild(Z);

% Setup busses
PF = 0.95;          % Power factor
Pbus1 = 0;
Pbus2 = -10e3;
Pbus3 = -10e3;
Pbus4 = -10e3;
Pbus5 = -10e3;
Pbus6 = -10e3;
Pbus7 = -10e3;
Pbus8 = -10e3;
Pbus9 = -10e3;

Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7 Pbus8 Pbus9];

Qin=Pin.*tan(acos(PF));

Vin = [400 zeros(1,nBus-1)]; % Bus1 = slack
% Solve load flow
[vOut,pSlack,qSlack,~]=pf.nrLoadFlow(Y,Pin',Qin')
