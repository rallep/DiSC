% Newton-Raphson Power flow example from: "Power System Analysis, second
% edition, Hadi Saadat" example 6.7 in book.
clc; clear; close all;
% Grid impedance data
% [from to R X l], from bus [-], to bus [-], resistance [ohm], reactance
% [ohm], length [km].
Z = [1 2 0.02 0.04 1;
     1 3 0.01 0.03 1;
     2 3 0.0125 0.025 1];
 
%% Setup power flow module
% Number of busses
numBus = max(max(Z(:,1)),max(Z(:,2)));
% Set type of each bus (bus 1 is slack bus)
param.type = [0 1 1];  % 0 = slack bus, 1 = PQ bus and 2 = PV bus
% Set base voltage of each bus
param.vBase = [1.05 1 1]; 
% Create power flow module
pFlow = powerFlow(param);
% Get bus Admittance matrix and take a copy of the original. This is used
% for the tap-changing transformers
Y = pFlow.lDataToY(Z);

P2_sch = -256.6/100;
Q2_sch = -110.2/100;
P3_sch = -138.6/100;
Q3_sch = -45.2/100;

Pin = [0 P2_sch P3_sch];
Qin = [0 Q2_sch Q3_sch];

[vOut,pSlack,qSlack,qOut] = pFlow.nrLoadFlow(Y,Pin',Qin')
%out = nrLoadFlow(Y,param.type,Pin',Qin',param.vBase,param.vBase,1e-7,100);
%out.Vout 