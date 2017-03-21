% Script for testing the wind power plant class
clc; close all; clear all;

% General setup
% Base parameters for converting to p.u. 
Sbase = 50e6;          % Complex power base [VA].
Vbase = 20e3;          % Voltage base [V].
Ts = 60;
onPU = false;
rng(1);

% Set parameters
param.sBase = Sbase;
param.vBase = Vbase;
param.pRated = 5e6;
param.sMax = 5e6;
param.wMin = 3;
param.wMax = 25;
param.wRated = 12;
param.Ts = Ts;
param.onPU = onPU;
param.numWt = 3;
param.z = 35;

% Construct WPP
WPP = wppAsset(param);

% Input wt1
dP = 0;             % Change in power reference
dPlim = -3e6;          % Power curtailment
qRef = 0;           % Reactive power reference
vRef = Vbase;       % Voltage reference for coupling point

WPP.setQmode(2);
WPP.setPF(0.9);

% Wind and voltage
wMean = 10;

N = 1000;
v = 20e3+0.1*20e3*sin((1:N)*0.008);

% Simulate
p = zeros(1,N);
q = zeros(1,N);
pAva = zeros(1,N);
for i=1:N
    [p(:,i), q(:,i), pAva(:,i)] = WPP.sample(wMean,v(i),dP,dPlim,qRef,vRef,i); 
end


%% Plotting
t = 0:N-1;
figure
plot(t,p/1e6,t,q/1e6,t,pAva/1e6)
legend('P','Q','P avail.')
ylabel('P/Q [MW/MVAR]')
xlabel('Time [min]')
