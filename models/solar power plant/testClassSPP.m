% Script for testing the solar power plant class
clc; close all; clear all;

% General setup
% Base parameters for converting to p.u. 
Sbase = 50e6;           % Complex power base [VA].
Vbase = 20e3;           % Voltage base [V].
Ts = 60;                % Sampling time [s].
onPU = false;           % Indicate if simulation is on per unit
rng(1);                 % Seed random generator
day = 160;              % Julian day of the year

% Set parameters
param.sBase = Sbase;
param.vBase = Vbase;
param.pRated = 5e3;
param.sMax = 5e3;
param.eta = 0.25;
param.A = 20;
param.numPv = 10;
param.Ts = Ts;
param.onPU = onPU;
param.lat = 56.889;         
param.t = 0.75;           
param.p = 100;

% Construct SPP
SPP = sppAsset(param);

% Input WPP
dP = 0;                 % Change in power reference
dPlim = 0;              % Power derate
qRef = 0;               % Reactive power reference
vRef = Vbase;           % Voltage reference 

SPP.setQmode(0);        % Set reactive control mode
SPP.setPF(0.99);         % Set power factor

% Cloud cover, number of samples and voltage
cc = 0;
N = 24*60*60/Ts;
v = 20e3+0.1*20e3*sin((1:N)*0.008);

% Simulate
p = zeros(1,N);
q = zeros(1,N);
pAva = zeros(1,N);
for i=1:N
    [p(:,i), q(:,i), pAva(:,i)] = SPP.sample(cc,v(i),dP,dPlim,qRef,vRef,i,day); 
end


%% Plotting
t = (0:N-1)/60*60/Ts;
figure
plot(t,p/1e3,t,q/1e3,t,pAva/1e3)
xlabel('Time [hrs]')
ylabel('P/Q [MW/MVAR]')
legend('P Output','Q Output','P Available')