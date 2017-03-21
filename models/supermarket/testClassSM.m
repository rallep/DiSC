% Script for testing the supermarket with flexible defrost cycles class
clc; close all; clear;

% Simulation parameters
Ts = 60;
N = 2*24*60*60/Ts; % Number of iterations
sBase = 50e6;

% Set parameters
param.sBase = sBase;
param.Ts = Ts;
param.pDefrost = [3e3 10e3];
param.onPU = false;

numDF = length(param.pDefrost);

sm = smAsset(param);

% Input
u = zeros(numDF,1);
u(1) = 1;

% Load data
pRS =10e3;

% Simulate
p = zeros(1,N);
q = zeros(1,N);
timeReady = zeros(numDF,N);
timeDefrost = zeros(numDF,N);

for k=1:N
    [p(k),q(k),timeReady(:,k),timeDefrost(:,k)] = sm.sample(u,pRS);

end

%% Plotting
tvec = 0:N-1;
figure
plot(tvec,p/1000,tvec,q/1000)
ylabel('Power [kW]')
xlabel('Time [min]')

figure
plot(tvec,timeReady,tvec,timeDefrost)
