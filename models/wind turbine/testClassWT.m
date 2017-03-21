% Script for testing the wind turbine asset
clc; clear all; close all;
rng(1)
% Base parameters for converting to p.u. 
Sbase = 50e6;          % Complex power base [VA].
Vbase = 20e3;          % Voltage base [V].
Ts = 60;

% Set parameters
param.sBase = Sbase;
param.vBase = Vbase;
param.pRated = 15e6;
param.sMax = 15e6;
param.wMin = 3;
param.wMax = 25;
param.wRated = 12;
param.Ts = Ts;
param.onPU = false;

% Object
wt1 = wtAsset(param);

% Input wt1
dP1 = 0;         % Change in power reference
dPlim1 = 0;      % Power curtailment
qRef1 = 0;       % Reactive power reference
vRef1 = Vbase;       % Voltage reference for coupling point
wt1.setQmode(2); % Qmode: 0 = constant power factor, 1 = follow Q reference, 2 = voltage droop control
wt1.setPF(1);    % Set wind turbine power factor

% Setup wind model
param.Ts = Ts;
param.z = 35;
ws = windSpeed(param);

% Wind and voltage
wMean = 10;

N = 1000;
v = 20e3+0.1*20e3*sin((1:N)*0.008);
% Simulate
p1 = zeros(1,N);
q1 = zeros(1,N);
p2 = zeros(1,N);
q2 = zeros(1,N);
w = zeros(1,N);
we = zeros(1,N);
pAva = zeros(1,N);
for i=1:N
    [w(i), we(i)] = ws.sample(i,wMean);
    [p1(i), q1(i), pAva(i)] = wt1.sample(w(i),v(i),dP1,dPlim1,qRef1,vRef1); 
end


%% Plotting
%t = (0:N-1)/(Ts/60);
t = 0:N-1;
figure
subplot(2,1,1)
plot(t*(Ts/60),w)
ylabel('Wind Speed [m/s]')
subplot(2,1,2)
plot(t*(Ts/60),pAva/1e6,t*(Ts/60),p1/1e6,t*(Ts/60),q1/1e6,'--')
ylabel('Power Output [MW]')
ylim([-16 16])
xlabel('Time [min]')
legend('P Avail.','P','Q','Location','SouthWest')
legend('boxoff')