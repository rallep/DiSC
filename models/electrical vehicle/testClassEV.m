% Script for testing the EV model
clc; close all; clear;

% Setup system
rng(1)                  % Seed random generator
% Base parameters for converting to p.u. 
Sbase = 400e3;          % Complex power base [VA].
Vbase = 0.4e3;          % Voltage base [V].
Ts = 60;                % Sample time [sec].

% Set parameters
param.sBase = Sbase;
param.vBase = Vbase;
param.sMax = 6e3;
param.pRated = 6e3;
param.eRated = 65e3*60*60;  % 65 kWh.
param.onPU = false;
param.Ts = Ts;
param.pRate = 20;

% Object
ev = evAsset(param);

% Input
day = 1;
pRef = 0.9*param.pRated;
qRef = 1000;
ev.setPmode(0);
ev.setQmode(2);
vRef = 400;


% Simulation time
N = 5*24*60*60/Ts;

% Simulate
p = zeros(1,N);
q = zeros(1,N);
e = zeros(1,N);
away = zeros(1,N);

v = 400+0.1*400*sin((1:N)*0.008);
for i=1:N
    % Itterate day
    if ~mod(i,24*60*60/Ts)
       if day == 365
           day = 1;
       else
           day = day +1;
       end
    end
    [p(i),q(i),e(i),away(i)] = ev.sample(i,day,v(i),pRef,qRef,vRef);
end

%% Plotting
t = (0:N-1)/(60*60/Ts);

figure
subplot(3,1,1)
plot(t,p/1e3,t,q/1e3)
ylabel('P/Q [kW/kVAR]')
subplot(3,1,2)
plot(t,e/(1e3*60*60))
ylabel('Energy [kWh]')
subplot(3,1,3)
plot(t,v/400)
ylabel('Voltage [PU]')
xlabel('Time [hrs]')