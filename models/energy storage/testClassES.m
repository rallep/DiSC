% Test of class
clc; clear; close all;

% Base parameters for converting to p.u. 
Sbase = 400e3;           % Complex power base [VA].
Vbase = 0.4e3;           % Voltage base [V].
Ts = 60;                 % Sampling time [sec].
N = 24*60*60/Ts;         % Number of samples [-].

% Set parameters
param.sBase = Sbase;
param.vBase = Vbase;
param.sMax = 3e3;
param.pRatedMax = 3e3;
param.pRatedMin = -3e3;
param.pRate = 50;
param.eRated = 60e3*60*60; % 60 kWh
param.Ts = Ts;
param.onPU = false;

% Create object
ES = esAsset(param);
ES.setPF(1);
ES.setDrain(1);
ES.setPmode(0);
ES.setQmode(1);

% Allocate memory
P = zeros(1,N);
Q = zeros(1,N);
E = zeros(1,N);
% Set references
Pref = [0*ones(1,100) -3000*ones(1,500) 3000*ones(1,200) 0*ones(1,N)];
qRef = 0;
vRef = 400*ones(1,N);

% Voltage
v = 400+0.1*400*sin((1:N)*0.008);

% Simulation
for i=1:N
    [P(i),Q(i),E(i)] = ES.sample(v(i),Pref(i),qRef,vRef(i));
end

%% Plotting
t = (0:(N-1));

figure
subplot(3,1,1)
plot(t,Pref(1:N)/1e3,t,P/1e3,t,Q)
grid
ylabel('P/Q [kW/kVAR]')
legend('Pref','P')
subplot(3,1,2)
plot(t,E/(1e3*60*60))
grid
ylabel('Energy [kWh]')
subplot(3,1,3)
plot(t,v/400)
ylabel('Voltage [PU]')
grid
xlabel('Time [hrs]')