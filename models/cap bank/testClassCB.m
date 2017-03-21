% Script for testing the capacitor bank asset
clc; clear; close all;

% General setup
% Base parameters for converting to p.u. 
Sbase = 400e3;           % Complex power base [VA].
Vbase = 0.4e3;           % Voltage base [V].
Ts = 60;                 % Sampling time [sec].
N = 24*2*60*60/Ts;       % Number of samples [-].
day = 1;                 % Julian day (1-365).

% Setup capacitor bank models
param.sBase = Sbase;
param.vBase = Vbase;
param.qMax = 50e3;
param.onPU = false;
param.Ts = Ts;
param.nSteps = 10;

% Create object
CB = cbAsset(param);
CB.setQmode(1);

% Allocate memory.
v = 400+400*0.1*sin((1:N)/300);
Q = zeros(1,N);
uRef = [0*ones(1,500) 10*ones(1,500) 0*ones(1,500) 1*ones(1,500) 3*ones(1,1000)];
vRef =400*ones(1,N);
for i=1:N
    if ~mod(i,24*60*60/Ts)
       if day == 365
           day = 1;
       else
           day = day +1;
       end
    end
    
    Q(i) = CB.sample(v(i),vRef(i),uRef(i),i,day);
end

%% Plotting
t = (0:(N-1))/(60*60/Ts);

figure
subplot(2,1,1)
plot(t,Q/1e3)
grid
ylabel('Reactive Power [VAR]')
subplot(2,1,2)
plot(t,v/400)
ylabel('Voltage [PU]')
xlabel('Time [hrs]')
grid
