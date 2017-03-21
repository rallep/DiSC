% Test class
clc; clear; close all;

% Base parameters for converting to p.u. 
Sbase = 400e3;          % Complex power base [VA].
Vbase = 0.4e3;          % Voltage base [V].
Ts = 1;                 % Sampling time [sec].
N = 24*60*60/Ts;        % Number of samples [-].

% Setup static compensator model
param.sBase = Sbase;
param.vBase = Vbase;
param.qRatedMax = 50e6;
param.qRatedMin = -20e6;
param.onPU = false;

% Create object
SC = scAsset(param);
SC.setQmode(1);         % Set reactive power control mode

% Voltage
v = 400+400*0.1*sin((1:N)/1000);

% Allocate memory
Q = zeros(1,N);
qRef = [50e6*ones(1,N/2) ones(1,N/2)*(-50e6)];
vRef =400*ones(1,N);
for i=1:N
    Q(i) = SC.sample(v(i),qRef(i),vRef(i));
end

%% Plotting
t = (0:(N-1))/(60*60/Ts);

figure
subplot(2,1,1)
plot(t,v)
ylabel('Voltage [V]')
subplot(2,1,2)
plot(t,Q/1e6)
grid
xlabel('Time [hrs]')
ylabel('Reactive Power [MVAR]')