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
param.pRate = 50;
param.onPU = false;

% Create object
flexA = flexAsset(param);

% Allocate memory
P = zeros(1,N);
Q = zeros(1,N);
pFlexUp = zeros(1,N);
pFlexDown = zeros(1,N);
qFlexUp = zeros(1,N);
qFlexDown = zeros(1,N);
satUp = zeros(1,N);
satDown = zeros(1,N);
% Set references
Pref = [0*ones(1,100) 3000*ones(1,500) -3000*ones(1,200) 0*ones(1,N)];
qRef = 200;

% Voltage

% Simulation
for i=1:N
    [P(i),Q(i),satUp(i),satDown(i),pFlexUp(i),pFlexDown(i),qFlexUp(i),qFlexDown(i)] = flexA.sample(Pref(i),qRef);
end

%% Plotting
t = (0:(N-1));

figure
subplot(3,1,1)
plot(t,P,t,Q)
grid
ylabel('P/Q [W/VAR]')
legend('P','Q')
subplot(3,1,2)
plot(t,pFlexUp,t,pFlexDown,t,qFlexUp,t,qFlexDown)
grid
ylabel('Power [W]')
legend('pFlexUp','pFlexDown','qFlexUp','qFlexDown')
subplot(3,1,3)
plot(t,satUp,t,satDown)
ylabel('Saturation [.]')
legend('Sat Up','Sat Down')
grid
xlabel('Time [hrs]')