% Script for testing the PV system model
clc; close all; clear;

% General input
rng(1)      % Seed random generator
Ts = 10;    % Sampling time

% Base parameters for converting to p.u. 
Sbase = 50e6;          % Complex power base [VA].
Vbase = 20e3;          % Voltage base [V].

% Setup PV object
param.sBase = Sbase;
param.vBase = Vbase;
param.pRated = 6e3;
param.sMax = 6e3;
param.eta = 0.25;
param.A = 27;
param.onPU = false;
param.Ts = Ts;

% Create PV object
PV = pvAsset(param);

% Create solar irradiance object
% Parameters
param.lat = 56.889;             % Latitude for Sørup (degrees)
day = 181;                      % Julian date, 1-365 (181=June 30)
param.t = 0.75;                 % Transmittance (unitless)
param.p = 100;                  % Air pressure (Kpa)
param.Ts = Ts;                  % Sampling time
SI = solarIrradiance(param);

% Setup simulation
N = 1*24*60*60/Ts;
% Cloud cover
cc = [0.5*ones(1,500*60/Ts) 0.5+0.5*sin(0.007/(60/Ts)*(1:300*60/Ts)-pi) 0.1*ones(1,800*60/Ts)...
    1*ones(1,700*60/Ts) 0.5+0.5*sin(0.007/(60/Ts)*(1:650*60/Ts)+pi/2) 0*ones(1,4000*60/Ts)...
    0.5+0.5*sin(0.005*60/Ts*(1:1000*60/Ts)-pi)];
cc = [1*ones(1,1000) 0.5+0.5*sin(0.004/(60/Ts)*(1:5740-1000)+pi/2) zeros(1,N)];


% Input to PV
dP = 0;                % Change in power reference
dPlim = 0;             % Derate power
qRef = 0;              % Reactive power reference
vRef = 400;            % Voltage reference
PV.setQmode(0);        % Qmode: 0 = constant power factor, 1 = follow Q reference, 2 = voltage droop control
PV.setPF(1);         % Set PV power factor

% Allocate memory
p = zeros(1,N);
q = zeros(1,N);
pAva = zeros (1,N);
Go = zeros(1,N);
for i=1:N
    Go(i) = SI.sample(i,day,cc(i));
    [p(i),q(i),pAva(i)] = PV.sample(Go(i),vRef,dP,dPlim,qRef,vRef);

end

%% Plotting
t = (0:N-1)/(60*60/Ts);
figure
subplot(2,1,1)
[h,hl1,hl2] = plotyy(t,Go,t,cc(1:length(t)));
set(hl2,'LineStyle','--','LineWidth',1.2)
ylabel(h(1),'Solar Irradiance [W/m^2]') % left y-axis
ylabel(h(2),'Cloud Cover [-]') % right y-axis
ylim([0 1000])
set(gca,'XTick',[0,6,12,18,24])
set(gca,'XTickLabel',{})

subplot(2,1,2)
plot(t,p/1e3,t,pAva/1e3,'--',t,q/1e3,'--')
ylabel('Power [kW]')
set(gca,'XTick',[0,6,12,18,24])
set(gca,'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
xlabel('Time of Day')
legend('Power Output','Available Power','Reactive Power','Location','SouthEast')

