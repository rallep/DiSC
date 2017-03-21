% Test solarIrradiance class
clc; close all; clear;
rng(1);

% Input
cc = 0;

% Parameters
lat = 56.889;       % Latitude for Sørup (degrees)
day = 50;           % Julian date, 1-365 (181=June 30)
t = 0.75;           % Transmittance (unitless)
S = 1367;           % Solar constant (w/m^2)
p = 100;            % Air pressure (Kpa)
Ts = 60;            % Sampling time
N = 24*1;           % Number of samples

% Setup solar irradiance model
param.lat = lat;
param.t = t;
param.p = p;
param.Ts = Ts;
SolarIrr = solarIrradiance(param);

% Time vector
t = 0:N*60*60/Ts-1;
% Allocate memory
Go = zeros(N*60*60/Ts-1,1);
figure
for d=1:365
    d
    for i=1:N*60*60/Ts
        Go(i) = SolarIrr.sample(i,d,cc); 
    end
    % Animate one year of solar irradiance
    plot(t/(60*60/Ts),Go)
    ylim([0 1000])
    xlim([0 24])
    drawnow
end

%% Plotting
figure
plot(t/(60*60/Ts),Go)
xlabel('Time [hrs]')
ylabel('Irradiance [W,m^2]')
