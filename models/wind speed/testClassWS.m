% Script for testing wind speed class
clc; clear; close all;
rng(1);

% Setup
param.Ts = 60*10;
param.z = 35;
vMean = 10;
ws = windSpeed(param);

%% Simulate
N = 3000;
w = zeros(1,N);
ws.setNumVdhSpec(20);

for i=1:N
    w(i) = ws.sample(i,vMean);
end

%% Plotting
figure
plot(w)
ylabel('Wind Speed [m/s]')
xlabel('Time [samples]')