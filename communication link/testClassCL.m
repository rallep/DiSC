% Script for testing the communication link class
clc; close all; clear;

% Control variables
% The communication link describes both the communication from e.g.,
% controller to asset and from asset to controller. However, the number of
% "packets" might differ. In other words the number of data entries received
% by the controller might not be the same as the ones send from the
% controlelr to the assets.
numLinksIn = 2;                           % Number communication links (data) from e.g., asset to controller
numLinksOut = 2;                          % Number of communication links (data) from e.g., controller to asset    
N=10;                                     % Number of samples                      

% Setup communication link
param.Ts = 60;                          % Sampling time
param.maxDelay = 5*60;                  % Maximum delay
param.numLinksIn = numLinksIn;          % Number of communication links/data packets one way
param.numLinksOut = numLinksOut;        % Number of communication links/data packets the other way
param.mu = 2*60;
param.sigma = 1*60;
% Create communication link object
cl = comLink(param);
cl.setInverseCDFdelay('uniform',param);
cl.setPrLoss(0.01);

Vout = rand(N,numLinksIn);
Vmes = zeros(N,numLinksIn);
Qref = zeros(N,numLinksOut);
%%
for i=1:N
    % Measurements to control communication
    Vmes(i,:) = cl.sampleIn(i,Vout(i,:));
    
    % Controller
    alpha = -0.15*eye(numLinksIn);
    Qsent = alpha*(abs(Vmes(i,:)')-1);
    % Saturation
    
    % Controller to asset communication
    Qref(i,:) = cl.sampleOut(i,Qsent);
end

%% Plotting
