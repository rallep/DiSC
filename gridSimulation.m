% Script for simulating the entire electrical grid, with both MV and LV
% benchmark grids.
%
% 19-09-2014, R. Pedersen, Aalborg University
clc; clear all; close all;

%% General setup
Ts = 60;        % Sampling time
N = 1440-1;     % Number of itterations 
day = 160;        % Julian day of the year (1-365)
onPU = false;   % To indicate if the system should be simulated on a per unit basis
rng(2); 


%% Setup the electrical grid
gridSetup

%% Simulation
% Input for assets
% Medium voltage tap changers
if onPU == true
    MV_tc1vRef = 1;
    LV_tc1vRef = 1;
else
    MV_tc1vRef = MV_vBase;
    LV_tc1vRef = LV_vBase;
end
MV_tc1uRef = 0;
MV_tc1.setMode(1);

% MV_tc2vRef = 1;
% MV_tc2uRef = 0;

% Low voltage tap changers
LV_tc1uRef = 0;
LV_tc1.setMode(0);

% Supermarket
sm1_u = [0 0];

% Solar Irradiance
MV_cc = 0.05;  % Cloud cover

% Solar PV power plant
MV_pv1vRef = MV_vBase;
MV_pv1dP = 0;          % Change in power reference
MV_pv1dPlim = 0;       % Power curtailment
MV_pv1qRef = 0;        % Reactive power reference
MV_pv1.setQmode(0);    % Set reactive control mode 

% Wind speed
MV_wMean = 11;

% Wind power plant
MV_wt1vRef = MV_vBase;
MV_wt1dP = 0;          % Change in power reference
MV_wt1dPlim = 0;       % Power curtailment
MV_wt1qRef = 0;        % Reactive power reference
MV_wt1.setQmode(0);    % Set reactive control mode 

% Allocate memory
vOut = ones(N+1,numBus);
pSlack = zeros(N,1);
qSlack = zeros(N,1);
MV_tapPos = zeros(N,1);
w = zeros(N,1);
we = zeros(N,1);

% For testing
MV_tapOld = 0;
tic
for i=1:N
    % Itterate day
    if ~mod(i,24*60*60/Ts)
       if day == 365
           day = 1;
       else
           day = day +1;
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set invironmental data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MV_wMean = 9+3*sin(i*0.005);    % 8 mean div 3 gives more than 10
    %MV_wMean = 8;
    MV_cc = 1+1*sin(i*0.005);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert assets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
    else
        % Tap changing transformers
        % Medium voltage
        [y1,y2,y3] = MV_tc1.sample(abs(vOut(i-1,MV_wt1Bus)),MV_tc1vRef,MV_tc1uRef,i,day);
        Y(MV_tc1BusFrom,MV_tc1BusFrom) = Yorg(MV_tc1BusFrom,MV_tc1BusFrom) + y1+y2;
        Y(MV_tc1BusTo,MV_tc1BusTo) = Yorg(MV_tc1BusTo,MV_tc1BusTo) +y1+y3;
        Y(MV_tc1BusFrom,MV_tc1BusTo) = -y1;
        Y(MV_tc1BusTo,MV_tc1BusFrom) = -y1;
        MV_tapPos(i) = MV_tc1.tapPos;
        if MV_tc1.tapPos ~= MV_tapOld
            disp('Tap has changed') 
            disp('Tap Position: ')
            MV_tc1.tapPos
            disp('Iteration')
            i
            disp('Number of taps pr day: ')
            MV_tc1.tapChDay
            MV_tapOld = MV_tc1.tapPos;
        end
        
%         [y1,y2,y3] = MV_tc2.sample(abs(vOut(i-1,MV_tc2BusTo)),MV_tc2vRef,MV_tc2uRef,i,day);
%         Y(MV_tc2BusFrom,MV_tc2BusFrom) = Yorg(MV_tc2BusFrom,MV_tc2BusFrom) + y1+y2;
%         Y(MV_tc2BusTo,MV_tc2BusTo) = Yorg(MV_tc2BusTo,MV_tc2BusTo) +y1+y3;
%         Y(MV_tc2BusFrom,MV_tc2BusTo) = -y1;
%         Y(MV_tc2BusTo,MV_tc2BusFrom) = -y1;
        
        % Low voltage
        [y1,y2,y3] = LV_tc1.sample(abs(vOut(i-1,LV_tc1BusTo)),LV_tc1vRef,LV_tc1uRef,i,day);
        Y(LV_tc1BusFrom,LV_tc1BusFrom) = Yorg(LV_tc1BusFrom,LV_tc1BusFrom) + y1+y2;
        Y(LV_tc1BusTo,LV_tc1BusTo) = Yorg(LV_tc1BusTo,LV_tc1BusTo) +y1+y3;
        Y(LV_tc1BusFrom,LV_tc1BusTo) = -y1;
        Y(LV_tc1BusTo,LV_tc1BusFrom) = -y1;
        
        % Solar PV power plant
        % Medium voltage
        MV_G = MV_si.sample(i,day,MV_cc);
        [p,q] = MV_pv1.sample(MV_G,abs(vOut(i-1,MV_pv1Bus)),MV_pv1dP,MV_pv1dPlim,MV_pv1qRef,MV_pv1vRef);
        Pin(i,MV_pv1Bus) = p*0;
        Qin(i,MV_pv1Bus) = q*0;
        
        % Wind power plant
        % Medium voltage
        [w(i), we(i)] = ws.sample(i,MV_wMean);
        [p,q,~] = MV_wt1.sample(w(i),abs(vOut(i-1,MV_wt1Bus)),MV_wt1dP,MV_wt1dPlim,MV_wt1qRef,MV_wt1vRef);
        [w(i), we(i)] = ws2.sample(i,MV_wMean);
        [p2,q2,~] = MV_wt2.sample(w(i),abs(vOut(i-1,MV_wt1Bus)),MV_wt1dP,MV_wt1dPlim,MV_wt1qRef,MV_wt1vRef);
        Pin(i,MV_wt1Bus) = p+p2;
        Qin(i,MV_wt1Bus) = q+q2;
        
        % Supermarket
        % Medium voltage
        [p,q,sm1_timeReady,sm1_timeDefrost] = sm1.sample(sm1_u,smData.pTs(i));
        Pin(i,sm1Bus) = p;
        Qin(i,sm1Bus) = q;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Electrical Grid  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out=nrLoadFlow(Y,type,Pin(i,:)',Qin(i,:)',Vin',tol,maxIte);
    vOut(i,:) = out.Vout;       % Complex Voltage
    pSlack(i) = out.Pslack;
    qSlack(i) = out.Qslack;
    out.nrIte;
    if out.nrIte == 100
        disp('ERROR nrIte = 100')
        break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurement + Delay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Substation Control 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Send references + Delay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end
toc

%% Plotting of result
tvec = (0:N-1)/60*60/Ts;
% MV grid
figure
plot(tvec,pSlack/1e6,tvec,qSlack/1e6)
xlabel('Time [1 min.]')
ylabel('Pslack / Qslack [MW / MVar]')
title('P/Q slack')

figure
subplot(2,1,1)
plot(tvec,abs(vOut(1:N,2))/MV_vBase,tvec,abs(vOut(1:N,4:12))/MV_vBase)
title('Voltages MV')
xlabel('Time [hrs]')
ylabel('Volatges [PU]')

subplot(2,1,2)
plot(tvec,abs(vOut(1:N,13:53))/LV_vBase)
title('Voltages LV')
xlabel('Time [hrs]')
ylabel('Voltage [PU]')

figure
subplot(2,1,1)
plot(tvec,w,tvec,we)
ylabel('Wind Speed [m/s]')
subplot(2,1,2)
plot(tvec,Pin(1:N,5)/1e6,tvec,Pin(1:N,MV_wt1Bus)/1e6)
title('Power WT and PV')
xlabel('Time [hrs]')
ylabel('Power [MW]')

figure
plot(tvec,MV_tapPos)
