% Script for simulating the low voltage scenario
clc; close all; clear all;

% General setup
Ts = 60;                        % Sampling time
N = floor((1439*60/Ts)*1);      % Number of iterations
day = 160;                        % Julian day of the year (1-365)
onPU = false;                   % To indicate if the system shoul be simulated on a per unit basis
rng(1);                         % Set random generator

withCtrl = true;

%% Setup electrical grid
sBase = 50e6;
vBase = 20e3;
zBase = vBase^2/sBase;

% Grid impedance data
% format: [from to Resistance(R[ohm]) Reactance(X[ohm])  length(l[km])]
Z = [1 2        inf     inf     1;
     2 3        0.5     0.5     10;
     2 6        0.13    0.09    5;
     3 4        0.5     0.5     7;
     3 6        0.9     0.6     10;
     4 5        0.5     0.5     7;
     5 8        1.5     0.9     15;
     6 7        0.9     1       5;
     7 8        0.5     0.5     10;
     5 9        inf     inf     1; 
     9 10       0.208   0.052   0.5;
     10 11      0.208   0.052   0.5;
     11 12      0.208   0.052   0.5;
     12 13      0.208   0.052   0.5;
     13 14      0.208   0.052   0.5;
     14 15      0.208   0.052   0.5;];
     
% Setup Newton-Raphson load flow solver
numBus = max(max(Z(:,2)));
param.type = [0 ones(1,numBus-1)]; % 0 = slack, 1 = PQ
param.vBase = [60e3 20e3*ones(1,7) 400*ones(1,7)];
pFlow = powerFlow(param);
% Construct grid admittance matrix
Y = pFlow.lDataToY(Z);
Yorg = Y;   % Make copy for trafo and switch implementation
 
% Setup transformers on MV and LV
% Medium voltage (Between bus 1 and 2)
param.zBase = zBase;
param.vBase = vBase;
param.nomTapRatio = 3;
param.Z = 2+1i*2;
param.Ts = Ts;
param.onPU = onPU;

MV_trafoBusFrom = 1;
MV_trafoBusTo = 2;

MV_trafo = tcAsset(param);

% Input trafo into admittance matrix
[y1,y2,y3] = MV_trafo.sample(1,1,0,1,1);
Y(MV_trafoBusFrom,MV_trafoBusFrom) = Yorg(MV_trafoBusFrom,MV_trafoBusFrom) + y1+y2;
Y(MV_trafoBusTo,MV_trafoBusTo) = Yorg(MV_trafoBusTo,MV_trafoBusTo) +y1+y3;
Y(MV_trafoBusFrom,MV_trafoBusTo) = -y1;
Y(MV_trafoBusTo,MV_trafoBusFrom) = -y1;

% Low voltage
% Low voltage (Between bus 5 and 9)
param.zBase = zBase;
param.vBase = 400;
param.onPU = onPU;
param.nomTapRatio = 50;
param.Z = 0.04+1i*0.04;
LV_trafoBusFrom = 5;
LV_trafoBusTo = 9;

LV_trafo = tcAsset(param);

% Input OLTC into admittance matrix first time
[y1,y2,y3] = LV_trafo.sample(1,1,0,1,1);
Y(LV_trafoBusFrom,LV_trafoBusFrom) = Yorg(LV_trafoBusFrom,LV_trafoBusFrom) + y1+y2;
Y(LV_trafoBusTo,LV_trafoBusTo) = Yorg(LV_trafoBusTo,LV_trafoBusTo) +y1+y3;
Y(LV_trafoBusFrom,LV_trafoBusTo) = -y1;
Y(LV_trafoBusTo,LV_trafoBusFrom) = -y1;

% Load consumption data and sample according to sampling time Ts
% Low voltage consumption profiles
HouseData = load('data/house1to900days31startSample1.mat');  % Mat file containing consumption data
% Interpolate 15 min. consumption data to match sampling
t = 0:length(HouseData.Data.HouseP)-1;
ti = 0:(length(HouseData.Data.HouseP)/(60*60*(1/Ts)*length(HouseData.Data.HouseP)/4)):length(HouseData.Data.HouseP)-1;
HouseData.Data.pTs = interp1(t,HouseData.Data.HouseP(:,:),ti);
HouseData.Data.pTs = HouseData.Data.pTs.*1000;

% Industry
induData = load('data/induPowerWinter');
t = 0:length(induData.p)-1;
ti = 0:(length(induData.p)/(60*60*(1/Ts)*length(induData.p))):length(induData.p)-1;
induData.pTs = interp1(t,induData.p(:,:),ti)';

% Agriculture
agriData = load('data/agriPowerWinter');
t = 0:length(agriData.p)-1;
ti = 0:(length(agriData.p)/(60*60*(1/Ts)*length(agriData.p))):length(agriData.p)-1;
agriData.pTs = interp1(t,agriData.p(:,:),ti)';

% Commercial
commData = load('data/commPowerWinter');
t = 0:length(commData.p)-1;
ti = 0:(length(commData.p)/(60*60*(1/Ts)*length(commData.p))):length(commData.p)-1;
commData.pTs = interp1(t,commData.p(:,:),ti)';

% Setup the busses
% Avtive and reactive power at each bus
pfIndu = 0.9;       % Power factor industry
pfAgri = 0.9;       % Power factor agriculture
pfResi = 0.97;      % Power factor Residential
pfComm = 0.95;      % Power factor commercial

Pbus1 = zeros(N,1);
Pbus2 = zeros(N,1);
Pbus3 = -sum(HouseData.Data.pTs(1:N,1:200),2);      % Residential
Pbus4 = -2*commData.pTs(1:N);                       % Commercial
Pbus5 = zeros(N,1);                                 % Residential (LV grid)
Pbus6 = zeros(N,1);
Pbus7 = -1*induData.pTs(1:N);                       % Industry
Pbus8 = -3*agriData.pTs(1:N);                       % Agriculture
Pbus9 = -2*sum(HouseData.Data.pTs(1:N,201:220),2);
Pbus10 = -sum(HouseData.Data.pTs(1:N,221:240),2);
Pbus11 = -sum(HouseData.Data.pTs(1:N,241:260),2);
Pbus12 = -sum(HouseData.Data.pTs(1:N,261:280),2);
Pbus13 = -sum(HouseData.Data.pTs(1:N,281:300),2);
Pbus14 = -sum(HouseData.Data.pTs(1:N,301:320),2);
Pbus15 = -sum(HouseData.Data.pTs(1:N,321:340),2);

Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7 Pbus8 Pbus9 Pbus10 ...
        Pbus11 Pbus12 Pbus13 Pbus14 Pbus15];

Qbus1 = zeros(N,1);
Qbus2 = zeros(N,1);
Qbus3 = Pbus3*tan(acos(pfResi));
Qbus4 = Pbus4*tan(acos(pfComm));
Qbus5 = zeros(N,1);
Qbus6 = zeros(N,1);
Qbus7 = Pbus7*tan(acos(pfIndu));
Qbus8 = Pbus8*tan(acos(pfAgri));
Qbus9 = Pbus9*tan(acos(pfResi));
Qbus10 = Pbus10*tan(acos(pfResi));
Qbus11 = Pbus11*tan(acos(pfResi));
Qbus12 = Pbus12*tan(acos(pfResi));
Qbus13 = Pbus13*tan(acos(pfResi));
Qbus14 = Pbus14*tan(acos(pfResi));
Qbus15 = Pbus15*tan(acos(pfResi));

Qin = [Qbus1 Qbus2 Qbus3 Qbus4 Qbus5 Qbus6 Qbus7 Qbus8 Qbus9 Qbus10 ...
        Qbus11 Qbus12 Qbus13 Qbus14 Qbus15];

%% Setup assets
% Medium voltage PV power plants
% % Wind speed
param.Ts = Ts;
param.z = 100;

MV_windSpeed = windSpeed(param);

% WT
% Wind power plant
param.Ts = Ts;
param.z = 100;
% WPP
MV_windTurbineBus = 6;
% Set parameters
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 2200e3;
param.sMax = 2200e3;
param.wMin = 3;
param.wMax = 25;
param.wRated = 12;
param.Ts = Ts;
param.onPU = onPU;
param.numWt = 10;
param.z = 100;
wpp1Prated = param.pRated*param.numWt;
% Create wind power plant object
MV_WT1 = wppAsset(param);

% MV_windTurbineBus = 6;
% param.sBase = sBase;
% param.vBase = vBase;
% param.pRated = 22e6;
% param.sMax = 22e6;
% param.wMin = 3;
% param.wMax = 25;
% param.wRated = 12;
% param.Ts = Ts;
% param.onPU =onPU;
% 
% MV_WT1 = wtAsset(param);

% Low voltage EVs
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 4e3;
param.eRated = 65e3*60*60;  % 65 kWh.
param.onPU = false;
param.Ts = Ts;
param.pRate = 100;

LV_ev1Bus = 9;
LV_ev1 = evAsset(param);
LV_ev2Bus = 13;
LV_ev2 = evAsset(param);
LV_ev3Bus = 13;
LV_ev3 = evAsset(param);
LV_ev4Bus = 13;
LV_ev4 = evAsset(param);
LV_ev5Bus = 14;
LV_ev5 = evAsset(param);
LV_ev6Bus = 14;
LV_ev6 = evAsset(param);
LV_ev7Bus = 14;
LV_ev7 = evAsset(param);
LV_ev8Bus = 14;
LV_ev8 = evAsset(param);
LV_ev9Bus = 15;
LV_ev9 = evAsset(param);
LV_ev10Bus = 15;
LV_ev10 = evAsset(param);
LV_ev11Bus = 15;
LV_ev11 = evAsset(param);
LV_ev12Bus = 15;
LV_ev12 = evAsset(param);
LV_ev13Bus = 15;
LV_ev13 = evAsset(param);

% Setup simulation
% Allocate Memory
vOut = ones(N+1,numBus);    % Voltages    
pSlack = zeros(N,1);        % Active power at slack bus
qSlack = zeros(N,1);        % Reactive power at slack bus
pWpp = zeros(N,1);          % Available WT power
w = zeros(N,1);             % Wind speed
LV_tapPos = zeros(N,1);
MV_tapPos = zeros(N,1);
MV_G = zeros(N,1);
pWpp = zeros(N,1);          % Available WT power
w = zeros(N,1);             % Wind speed

% Input to Assets
% Medium voltage
% Mean wind speed
wMean = 14;                                                                 % Wind Speed

% Wind turbine asset
wt_vRef = vBase;            % Voltage reference
wt_dP = 0;                  % Curtail power
if withCtrl == true
    wt_dPlim = 0;             % Derate power
else
    wt_dPlim = 0;
end
wt_qRef = 0;                % Reactive power reference
MV_WT1.setQmode(0);         % Set reactive control mode

% Transformer
MV_trafo_vRef = vBase;
MV_trafo_uRef = 0;
MV_trafo.setMode(1);
tapOld = 0;
MV_trafo.setTapSpec(15,15,50,5*60,0.0125,0.06);  

% Low voltage
% Transformer
LV_trafo_vRef = 400;
LV_trafo_uRef = 0;
LV_trafo.setMode(1);
LV_trafo.setTapSpec(15,15,50,5*60,0.0125,0.04);

LV_cc = 0;

% EVs
if withCtrl == true
    EV_pRef = param.pRated/2;
else
    EV_pRef = param.pRated;
end
LV_ev1_pRef = EV_pRef;
LV_ev1.setPmode(1);

LV_ev2_pRef = EV_pRef;
LV_ev2.setPmode(1);

LV_ev3_pRef = EV_pRef;
LV_ev3.setPmode(1);

LV_ev4_pRef = EV_pRef;
LV_ev4.setPmode(1);

LV_ev5_pRef = EV_pRef;
LV_ev5.setPmode(1);

LV_ev6_pRef = EV_pRef;
LV_ev6.setPmode(1);

LV_ev7_pRef = EV_pRef;
LV_ev7.setPmode(1);

LV_ev8_pRef = EV_pRef;
LV_ev8.setPmode(1);

LV_ev9_pRef = EV_pRef;
LV_ev9.setPmode(1);

LV_ev10_pRef = EV_pRef;
LV_ev10.setPmode(1);

LV_ev11_pRef = EV_pRef;
LV_ev11.setPmode(1);

LV_ev12_pRef = EV_pRef;
LV_ev12.setPmode(1);

LV_ev13_pRef = EV_pRef;
LV_ev13.setPmode(1);

temp=9;
%% Run simulation
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
    wMean = temp+3*sin(i*0.005); 
    
    if i==1260
        temp=8;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert assets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
    else
        % Tap changing transformers
        [y1,y2,y3] = MV_trafo.sample(abs(vOut(i-1,MV_windTurbineBus)),MV_trafo_vRef,MV_trafo_uRef,i,day);
        Y(MV_trafoBusFrom,MV_trafoBusFrom) = Yorg(MV_trafoBusFrom,MV_trafoBusFrom) + y1+y2;
        Y(MV_trafoBusTo,MV_trafoBusTo) = Yorg(MV_trafoBusTo,MV_trafoBusTo) +y1+y3;
        Y(MV_trafoBusFrom,MV_trafoBusTo) = -y1;
        Y(MV_trafoBusTo,MV_trafoBusFrom) = -y1;
        
        MV_tapPos(i) = MV_trafo.tapPos;
        
        [y1,y2,y3] = LV_trafo.sample(abs(vOut(i-1,11)),LV_trafo_vRef,LV_trafo_uRef,i,day);
        Y(LV_trafoBusFrom,LV_trafoBusFrom) = Yorg(LV_trafoBusFrom,LV_trafoBusFrom) + y1+y2;
        Y(LV_trafoBusTo,LV_trafoBusTo) = Yorg(LV_trafoBusTo,LV_trafoBusTo) +y1+y3;
        Y(LV_trafoBusFrom,LV_trafoBusTo) = -y1;
        Y(LV_trafoBusTo,LV_trafoBusFrom) = -y1;
        
        LV_tapPos(i) = LV_trafo.tapPos;
        if LV_trafo.tapPos ~= tapOld
            disp('Tap has changed') 
            disp('Tap Position: ')
            LV_trafo.tapPos
            disp('Iteration')
            i
            disp('Number of taps pr day: ')
            LV_trafo.tapChDay
            tapOld = LV_trafo.tapPos;
        end
        % Medium Voltage
        % Wind power plant
%         [w(i), we] = MV_windSpeed.sample(i,wMean);
        [p,q,pWpp(i)] = MV_WT1.sample(wMean,abs(vOut(i-1,MV_windTurbineBus)),wt_dP,wt_dPlim,wt_qRef,wt_vRef,i);
        Pin(i,MV_windTurbineBus) = p;
        Qin(i,MV_windTurbineBus) = q;
        
        % EVs
        [pEv1,qEv1,e,away] = LV_ev1.sample(i,day,400,LV_ev1_pRef,0,400);
        Pin(i,LV_ev1Bus) = Pin(i,LV_ev1Bus)+pEv1;
        Qin(i,LV_ev1Bus) = Qin(i,LV_ev1Bus)+qEv1;
        
        [pEv2,qEv2,e,away] = LV_ev2.sample(i,day,400,LV_ev2_pRef,0,400);
        
        [pEv3,qEv3,e,away] = LV_ev3.sample(i,day,400,LV_ev3_pRef,0,400);
        
        [pEv4,qEv4,e,away] = LV_ev4.sample(i,day,400,LV_ev4_pRef,0,400);
        Pin(i,LV_ev4Bus) = Pin(i,LV_ev4Bus) + pEv2+pEv3+pEv4;
        Qin(i,LV_ev4Bus) = Qin(i,LV_ev4Bus) + qEv2+qEv3+qEv4;
        
        [pEv5,qEv5,e,away] = LV_ev5.sample(i,day,400,LV_ev5_pRef,0,400);
        
        [pEv6,qEv6,e,away] = LV_ev6.sample(i,day,400,LV_ev6_pRef,0,400);
        
        [pEv7,qEv7,e,away] = LV_ev7.sample(i,day,400,LV_ev7_pRef,0,400);
        
        [pEv8,qEv8,e,away] = LV_ev8.sample(i,day,400,LV_ev8_pRef,0,400);
        Pin(i,LV_ev8Bus) = Pin(i,LV_ev8Bus)+pEv5+pEv6+pEv7+pEv8;
        Qin(i,LV_ev8Bus) = Qin(i,LV_ev8Bus)+qEv5+qEv6+qEv7+qEv8;
        
        [pEv9,qEv9,e,away] = LV_ev9.sample(i,day,400,LV_ev9_pRef,0,400);
        
        [pEv10,qEv10,e,away] = LV_ev10.sample(i,day,400,LV_ev10_pRef,0,400);
        
        [pEv11,qEv11,e,away] = LV_ev11.sample(i,day,400,LV_ev11_pRef,0,400);
        
        [pEv12,qEv12,e,away] = LV_ev12.sample(i,day,400,LV_ev12_pRef,0,400);

        [pEv13,qEv13,e,away] = LV_ev13.sample(i,day,400,LV_ev13_pRef,0,400);
        Pin(i,LV_ev13Bus) = Pin(i,LV_ev13Bus)+pEv9+pEv10+pEv11+pEv12+pEv13;
        Qin(i,LV_ev13Bus) = Qin(i,LV_ev13Bus)+qEv9+qEv10+qEv11+qEv12+qEv13;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Electrical Grid  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [vOut(i,:),pSlack(i),qSlack(i),~]=pFlow.nrLoadFlow(Y,Pin(i,:)',Qin(i,:)');
    if pFlow.nIte == 100
        disp('ERROR nrIte = 100')
        break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if withCtrl == true
        if i > 0 && i < 660
            wt_dPlim = 8e6;
        elseif i >= 660 && i < 1320
            wt_dPlim = 0;
        elseif i >= 1320
            wt_dPlim = 0;
        end
    end  
end
%% Plotting
close all;
tvec = (0:N-1)/60*60/Ts;

% P and Q slack
figure
plot(tvec,pSlack/1e3,tvec,qSlack/1e3)
xlabel('Time [1 min.]')
ylabel('Pslack / Qslack [kW / kVar]')
title('P/Q slack')

% Voltages
figure
plot(tvec,abs(vOut(1:N,2:8))/vBase)
title('Voltages MV')
xlabel('Time [hrs]')
ylabel('Volatges [PU]')

% Voltages
figure
plot(tvec,abs(vOut(1:N,9:end))/0.4e3)
title('Voltages LV')
xlabel('Time [hrs]')
ylabel('Volatges [PU]')
ylim([0.8 1.2])

% Wind turbine
figure
plot(tvec,pWpp/1e6,tvec,Pin(1:N,MV_windTurbineBus)/1e6,tvec,Qin(1:N,MV_windTurbineBus)/1e6)

% Without Control
if withCtrl == false
    figure
    [hAx,hLine1,hLine2] = plotyy(tvec,pWpp/1e6,tvec,[MV_tapPos LV_tapPos]);
    set(hLine2,'LineWidth',2)
    set(hLine1,'color','b')
    set(hLine1,'LineStyle','-')
    set(hLine2(1),'color',[0 0.5 0])
    set(hLine2(2),'color','r')

    set(hAx(1),'ycolor','k')
    set(hAx(1),'ylim',[0 23],'ytick',[0 5 10 15 20 25]);
    set(hAx(2),'ylim',[-6.1 8.1],'ytick',[-7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9]);
    set(hAx(1),'xlim',[0 24],'xtick',[0 6 12 18 24],'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
    set(hAx(2),'xlim',[0 24],'xtick',[])
    xlabel('Time of Day')
    ylabel(hAx(1),'Power [MW]') % left y-axis
    ylabel(hAx(2),'Tap Position [-]')
    ylabh = get(hAx(2),'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [-0.8 0 0])
    box off
    leg = legend('Produced','HV/MV Trafo.','MV/LV Trafo.','location',[0.4 0.905 0.1 0.1],'orientation','horizontal');
    legend('boxoff')
    set(leg,'Interpreter','latex','Fontsize',8)

    SAVE_PATH = [cd '/figures'];
    SAVE_MF = true;
    SIZE = [18 10];
    my_save_matlabfrag_papers('MVscenario_noCtrl',SAVE_MF,SIZE,SAVE_PATH);
end

% With control
if withCtrl == true
    figure
    [hAx,hLine1,hLine2] = plotyy(tvec,[pWpp/1e6 Pin(1:N,MV_windTurbineBus)/1e6],tvec,[MV_tapPos LV_tapPos]);
    set(hLine2,'LineWidth',2)
    set(hLine1(1),'LineStyle',':')
    set(hLine1(1),'color','b')
    set(hLine1(2),'color','b')
    set(hLine2(1),'color',[0 0.5 0])
    set(hLine2(2),'color','r')
    set(hAx(1),'ylim',[0 23],'ytick',[0 5 10 15 20 25]);
    set(hAx(2),'ylim',[-3.1 5.1],'ytick',[-4 -3 -2 -1 0 1 2 3 4 5]);
    set(hAx(1),'xlim',[0 24],'xtick',[0 6 12 18 24],'XTickLabel',{'0:00','6:00','12:00','18:00','24:00'})
    set(hAx(2),'xlim',[0 24],'xtick',[])
    xlabel('Time of Day')
    ylabel(hAx(1),'Power [MW]')
    ylabel(hAx(2),'Tap Position [-]')
    ylabh = get(hAx(2),'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [-0.8 0 0])
    box off
    leg = legend('Available','Produced','HV/MV Trafo.','MV/LV Trafo.','location',[0.4 0.905 0.1 0.1],'orientation','horizontal');
    legend('boxoff')
    set(leg,'Interpreter','latex','Fontsize',8)
    
    SAVE_PATH = [cd '/figures'];
    SAVE_MF = true;
    SIZE = [18 10];
    my_save_matlabfrag_papers('MVscenario_Ctrl',SAVE_MF,SIZE,SAVE_PATH);
end


