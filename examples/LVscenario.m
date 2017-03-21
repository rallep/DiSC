% Script for simulating the low voltage scenario
clc; close all; clear all;

% General setup
Ts = 60;                        % Sampling time
N = floor((1439*60/Ts)*1);      % Number of iterations
day = 160;                        % Julian day of the year (1-365)
onPU = false;                   % To indicate if the system shoul be simulated on a per unit basis
rng(1);                         % Set random generator

withCtrl = false;

%% Setup electrical grid
sBase = 50e6;
vBase = 20e3;
zBase = vBase^2/sBase;

% Grid impedance data
% format: [from to Resistance(R[ohm]) Reactance(X[ohm])  length(l[km])]
Z = [1 2    inf     inf     1;
     2 3    0.5     0.5     10;
     2 6    0.13    0.09     5;
     3 4    0.5     0.5       7;
     3 6    0.9     0.6     10;
     4 5    0.5     0.5     7;
     5 8    1.5      0.9       15;
     6 7    0.9     1     5;
     7 8    0.5     0.5    10;
     5 9    inf     inf     1; 
     9 10    0.208   0.052   0.5;
     10 11    0.208   0.052   0.5;
     11 12    0.208   0.052   0.5;
     12 13    0.208   0.052   0.5;
     13 14    0.208   0.052   0.5;
     14 15    0.208   0.052   0.5;
     ];

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
param.vBase = vBase;
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
numBus = length(Y);
% Type for each bus
type = [0 ones(1,numBus-1)]; % 0 = slack, 1 = PQ
% Voltages in per unit
Vin = [60e3 zeros(1,numBus-1)];
Vbase = [60e3 20e3*ones(1,7) 400*ones(1,7)];
% Avtive and reactive power at each bus
pfIndu = 0.9;       % Power factor industry
pfAgri = 0.9;       % Power factor agriculture
pfResi = 0.97;      % Power factor Residential
pfComm = 0.95;      % Power factor commercial

Pbus1 = zeros(N,1);
Pbus2 = zeros(N,1);
Pbus3 = -sum(HouseData.Data.pTs(1:N,1:200),2);      % Residential
Pbus4 = -4*commData.pTs(1:N);                       % Commercial
Pbus5 = zeros(N,1);                                 % Residential (LV grid)
Pbus6 = zeros(N,1);
Pbus7 = -2*induData.pTs(1:N);                       % Industry
Pbus8 = -5*agriData.pTs(1:N);                       % Agriculture
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
% Solar irradiance
param.lat = 56.889;     % Latitude for Sørup (degrees)
param.t = 0.75;         % Transmittance (unitless)
param.S = 1367;         % Solar constant (w/m^2)
param.p = 100;          % Air pressure (Kpa)
param.Ts = Ts;

MV_si = solarIrradiance(param);

% PV
MV_pvBus = 6;
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 6e6;
param.sMax = 6e6;
param.eta = 0.25;
param.A = 13000;
param.onPU = onPU;

MV_pv = pvAsset(param);

% Low voltage PV
% Solar irradiance
LV_si_Bus10 = solarIrradiance(param);
LV_si_Bus11 = solarIrradiance(param);
LV_si_Bus14 = solarIrradiance(param);
LV_si_Bus15 = solarIrradiance(param);
% PV
LV_pv1Bus = 10;
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 2*6e3;
param.sMax = 2*6e3;
param.eta = 0.25;
param.A = 2*25;
param.onPU = onPU;

LV_pv1 = pvAsset(param);

LV_pv2Bus = 11;
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 2*6e3;
param.sMax = 2*6e3;
param.eta = 0.25;
param.A = 2*25;
param.onPU = onPU;

LV_pv2 = pvAsset(param);

LV_pv3Bus = 14;
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 3*6e3;
param.sMax = 3*6e3;
param.eta = 0.25;
param.A = 3*25;
param.onPU = onPU;

LV_pv3= pvAsset(param);

LV_pv4Bus = 15;
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 5*6e3;
param.sMax = 5*6e3;
param.eta = 0.25;
param.A = 5*25;
param.onPU = onPU;

LV_pv4= pvAsset(param);

% Low voltage Energy Storage
LV_es1Bus = 12;
LV_es2Bus = 15;
param.sBase = sBase;
param.vBase = 400;
param.pRatedMax = 10e3;
param.pRatedMin = -10e3;
param.pRate = 50;
param.eRated = 100e3*60*60; % 100 kWh
param.Ts = Ts;
param.onPU = onPU;

LV_es1 = esAsset(param);
LV_es2 = esAsset(param);

% Low Voltage EV
param.sBase = sBase;
param.vBase = vBase;
param.pRated = 4e3;
param.eRated = 65e3*60*60;  % 65 kWh.
param.onPU = false;
param.Ts = Ts;

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



%% Setup simulation
% Allocate Memory
vOut = ones(N+1,numBus);    % Voltages    
pSlack = zeros(N,1);        % Active power at slack bus
qSlack = zeros(N,1);        % Reactive power at slack bus
pWpp = zeros(N,1);          % Available WT power
w = zeros(N,1);             % Wind speed
tapPos = zeros(N,1);
MV_G = zeros(N,1);
LV_es1_e = zeros(N,1);
LV_es2_e = zeros(N,1);
pEs1 = zeros(N,1);
qEs1 = zeros(N,1);
pEs2 = zeros(N,1);
qEs2 = zeros(N,1);

% Input to Assets
% Medium voltage
MV_cc = 0.00;
% Solar PV power plant
MV_pv_Vref = vBase;
MV_pv_dP = 0;          % Change in power reference
MV_pv_dPlim = 0;       % Power curtailment
MV_pv_qRef = 0;        % Reactive power reference
MV_pv.setQmode(0);    % Set reactive control mode 

% Transformer
MV_trafo_vRef = vBase;
MV_trafo_uRef = 0;
MV_trafo.setMode(1);
tapOld = 0;
MV_trafo.setTapSpec(10,10,10,5*60,0.0125,0.05);  

% Low voltage
% Transformer
LV_trafo_vRef = 400;
LV_trafo_uRef = 0;
LV_trafo.setMode(1);
LV_trafo.setTapSpec(10,10,10,5*60,0.0125,0.04);

LV_cc = 0.00;
% Solar PV power plant
LV_pv1_Vref = 400;
LV_pv1_dP = 0;          % Change in power reference
LV_pv1_dPlim = 0;       % Power curtailment
LV_pv1_qRef = 0;        % Reactive power reference
LV_pv1.setQmode(0);    % Set reactive control mode 

LV_pv2_Vref = 400;
LV_pv2_dP = 0;          % Change in power reference
LV_pv2_dPlim = 0;       % Power curtailment
LV_pv2_qRef = 0;        % Reactive power reference
LV_pv2.setQmode(0);    % Set reactive control mode 

LV_pv3_Vref = 400;
LV_pv3_dP = 0;          % Change in power reference
LV_pv3_dPlim = 0;       % Power curtailment
LV_pv3_qRef = 0;        % Reactive power reference
LV_pv3.setQmode(0);    % Set reactive control mode 

LV_pv4_Vref = 400;
LV_pv4_dP = 0;          % Change in power reference
LV_pv4_dPlim = 0;       % Power curtailment
LV_pv4_qRef = 0;        % Reactive power reference
LV_pv4.setQmode(0);    % Set reactive control mode 

% Energy storage
LV_es1.setPF(0.97);
LV_es2.setPF(0.97);
LV_es1.setDrain(1);
LV_es2.setDrain(1);
if withCtrl == true
    LV_es1.setPmode(1);
    LV_es2.setPmode(1);
else
    LV_es1.setPmode(0);
    LV_es2.setPmode(0);
end
LV_es1_vRef = 400;
LV_es2_vRef = 400;
LV_es1_pRef = 0;
LV_es2_pRef = 0;

% EV
if withCtrl == true
    EV_pRef = param.pRated/3;
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
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert assets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
    else
        % Tap changing transformers
        [y1,y2,y3] = MV_trafo.sample(abs(vOut(i-1,MV_trafoBusTo)),MV_trafo_vRef,MV_trafo_uRef,i,day);
        Y(MV_trafoBusFrom,MV_trafoBusFrom) = Yorg(MV_trafoBusFrom,MV_trafoBusFrom) + y1+y2;
        Y(MV_trafoBusTo,MV_trafoBusTo) = Yorg(MV_trafoBusTo,MV_trafoBusTo) +y1+y3;
        Y(MV_trafoBusFrom,MV_trafoBusTo) = -y1;
        Y(MV_trafoBusTo,MV_trafoBusFrom) = -y1;
        
        tapPos(i) = MV_trafo.tapPos;
        if MV_trafo.tapPos ~= tapOld
            disp('Tap has changed') 
            disp('Tap Position: ')
            MV_trafo.tapPos
            disp('Iteration')
            i
            disp('Number of taps pr day: ')
            MV_trafo.tapChDay
            tapOld = MV_trafo.tapPos;
        end
        
        [y1,y2,y3] = LV_trafo.sample(abs(vOut(i-1,9)),LV_trafo_vRef,LV_trafo_uRef,i,day);
        Y(LV_trafoBusFrom,LV_trafoBusFrom) = Yorg(LV_trafoBusFrom,LV_trafoBusFrom) + y1+y2;
        Y(LV_trafoBusTo,LV_trafoBusTo) = Yorg(LV_trafoBusTo,LV_trafoBusTo) +y1+y3;
        Y(LV_trafoBusFrom,LV_trafoBusTo) = -y1;
        Y(LV_trafoBusTo,LV_trafoBusFrom) = -y1;
        
        LV_tapPos(i) = LV_trafo.tapPos;
        
        % Solar PV power plant
        % Medium voltage
        MV_G(i) = MV_si.sample(i,day,MV_cc);
        [p,q] = MV_pv.sample(MV_G(i),abs(vOut(i-1,MV_pvBus)),MV_pv_dP,MV_pv_dPlim,MV_pv_qRef,MV_pv_Vref);
        Pin(i,MV_pvBus) = p;
        Qin(i,MV_pvBus) = q;
        
        % Low voltage
        % PV
        LV_G = LV_si_Bus10.sample(i,day,LV_cc);
        [p,q] = LV_pv1.sample(LV_G,abs(vOut(i-1,LV_pv1Bus)),LV_pv1_dP,LV_pv1_dPlim,LV_pv1_qRef,LV_pv1_Vref);
        Pin(i,LV_pv1Bus) = Pin(i,LV_pv1Bus)+p;
        Qin(i,LV_pv1Bus) = Qin(i,LV_pv1Bus)+q;
      
        LV_G = LV_si_Bus11.sample(i,day,LV_cc);
        [p,q] = LV_pv2.sample(LV_G,abs(vOut(i-1,LV_pv2Bus)),LV_pv2_dP,LV_pv2_dPlim,LV_pv2_qRef,LV_pv2_Vref);
        Pin(i,LV_pv2Bus) = Pin(i,LV_pv2Bus)+p;
        Qin(i,LV_pv2Bus) = Qin(i,LV_pv2Bus)+q;
        
        LV_G = LV_si_Bus14.sample(i,day,LV_cc);
        [p,q] = LV_pv3.sample(LV_G,abs(vOut(i-1,LV_pv3Bus)),LV_pv3_dP,LV_pv3_dPlim,LV_pv3_qRef,LV_pv3_Vref);
        Pin(i,LV_pv3Bus) = Pin(i,LV_pv3Bus)+p;
        Qin(i,LV_pv3Bus) = Qin(i,LV_pv3Bus)+q;
        
        LV_G = LV_si_Bus15.sample(i,day,LV_cc);
        [p,q] = LV_pv4.sample(LV_G,abs(vOut(i-1,LV_pv4Bus)),LV_pv4_dP,LV_pv4_dPlim,LV_pv4_qRef,LV_pv4_Vref);
        Pin(i,LV_pv4Bus) = Pin(i,LV_pv4Bus)+p;
        Qin(i,LV_pv4Bus) = Qin(i,LV_pv4Bus)+q;
        
        % ES
        [pEs1(i),qEs1(i),LV_es1_e(i)] = LV_es1.sample(abs(vOut(i-1,LV_es1Bus)),LV_es1_pRef,0,LV_es1_vRef);
        Pin(i,LV_es1Bus) = Pin(i,LV_es1Bus)+pEs1(i);
        Qin(i,LV_es1Bus) = Qin(i,LV_es1Bus)+qEs1(i);
        
        [pEs2(i),qEs2(i),LV_es2_e(i)] = LV_es2.sample(abs(vOut(i-1,LV_es2Bus)),LV_es2_pRef,0,LV_es2_vRef);
        Pin(i,LV_es2Bus) = Pin(i,LV_es2Bus)+pEs2(i);
        Qin(i,LV_es2Bus) = Qin(i,LV_es2Bus)+qEs2(i);
        
        % EV
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

% Voltages during a day without control
if withCtrl==false
    figure
    plot(9:15,abs(vOut(735,9:end))/400,'-o',9:15,abs(vOut(1380,9:end))/400,'-s',9:15,ones(1,7)*1.1,'--r',9:15,ones(1,7)*0.9,'--r')
    box off
    xlabel('Bus Number [-]')
    ylabel('Voltage [PU]')
    l=legend('Noon','Evening','Constraints','Orientation','horizontal','Location',[0.44 0.85 0.2 0.2]);
    legend('boxoff') 
    set(l,'Interpreter','latex','Fontsize',8)
    xlim([9 15])
    ylim([0.85 1.13])
    set(gca,'XTick',[9,10,11,12,13,14,15])
    set(gca,'YTick',[0.85,0.9,0.95,1,1.05,1.1,1.15])
    set(gca,'YTickLabel',{'0.85','0.90','0.95','1.00','1.05','1.10','1.15'})

    SAVE_PATH = [cd '/figures'];
    SAVE_MF = true;
    SIZE = 'paper';
    my_save_matlabfrag_papers('LVscenario_noCtrl',SAVE_MF,SIZE,SAVE_PATH);
else
    figure
    plot(9:15,abs(vOut(735,9:end))/400,'-o',9:15,abs(vOut(1380,9:end))/400,'-s',9:15,ones(1,7)*1.1,'--r',9:15,ones(1,7)*0.9,'--r')
    box off
    xlabel('Bus Number [-]')
    ylabel('Voltage [PU]')
    l=legend('Noon','Evening','Constraints','Orientation','horizontal','Location',[0.44 0.85 0.2 0.2]);
    legend('boxoff') 
    set(l,'Interpreter','latex','Fontsize',8)
    xlim([9 15])
    ylim([0.88 1.12])
    set(gca,'XTick',[9,10,11,12,13,14,15])
    set(gca,'YTick',[0.85,0.9,0.95,1,1.05,1.1,1.15])
    set(gca,'YTickLabel',{'0.85','0.90','0.95','1.00','1.05','1.10','1.15'})

    SAVE_PATH = [cd '/figures'];
    SAVE_MF = true;
    SIZE = 'paper';
    my_save_matlabfrag_papers('LVscenario_Ctrl',SAVE_MF,SIZE,SAVE_PATH);
end

figure
subplot(2,1,1)
plot(tvec,LV_es1_e,tvec,LV_es2_e)
subplot(2,1,2)
plot(tvec,pEs1,tvec,qEs1,tvec,pEs2,tvec,qEs2)
