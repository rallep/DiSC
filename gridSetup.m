% Script for setting up the benchmark grid
% R.Pedersen 10-09-2014, Aalborg University

%% Setup the grid
% MV grid base parameters
MV_sBase = 50e6;                        % Complex power base [VAR]
MV_vBase = 20e3;                        % Voltage base [V]
MV_zBase = MV_vBase^2/MV_sBase;         % Impedance base [Ohm]

% LV grid
LV_sBase = 400e3;                       % Complex power base [VAR]
LV_vBase = 0.4e3;                       % Voltage base [V]
LV_zBase = LV_vBase^2/LV_sBase;         % Impedance base [Ohm]

% Load grid admittance matrix
Y = benchmarkAdmittanceMatrix('both',true,true,MV_zBase,LV_zBase,onPU);
Yorg = Y; % Used for tap changing transformers and switches

% Power facotrs
pfIndu = 0.9;       % Power factor industry
pfAgri = 0.9;       % Power factor agriculture
pfResi = 0.97;      % Power factor Residential
pfComm = 0.95;      % Power factor commercial

%% Setup Inflexible Consumption
% Medium voltage consumption profiles
% Supermarket
smData = load('mv_grid/consData/smPowerJan2013.mat');  % Mat file containing consumption data
t = 0:length(smData.power)-1;
ti = 0:(length(smData.power)/(60*60*(1/Ts)*length(smData.power)/60)):length(smData.power)-1;
smData.pTs = interp1(t,smData.power(:,:),ti);

% Industry
induData = load('mv_grid/consData/induPowerWinter');
t = 0:length(induData.p)-1;
ti = 0:(length(induData.p)/(60*60*(1/Ts)*length(induData.p))):length(induData.p)-1;
induData.pTs = interp1(t,induData.p(:,:),ti)';

% Agriculture
agriData = load('mv_grid/consData/agriPowerWinter');
t = 0:length(agriData.p)-1;
ti = 0:(length(agriData.p)/(60*60*(1/Ts)*length(agriData.p))):length(agriData.p)-1;
agriData.pTs = interp1(t,agriData.p(:,:),ti)';

% Commercial
commData = load('mv_grid/consData/commPowerWinter');
t = 0:length(commData.p)-1;
ti = 0:(length(commData.p)/(60*60*(1/Ts)*length(commData.p))):length(commData.p)-1;
commData.pTs = interp1(t,commData.p(:,:),ti)';

% Low voltage consumption profiles
HouseData = load('lv_grid/consData/house1to120days31startSample17280.mat');  % Mat file containing consumption data
% Interpolate 15 min. consumption data to match sampling
t = 0:length(HouseData.Data.HouseP)-1;
ti = 0:(length(HouseData.Data.HouseP)/(60*60*(1/Ts)*length(HouseData.Data.HouseP)/4)):length(HouseData.Data.HouseP)-1;
HouseData.Data.pTs = interp1(t,HouseData.Data.HouseP(:,:),ti);
HouseData.Data.pTs = HouseData.Data.pTs.*1000;

% Set number of samples
findMinSample = [length(HouseData.Data.pTs) length(smData.pTs) length(induData.pTs)];
numSamples = min(findMinSample);         % Maximum number of samples

%% Setup Busses MV = bus 1-12, LV = bus 13-53
numBus = length(Y);
% Type for each bus
type = [0 ones(1,numBus-1)]; % 0 = slack, 1 = PQ
% Voltages in per unit
if onPU == true
    Vin = [1 zeros(1,numBus-1)]; % Bus one is the slack bus
else
    Vin = [60e3 zeros(1,numBus-1)];
end

% Active power input at each bus
% Medium voltage grid
Pbus1 = zeros(numSamples,1);
Pbus2 = zeros(numSamples,1);
Pbus3 = zeros(numSamples,1);
Pbus4 = -5*induData.pTs(1:numSamples);      % Industry
Pbus5 = zeros(numSamples,1);                % PV plant (From object)
Pbus6 = -6*commData.pTs(1:numSamples);      % Commercial
Pbus7 = zeros(numSamples,1);        
Pbus8 = -6*agriData.pTs(1:numSamples);      % Agriculture
Pbus9 = zeros(numSamples,1);                % Supermarket (From object)
Pbus10 = zeros(numSamples,1);               % Residential (The low voltage grid)
Pbus11 = zeros(numSamples,1);
Pbus12 = zeros(numSamples,1);               % Wind turbine (From object)

% Low voltage grid
Pbus13 = zeros(numSamples,1);
Pbus14 = -sum(HouseData.Data.pTs(1:numSamples,1:8),2);
Pbus15 = -sum(HouseData.Data.pTs(1:numSamples,9:12),2);
Pbus16 = -sum(HouseData.Data.pTs(1:numSamples,13:15),2);
Pbus17 = -sum(HouseData.Data.pTs(1:numSamples,16:19),2);
Pbus18 = -sum(HouseData.Data.pTs(1:numSamples,20:22),2);
Pbus19 = -sum(HouseData.Data.pTs(1:numSamples,23:27),2);
Pbus20 = zeros(numSamples,1);
Pbus21 = -sum(HouseData.Data.pTs(1:numSamples,28:31),2);
Pbus22 = -sum(HouseData.Data.pTs(1:numSamples,32:35),2);
Pbus23 = -sum(HouseData.Data.pTs(1:numSamples,36:37),2);
Pbus24 = -sum(HouseData.Data.pTs(1:numSamples,38:38),2);
Pbus25 = -sum(HouseData.Data.pTs(1:numSamples,39:40),2);
Pbus26 = -sum(HouseData.Data.pTs(1:numSamples,41:45),2);
Pbus27 = -sum(HouseData.Data.pTs(1:numSamples,46:49),2);
Pbus28 = -sum(HouseData.Data.pTs(1:numSamples,50:52),2);
Pbus29 = -sum(HouseData.Data.pTs(1:numSamples,53:54),2);
Pbus30 = -sum(HouseData.Data.pTs(1:numSamples,55:58),2);
Pbus31 = -sum(HouseData.Data.pTs(1:numSamples,59:61),2);
Pbus32 = -sum(HouseData.Data.pTs(1:numSamples,62:63),2);
Pbus33 = -sum(HouseData.Data.pTs(1:numSamples,64:67),2);
Pbus34 = -sum(HouseData.Data.pTs(1:numSamples,68:70),2);
Pbus35 = -sum(HouseData.Data.pTs(1:numSamples,71:72),2);
Pbus36 = -sum(HouseData.Data.pTs(1:numSamples,73:75),2);
Pbus37 = -sum(HouseData.Data.pTs(1:numSamples,76:77),2);
Pbus38 = -sum(HouseData.Data.pTs(1:numSamples,78:80),2);
Pbus39 = -sum(HouseData.Data.pTs(1:numSamples,81:81),2);
Pbus40 = -sum(HouseData.Data.pTs(1:numSamples,82:83),2);
Pbus41 = -sum(HouseData.Data.pTs(1:numSamples,84:85),2);
Pbus42 = -sum(HouseData.Data.pTs(1:numSamples,86:87),2);
Pbus43 = -sum(HouseData.Data.pTs(1:numSamples,88:90),2);
Pbus44 = zeros(numSamples,1);
Pbus45 = -sum(HouseData.Data.pTs(1:numSamples,91:94),2);
Pbus46 = -sum(HouseData.Data.pTs(1:numSamples,95:97),2);
Pbus47 = -sum(HouseData.Data.pTs(1:numSamples,98:104),2);
Pbus48 = -sum(HouseData.Data.pTs(1:numSamples,105:106),2);
Pbus49 = -sum(HouseData.Data.pTs(1:numSamples,107:108),2);
Pbus50 = -sum(HouseData.Data.pTs(1:numSamples,109:109),2);
Pbus51 = -sum(HouseData.Data.pTs(1:numSamples,110:111),2);
Pbus52 = -sum(HouseData.Data.pTs(1:numSamples,112:113),2);
Pbus53 = -sum(HouseData.Data.pTs(1:numSamples,114:116),2);

% Collect into a matrix
if onPU == true
    Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7 Pbus8 Pbus9 Pbus10 Pbus11...
            Pbus12 Pbus13 Pbus14 Pbus15 Pbus16 Pbus17 Pbus18 Pbus19 Pbus20...
            Pbus21 Pbus22 Pbus23 Pbus24 Pbus25 Pbus26 Pbus27 Pbus28 Pbus29...
            Pbus30 Pbus31 Pbus32 Pbus33 Pbus34 Pbus35 Pbus36 Pbus37 Pbus38...
            Pbus39 Pbus40 Pbus41 Pbus42 Pbus43 Pbus44 Pbus45 Pbus46 Pbus47...
            Pbus48 Pbus49 Pbus50 Pbus51 Pbus52 Pbus53];
    Pin(:,1:12) = Pin(:,1:12)/MV_sBase;
    Pin(:,13:end) = Pin(:,13:end)/LV_sBase;
else
    Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7 Pbus8 Pbus9 Pbus10 Pbus11...
            Pbus12 Pbus13 Pbus14 Pbus15 Pbus16 Pbus17 Pbus18 Pbus19 Pbus20...
            Pbus21 Pbus22 Pbus23 Pbus24 Pbus25 Pbus26 Pbus27 Pbus28 Pbus29...
            Pbus30 Pbus31 Pbus32 Pbus33 Pbus34 Pbus35 Pbus36 Pbus37 Pbus38...
            Pbus39 Pbus40 Pbus41 Pbus42 Pbus43 Pbus44 Pbus45 Pbus46 Pbus47...
            Pbus48 Pbus49 Pbus50 Pbus51 Pbus52 Pbus53];
end

% Reactive power input at each bus
% Medium voltage
Qbus1 = zeros(numSamples,1);
Qbus2 = zeros(numSamples,1);
Qbus3 = zeros(numSamples,1);
Qbus4 = Pbus4*tan(acos(pfIndu));        % Industri
Qbus5 = zeros(numSamples,1);            % PV plant (From object)
Qbus6 = Pbus6*tan(acos(pfComm));        % Commercial
Qbus7 = zeros(numSamples,1);
Qbus8 = Pbus8*tan(acos(pfAgri));        % Aggriculture
Qbus9 = zeros(numSamples,1);            % Supermarket (From object)
Qbus10 = zeros(numSamples,1);           % Residential (low voltage grid)
Qbus11 = zeros(numSamples,1);
Qbus12 = zeros(numSamples,1);           % Wind farm (From object)

% Low voltage 
Qbus13 = zeros(numSamples,1);
Qbus14 = Pbus14*tan(acos(pfResi));
Qbus15 = Pbus15*tan(acos(pfResi));
Qbus16 = Pbus16*tan(acos(pfResi));
Qbus17 = Pbus17*tan(acos(pfResi));
Qbus18 = Pbus18*tan(acos(pfResi));
Qbus19 = Pbus19*tan(acos(pfResi));
Qbus20 = Pbus20*tan(acos(pfResi));
Qbus21 = Pbus21*tan(acos(pfResi));
Qbus22 = Pbus22*tan(acos(pfResi));
Qbus23 = Pbus23*tan(acos(pfResi));
Qbus24 = Pbus24*tan(acos(pfResi));
Qbus25 = Pbus25*tan(acos(pfResi));
Qbus26 = Pbus26*tan(acos(pfResi));
Qbus27 = Pbus27*tan(acos(pfResi));
Qbus28 = Pbus28*tan(acos(pfResi));
Qbus29 = Pbus29*tan(acos(pfResi));
Qbus30 = Pbus30*tan(acos(pfResi));
Qbus31 = Pbus31*tan(acos(pfResi));
Qbus32 = Pbus32*tan(acos(pfResi));
Qbus33 = Pbus33*tan(acos(pfResi));
Qbus34 = Pbus34*tan(acos(pfResi));
Qbus35 = Pbus35*tan(acos(pfResi));
Qbus36 = Pbus36*tan(acos(pfResi));
Qbus37 = Pbus37*tan(acos(pfResi));
Qbus38 = Pbus38*tan(acos(pfResi));
Qbus39 = Pbus39*tan(acos(pfResi));
Qbus40 = Pbus40*tan(acos(pfResi));
Qbus41 = Pbus41*tan(acos(pfResi));
Qbus42 = Pbus42*tan(acos(pfResi));
Qbus43 = Pbus43*tan(acos(pfResi));
Qbus44 = Pbus44*tan(acos(pfResi));
Qbus45 = Pbus45*tan(acos(pfResi));
Qbus46 = Pbus46*tan(acos(pfResi));
Qbus47 = Pbus47*tan(acos(pfResi));
Qbus48 = Pbus48*tan(acos(pfResi));
Qbus49 = Pbus49*tan(acos(pfResi));
Qbus50 = Pbus50*tan(acos(pfResi));
Qbus51 = Pbus51*tan(acos(pfResi));
Qbus52 = Pbus52*tan(acos(pfResi));
Qbus53 = Pbus53*tan(acos(pfResi));

% Collect into a matrix
if onPU == true
    Qin = [Qbus1 Qbus2 Qbus3 Qbus4 Qbus5 Qbus6 Qbus7 Qbus8 Qbus9 Qbus10 Qbus11...
            Qbus12 Qbus13 Qbus14 Qbus15 Qbus16 Qbus17 Qbus18 Qbus19 Qbus20...
            Qbus21 Qbus22 Qbus23 Qbus24 Qbus25 Qbus26 Qbus27 Qbus28 Qbus29...
            Qbus30 Qbus31 Qbus32 Qbus33 Qbus34 Qbus35 Qbus36 Qbus37 Qbus38...
            Qbus39 Qbus40 Qbus41 Qbus42 Qbus43 Qbus44 Qbus45 Qbus46 Qbus47...
            Qbus48 Qbus49 Qbus50 Qbus51 Qbus52 Qbus53];

    Qin(:,1:12) = Qin(:,1:12)/MV_sBase;
    Qin(:,13:end) = Qin(:,13:end)/LV_sBase;
else
    Qin = [Qbus1 Qbus2 Qbus3 Qbus4 Qbus5 Qbus6 Qbus7 Qbus8 Qbus9 Qbus10 Qbus11...
        Qbus12 Qbus13 Qbus14 Qbus15 Qbus16 Qbus17 Qbus18 Qbus19 Qbus20...
        Qbus21 Qbus22 Qbus23 Qbus24 Qbus25 Qbus26 Qbus27 Qbus28 Qbus29...
        Qbus30 Qbus31 Qbus32 Qbus33 Qbus34 Qbus35 Qbus36 Qbus37 Qbus38...
        Qbus39 Qbus40 Qbus41 Qbus42 Qbus43 Qbus44 Qbus45 Qbus46 Qbus47...
        Qbus48 Qbus49 Qbus50 Qbus51 Qbus52 Qbus53];
end
    
%% Setup Assets
% On-Load Tap Changing Transformers
% Medium voltage (Between bus 1 and 2, and between bus 1 and 3)
param.zBase = MV_zBase;
param.vBase = MV_vBase;
param.nomTapRatio = 3;
param.Z = 0.6-1i*0.5;
param.Ts = Ts;
param.onPU = onPU;

MV_tc1BusFrom = 1;
MV_tc1BusTo = 2;

MV_tc1 = tcAsset(param);
MV_tc1.setTapSpec(15,15,50,5*60,0.0125,0.05);

% MV_tc2BusFrom = 1;
% MV_tc2BusTo = 3;
% 
% MV_tc2 = tcAsset(param);

% Input OLTC into admittance matrix first time
[y1,y2,y3] = MV_tc1.sample(1,1,0,1,1);
Y(MV_tc1BusFrom,MV_tc1BusFrom) = Yorg(MV_tc1BusFrom,MV_tc1BusFrom) + y1+y2;
Y(MV_tc1BusTo,MV_tc1BusTo) = Yorg(MV_tc1BusTo,MV_tc1BusTo) +y1+y3;
Y(MV_tc1BusFrom,MV_tc1BusTo) = -y1;
Y(MV_tc1BusTo,MV_tc1BusFrom) = -y1;

% [y1,y2,y3] = MV_tc2.sample(1,1,0,1,1);
% Y(MV_tc2BusFrom,MV_tc2BusFrom) = Yorg(MV_tc2BusFrom,MV_tc2BusFrom) + y1+y2;
% Y(MV_tc2BusTo,MV_tc2BusTo) = Yorg(MV_tc2BusTo,MV_tc2BusTo) +y1+y3;
% Y(MV_tc2BusFrom,MV_tc2BusTo) = -y1;
% Y(MV_tc2BusTo,MV_tc2BusFrom) = -y1;

% Low voltage (Between bus 10 and 13)
param.zBase = LV_zBase;
param.vBase = LV_vBase;
param.onPU = onPU;
param.nomTapRatio = 50;
param.Z = 0.04-1i*0.04;
LV_tc1BusFrom = 10;
LV_tc1BusTo = 13;

LV_tc1 = tcAsset(param);

% Input OLTC into admittance matrix first time
[y1,y2,y3] = LV_tc1.sample(1,1,0,1,1);
Y(LV_tc1BusFrom,LV_tc1BusFrom) = Yorg(LV_tc1BusFrom,LV_tc1BusFrom) + y1+y2;
Y(LV_tc1BusTo,LV_tc1BusTo) = Yorg(LV_tc1BusTo,LV_tc1BusTo) +y1+y3;
Y(LV_tc1BusFrom,LV_tc1BusTo) = -y1;
Y(LV_tc1BusTo,LV_tc1BusFrom) = -y1;

% Solar PV power plant
% Medium voltage (bus 5)
% Solar irradiance
param.lat = 56.889;     % Latitude for Sørup (degrees)
param.t = 0.75;         % Transmittance (unitless)
param.S = 1367;         % Solar constant (w/m^2)
param.p = 100;          % Air pressure (Kpa)
param.Ts = Ts;

MV_si = solarIrradiance(param);

% PV
MV_pv1Bus = 5;
param.sBase = MV_sBase;
param.vBase = MV_vBase;
param.pRated = 7e6;
param.sMax = 7e6;
param.eta = 0.25;
param.A = 30000;

MV_pv1 = pvAsset(param);

% Wind power plant
% Medium Voltage (bus 12)
% Wind speed
param.Ts = Ts;
param.z = 100;

ws = windSpeed(param);
ws2 = windSpeed(param);

% WT
MV_wt1Bus = 12;
param.sBase = MV_sBase;
param.vBase = MV_vBase;
param.pRated = 7.5e6;
param.sMax = 7.5e6;
param.wMin = 3;
param.wMax = 25;
param.wRated = 12;
param.Ts = Ts;
param.onPU = false;

MV_wt1 = wtAsset(param);
MV_wt2 = wtAsset(param);

% Supermarkets
% Medium voltage (bus 9)
sm1Bus = 9;
param.sBase = MV_sBase;
param.Ts = Ts;
param.onPU = onPU;
param.pDefrost = [12e3 3.9e3];

sm1 = smAsset(param);

%% Setup Newton-Raphson load flow solver
maxIte = 100;
if onPU == true
    tol = 1e-7;
else
    tol = 1e-5;
end