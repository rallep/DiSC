% Script for setting up the LV benchmark grid
%
%
% R. Pedersen 5-26-2014, Aalborg University

% Base parameters for converting to p.u. 
Sbase = 400e3;          % Complex power base [VA].
Vbase = 0.4e3;          % Voltage base [V].
Zbase = Vbase^2/Sbase;  % Impedance base [Ohm].

% Line impedance matrix, with trafo included
Z = benchmarkImpedanceMatrix('LV');

% Normalize to p.u. 
Z(:,3:4) = Z(:,3:4)./Zbase;

% Form Bus admittance matrix
Y = ztoybus(Z);

%% Load consumption data
HouseData = load('consData/house1to120days2.mat');  % Mat file containing consumption data

% Interpolate 15 min. consumption data to match sampling
t = 0:length(HouseData.Data.HouseP)-1;
ti = 0:(length(HouseData.Data.HouseP)/(60*60*(1/Ts)*length(HouseData.Data.HouseP)/4)):length(HouseData.Data.HouseP)-1;
HouseData.Data.pTs = interp1(t,HouseData.Data.HouseP(:,:),ti);

numSamples = length(HouseData.Data.pTs);         % 

%% Setup Busses
numBus = length(Y);
% Type
type = [0 ones(1,numBus-1)]; % 0=slack, 1=PQ 
% Power factor
PF = 0.97;
% Voltages in p.u.
Vin = [1 zeros(1,numBus-1)]; % Bus1 = slack
% P at each bus
Pbus1 = zeros(numSamples,1);
Pbus2 = zeros(numSamples,1);
Pbus3 = -sum(HouseData.Data.pTs(:,1:8),2).*1000;
Pbus4 = -sum(HouseData.Data.pTs(:,9:12),2).*1000;
Pbus5 = -sum(HouseData.Data.pTs(:,13:15),2).*1000;
Pbus6 = -sum(HouseData.Data.pTs(:,16:19),2).*1000;
Pbus7 = -sum(HouseData.Data.pTs(:,20:22),2).*1000;
Pbus8 = -sum(HouseData.Data.pTs(:,23:27),2).*1000;
Pbus9 = zeros(numSamples,1);
Pbus10 = -sum(HouseData.Data.pTs(:,28:31),2).*1000;
Pbus11 = -sum(HouseData.Data.pTs(:,32:35),2).*1000;
Pbus12 = -sum(HouseData.Data.pTs(:,36:37),2).*1000;
Pbus13 = -sum(HouseData.Data.pTs(:,38:38),2).*1000;
Pbus14 = -sum(HouseData.Data.pTs(:,39:40),2).*1000;
Pbus15 = -sum(HouseData.Data.pTs(:,41:45),2).*1000;
Pbus16 = -sum(HouseData.Data.pTs(:,46:49),2).*1000;
Pbus17 = -sum(HouseData.Data.pTs(:,50:52),2).*1000;
Pbus18 = -sum(HouseData.Data.pTs(:,53:54),2).*1000;
Pbus19 = -sum(HouseData.Data.pTs(:,55:58),2).*1000;
Pbus20 = -sum(HouseData.Data.pTs(:,59:61),2).*1000;
Pbus21 = -sum(HouseData.Data.pTs(:,62:63),2).*1000;
Pbus22 = -sum(HouseData.Data.pTs(:,64:67),2).*1000;
Pbus23 = -sum(HouseData.Data.pTs(:,68:70),2).*1000;
Pbus24 = -sum(HouseData.Data.pTs(:,71:72),2).*1000;
Pbus25 = -sum(HouseData.Data.pTs(:,73:75),2).*1000;
Pbus26 = -sum(HouseData.Data.pTs(:,76:77),2).*1000;
Pbus27 = -sum(HouseData.Data.pTs(:,78:80),2).*1000;
Pbus28 = -sum(HouseData.Data.pTs(:,81:81),2).*1000;
Pbus29 = -sum(HouseData.Data.pTs(:,82:83),2).*1000;
Pbus30 = -sum(HouseData.Data.pTs(:,84:85),2).*1000;
Pbus31 = -sum(HouseData.Data.pTs(:,86:87),2).*1000;
Pbus32 = -sum(HouseData.Data.pTs(:,88:90),2).*1000;
Pbus33 = zeros(numSamples,1);
Pbus34 = -sum(HouseData.Data.pTs(:,91:94),2).*1000;
Pbus35 = -sum(HouseData.Data.pTs(:,95:97),2).*1000;
Pbus36 = -sum(HouseData.Data.pTs(:,98:104),2).*1000;
Pbus37 = -sum(HouseData.Data.pTs(:,105:106),2).*1000;
Pbus38 = -sum(HouseData.Data.pTs(:,107:108),2).*1000;
Pbus39 = -sum(HouseData.Data.pTs(:,109:109),2).*1000;
Pbus40 = -sum(HouseData.Data.pTs(:,110:111),2).*1000;
Pbus41 = -sum(HouseData.Data.pTs(:,112:113),2).*1000;
Pbus42 = -sum(HouseData.Data.pTs(:,114:116),2).*1000;

Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7 Pbus8 Pbus9 Pbus10 Pbus11 Pbus12 Pbus13 Pbus14...
    Pbus15 Pbus16 Pbus17 Pbus18 Pbus19 Pbus20 Pbus21 Pbus22 Pbus23 Pbus24 Pbus25 Pbus26 Pbus27...
    Pbus28 Pbus29 Pbus30 Pbus31 Pbus32 Pbus33 Pbus34 Pbus35 Pbus36 Pbus37 Pbus38 Pbus39 Pbus40...
    Pbus41 Pbus42]/Sbase;

% Q at each bus
Qin = Pin.*tan(acos(PF));

%% Setup PV systems
randPlace = true;
numPV = 40;
busInt = [3 42];
sMaxInt = [3e3 8e3];
pRatedIntMax = [3e3 8e3];
areaInt = [25 35];
etaInt = [0.18 0.25];

[pvSys, pvBus] = setupPVsystems(randPlace,numPV,Sbase,Vbase,busInt,sMaxInt,pRatedIntMax,areaInt,etaInt);

% Set reactive control mode and power factor
for i=1:numPV
    pvSys(i).setQmode(0);
    pvSys(i).setPF(1);
end

%% Setup solar irradiance objects
% Cloud cover
wEnvS = sin(Ts/30000*(0:1:numSamples-1)');
% Solar irradiance
numSI = numPV;
latitude = 56.889;      % Latitude for Sørup (degrees)
pAirInt = [90 110];     % Air pressure interval
transInt = [0.65 0.85]; % Transmittance interval

siSys = setupSolarIrradiance(Ts,numSI,latitude,pAirInt,transInt);

%% Setup energy storage systems
randPlace = false;
numES = 1;
busInt = [2];
sMaxInt = [3e3 3e3];
pRatedIntMax = [3e3 3e3];
pRatedIntMin = [-3e3 -3e3];
eRatedInt = [3e6 3e6]*60*60/Ts;        % [MWh]

[esSys, esBus] = setupESsystems(Ts,randPlace,numES,Sbase,Vbase,busInt,sMaxInt,pRatedIntMax,pRatedIntMin,eRatedInt);

for i=1:numES
    esSys(i).setPmode(0);
    esSys(i).setPF(0.95);
end

%% Setup Newton-Raphson load flow solver
maxIte = 100;
tol = 1e-6;

