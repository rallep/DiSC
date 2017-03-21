% Example to illustrate how DiSC works
clc; close all; clear;
% This example illustrates some of the features of DiSC. All information is
% collected in this script to illustrate which steps must be taken to carry
% out a simulation. Some of this can be automated and will be in future
% releases.
%
% The progress bar seen when running the simulation degrades simulation
% speed, but can be removed for speed.
%% General setup
Ts = 60;                % [s]. Sampling time
N = (24*60*60/Ts)-1;    % [-]. Number of itterations (only one day of data) 
day = 60;               % [day]. Julian day of the year (1-365)
onPU = false;           % To indicate if the system should be simulated on a per unit basis (only works with on voltage level at the moment)
rng(1);                 % Seed random generator to replicate results

%% Input to simulation
% Here you can change the values to see what happens
% Mean wind speed
wMean = 9;
% Mean cloud cover
ccMean = 0.5;
% Grid Controller for keeping voltage at bus 5 at a reference value
Kp = -3000;             % [-]. Proportional gain (-3000 works)
Ki = -100;              % [-]. Integral gain (-100 works)
Ctrl_vRef = 20e3;       % [V]. Grid controller voltage reference
lossPr = 0.1;           % [-]. Loss probability
maxDelay = 10*60;       % [s]. Maximum allowable delay
delayDist = 'normal';   % ['uniform','normal']. Delay distribution
mu = 3*60;              % [s]. Mean of normal distribution
sigma = 2*60;           % [s]. Standard deviation of normal distribution
% Wind power plant input
wppQmode = 0;           % [0,1,2]. Wind power plant (WPP) reactive power control mode,
                        % 0 = constant power factor, 1 = follow reactive power
                        % reference and 2 = voltage droop control
                        % (NOTE: The system will become unstable with the default droop control. This is because of the size of the WPP compared to gain of droop controller).
wpp_qRef = 0;           % [VAR]. WPP reactive power reference


%% Setup the electrical grid
% Grid impedance data (inf indicates that a transformer will be placed between the two busses, this is done later in this script)
% format:   [from   to  Resistance(R[ohm])  Reactance(X[ohm])   length(l[km])]
Z =         [1      2   inf                 inf                 1;
             2      3   0.1                 0.1                 1;
             2      6   0.13                0.09                5;
             3      4   0.13                0.09                10;
             3      6   0.32                0.15                10;
             4      5   0.13                0.09                7;
             5      8   0.32                0.15                15;
             6      7   0.1                 0.1                 5;
             7      8   0.05                0.5                 10;
             5      9   inf                 inf                 1; 
             9      10  0.208               0.052               0.5;
             10     11  0.208               0.052               0.5;
             11     12  0.208               0.052               0.5;
             12     13  0.208               0.052               0.5;
             13     14  0.208               0.052               0.5;
             14     15  0.208               0.052               0.5];
 
%% Setup power flow module
% Number of busses
numBus = max(max(Z(:,1)),max(Z(:,2)));
% Set type of each bus (bus 1 is slack bus)
param.type = [0 ones(1,numBus-1)];  % 0 = slack bus, 1 = PQ bus and 2 = PV bus
% Set base voltage of each bus
param.vBase = [60e3 20e3*ones(1,7) 400*ones(1,7)]; % Bus 1 is HV grid, bus 2-8 is MV grid and bus 9-15 is LV grid 
% Create power flow module
pFlow = powerFlow(param);
% Get bus Admittance matrix and take a copy of the original. This is used
% for the tap-changing transformers
Y = pFlow.lDataToY(Z);
Yorg = Y; % Yorg is used in the tap-changing transformers, which alters the admittance matrix

%% Setup transformers
% From HV to MV grid (bus 1 to bus 2)
param.nomTapRatio = 3;      % Transformer nominal tap ratio
param.Z = 3+1i*13;           % Impedance of transformer
param.Ts = Ts;              % Sampling time
param.onPU = onPU;          % Indicate if simulation is on per unit
param.vBase = 20e3;         % Base voltage on secondary side of the transformer
MV_trafoBusFrom = 1;        % Bus where the transformer is going from
MV_trafoBusTo = 2;          % Bus where the transformer is going to
MV_trafo = tcAsset(param);  % Create transformer object
% Input transformer into admittance matrix for the first time
[y1,y2,y3] = MV_trafo.sample(1,1,0,1,1);
Y(MV_trafoBusFrom,MV_trafoBusFrom) = Yorg(MV_trafoBusFrom,MV_trafoBusFrom) + y1+y2;
Y(MV_trafoBusTo,MV_trafoBusTo) = Yorg(MV_trafoBusTo,MV_trafoBusTo) +y1+y3;
Y(MV_trafoBusFrom,MV_trafoBusTo) = -y1;
Y(MV_trafoBusTo,MV_trafoBusFrom) = -y1;

% From MV to LV (bus 5 to 9)
param.nomTapRatio = 50;     % Transformer nominal tap ratio
param.Z = 0.04+1i*0.04;     % Impedance of transformer
param.Ts = Ts;              % Sampling time
param.onPU = onPU;          % Indicate if simulation is on per unit
param.vBase = 400;          % Base voltage on secondary side of the transformer
LV_trafoBusFrom = 5;        % Bus where the transformer is going from
LV_trafoBusTo = 9;          % Bus where the transformer is going to
LV_trafo = tcAsset(param);  % Create transformer object
% Input transformer into admittance matrix for the first time
[y1,y2,y3] = LV_trafo.sample(1,1,0,1,1);
Y(LV_trafoBusFrom,LV_trafoBusFrom) = Yorg(LV_trafoBusFrom,LV_trafoBusFrom) + y1+y2;
Y(LV_trafoBusTo,LV_trafoBusTo) = Yorg(LV_trafoBusTo,LV_trafoBusTo) +y1+y3;
Y(LV_trafoBusFrom,LV_trafoBusTo) = -y1;
Y(LV_trafoBusTo,LV_trafoBusFrom) = -y1;

%% Load consumption data and upsample to sampling time (Ts)
% Residential consumption profiles
HouseData = load('../data/houseData400_oneDay.mat');  % Mat file containing consumption data
% Interpolate 15 min. consumption data to match sampling
t = 0:size(HouseData.Data.HouseP,1)-1;
ti = 0:(size(HouseData.Data.HouseP,1)/(60*60*(1/Ts)*size(HouseData.Data.HouseP,1)/4)):size(HouseData.Data.HouseP,1)-1;
HouseData.Data.pTs = interp1(t,HouseData.Data.HouseP(:,:),ti);
HouseData.Data.pTs = HouseData.Data.pTs.*1000;      % Data is originally in kW

% Industry
induData = load('../data/induPower');
t = 0:length(induData.p)-1;
ti = 0:(length(induData.p)/(60*60*(1/Ts)*length(induData.p))):length(induData.p)-1;
induData.pTs = interp1(t,induData.p(:,:),ti)';

% Agriculture
agriData = load('../data/agriPower');
t = 0:length(agriData.p)-1;
ti = 0:(length(agriData.p)/(60*60*(1/Ts)*length(agriData.p))):length(agriData.p)-1;
agriData.pTs = interp1(t,agriData.p(:,:),ti)';

% Commercial
commData = load('../data/commPower');
t = 0:length(commData.p)-1;
ti = 0:(length(commData.p)/(60*60*(1/Ts)*length(commData.p))):length(commData.p)-1;
commData.pTs = interp1(t,commData.p(:,:),ti)';

%% Setup all busses
% Power factors
pfIndu = 0.9;       % Power factor industry
pfAgri = 0.9;       % Power factor agriculture
pfResi = 0.97;      % Power factor Residential
pfComm = 0.95;      % Power factor commercial
% Setup active power on busses (negative means consumption)
Pbus1 = zeros(N,1);                                 % High voltage grid (slack bus)
Pbus2 = zeros(N,1);                                 % No production or consumption
Pbus3 = -sum(HouseData.Data.pTs(1:N,1:200),2);      % Aggregated 200 residential loads
Pbus4 = -4*commData.pTs(1:N);                       % 5 times commercial load
Pbus5 = zeros(N,1);                                 % Residential (LV grid, bus 9-15)
Pbus6 = zeros(N,1);                                 % Wind power plant (added later)
Pbus7 = -2*2*induData.pTs(1:N);                       % 5 times industrial load
Pbus8 = -3*agriData.pTs(1:N);                       % 6 times agricultural load
Pbus9 = -sum(HouseData.Data.pTs(1:N,201:220),2);    % 20 residential loads
Pbus10 = -sum(HouseData.Data.pTs(1:N,221:240),2);   % 20 residential loads
Pbus11 = -sum(HouseData.Data.pTs(1:N,241:260),2);   % 20 residential loads
Pbus12 = -sum(HouseData.Data.pTs(1:N,261:280),2);   % 20 residential loads
Pbus13 = -sum(HouseData.Data.pTs(1:N,281:300),2);   % 20 residential loads
Pbus14 = -sum(HouseData.Data.pTs(1:N,301:320),2);   % 20 residential loads
Pbus15 = -sum(HouseData.Data.pTs(1:N,321:340),2);   % 20 residential loads
% Form Pin matrix used for power flow simulation (we are allocating memory to speed up simulation)
Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7 Pbus8 Pbus9 Pbus10 ...
        Pbus11 Pbus12 Pbus13 Pbus14 Pbus15];
% Setup reactive power on busses (based on power factor)
Qbus1 = zeros(N,1);                                 % High voltage grid (slack bus)
Qbus2 = zeros(N,1);                                 % No production or consumption
Qbus3 = Pbus3*tan(acos(pfResi));                    % Aggregated 200 residential loads
Qbus4 = Pbus4*tan(acos(pfComm));                    % 2 times commercial load
Qbus5 = zeros(N,1);                                 % Residential (LV grid, bus 9-15)
Qbus6 = zeros(N,1);                                 % Wind power plant (added later)
Qbus7 = Pbus7*tan(acos(pfIndu));                    % Industrial load
Qbus8 = Pbus8*tan(acos(pfAgri));                    % 3 times agricultural load
Qbus9 = Pbus9*tan(acos(pfResi));                    % 20 residential loads
Qbus10 = Pbus10*tan(acos(pfResi));                  % 20 residential loads
Qbus11 = Pbus11*tan(acos(pfResi));                  % 20 residential loads
Qbus12 = Pbus12*tan(acos(pfResi));                  % 20 residential loads
Qbus13 = Pbus13*tan(acos(pfResi));                  % 20 residential loads
Qbus14 = Pbus14*tan(acos(pfResi));                  % 20 residential loads
Qbus15 = Pbus15*tan(acos(pfResi));                  % 20 residential loads
% Form Qin matrix used for power flow simulation (we are allocating memory to speed up simulation)
Qin = [Qbus1 Qbus2 Qbus3 Qbus4 Qbus5 Qbus6 Qbus7 Qbus8 Qbus9 Qbus10 ...
        Qbus11 Qbus12 Qbus13 Qbus14 Qbus15];

%% Setup communication links
% Grid controller to wind power plant
param.Ts = Ts;
param.maxDelay = maxDelay;
param.numLinks = 1;
param.mu = mu;
param.sigma = sigma;
clinkCtrlToWpp = comLink(param);
clinkCtrlToWpp.setPrLoss(lossPr);
clinkCtrlToWpp.setInverseCDFdelay(delayDist,param);
    
%% Setup assets
% On the MV grid
% Wind power plant (Wind farm consisting of 10 2MW wind turbines)
MV_wpp_Bus = 5;         % Placement of wind power plant
param.numWt = 10;       % Number of wind turbines in the wind power plant
param.onPU = onPU;      % Indicate if simulation is on per unit
param.vBase = 20e3;     % Base voltage of the connection point
param.pRated = 2e6;     % Rated power of each wind turbine
param.sMax = 2e6;       % Maximum apparent power of each wind turbine
param.wMin = 3;         % Cut-in wind speed
param.wMax = 25;        % Cut-out wind speed
param.wRated = 12;      % Rated wind speed
param.Ts = Ts;          % Simulation time
param.z = 35;           % Height of each wind turbine from ground
% Create wind power plant object
MV_wpp = wppAsset(param);

% On the LV grid
% PV systems (4 on bus 10, 4 on bus 11, 4 on bus 14 and 6 on bus 15)
% Set parameters for PV system on bus 10 and 11
param.numPv = 4;        % Number of PVs
param.onPU = onPU;      % Indicate if simulation is on per unit
param.vBase = 400;      % Voltage base
param.pRated = 6e3;     % Rated power of inverter
param.sMax = 6e3;       % Maximum apparent power of inverter
param.eta = 0.25;       % Efficiency of solar cells
param.A = 25;           % Area of solar cells
param.Ts = Ts;          % Sampling time
param.lat = 57;         % Latitude of PV system (57 degrees is Aalborg, Denmark)            
param.t = 0.75;         % Transmittance at location        
param.p = 100;          % Air pressure at location
% Create PV systems on bus 10 and 11
LV_PV10_Bus = 10;       % Bus where the system is placed
LV_PV10 = sppAsset(param);
LV_PV11_Bus = 11;       % Bus where the system is placed
LV_PV11 = sppAsset(param);
% Set parameters for PV systems on bus 14
param.numPv = 4;        % Number of PVs
% Create PV systems on bus 14
LV_PV14_Bus = 14;       % Bus where the system is placed
LV_PV14 = sppAsset(param);
% Set parameters for PV systems on bus 15
param.numPv = 6;        % Number of PVs
% Create PV systems on bus 15
LV_PV15_Bus = 15;       % Bus where the system is placed
LV_PV15 = sppAsset(param);

%% Setup simulation
% Input to assets
% MV grid
% Transformer
MV_trafo_vRef = 20e3;       % Voltage reference for local voltage control on secondary side
MV_trafo_uRef = 0;          % Tap position reference
MV_trafo.setMode(1);        % Control mode of tap-changing transformer (0 = follow tap reference, 1 = voltage control)
MV_trafo.setTapSpec(10,10,50,5*60,0.0125,0.04); % Set transformer parameters (see tcAsset documentation)
% Wind power plant
MV_wpp_vRef = 20e3;         % Voltage reference
MV_wpp_dP = 0;              % Power reference
MV_wpp_dPlim = 0;           % Derate power
MV_wpp_qRef = wpp_qRef;     % Reactive power reference
MV_wpp.setQmode(wppQmode);  % Set reactive control mode (0 = constant power factor, 1 = follow reference, 2 = voltage droop control)

% LV grid
% Transformer
LV_trafo_vRef = 400;        % Voltage reference for local voltage control on secondary side
LV_trafo_uRef = 0;          % Tap position reference
LV_trafo.setMode(1);        % Control mode of tap-changing transformer (0 = follow tap reference, 1 = voltage control)
LV_trafo.setTapSpec(10,10,15,5*60,0.0125,0.03);  % Set transformer parameters (see tcAsset documentation)
% Solar PV systems input (assumes identical inputs to all systems)
LV_spp_vRef = 400;          % Voltage reference
LV_spp_dP = 0;              % Power reference
LV_spp_dPlim = 0;           % Derate power
LV_spp_qRef = 0;            % Reactive power reference

% Grid Control
errI = zeros(N,1);
v=20e3;
pRef = 0;
qTest = 0;

% Allocate memory
vOut = ones(N+1,numBus);    % Voltage output from power flow solution    
pSlack = zeros(N,1);        % Active power at slack bus
qSlack = zeros(N,1);        % Reactive power at slack bus
pWpp = zeros(N,1);          % Available wind power plant active power
MV_wpp_dP = zeros(N,1);     % Reference for the wind power plant
vMeas = zeros(N,1);         % Measured voltage (used for the grid controller)
MV_wpp_ref = zeros(N,1);    % Reference for the WPP active power
pWppOut = zeros(N,1);       % Active power output of the WPP
qWppOut = zeros(N,1);       % Reactive power output of the WPP

%% Run simulation
h = waitbar(0,'Simulation is 0% finished');
tic % Time
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
    wIn = wMean+3*sin(i*0.005);     % Mean wind changes during the day
    ccIn = ccMean+0.5*sin(i*0.009); % Mean cloud cover changes during the day 
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert assets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
    else
        % Tap changing transformers
        % MV grid
        [y1,y2,y3] = MV_trafo.sample(abs(vOut(i-1,MV_trafoBusTo)),MV_trafo_vRef,MV_trafo_uRef,i,day);
        Y(MV_trafoBusFrom,MV_trafoBusFrom) = Yorg(MV_trafoBusFrom,MV_trafoBusFrom) + y1+y2;
        Y(MV_trafoBusTo,MV_trafoBusTo) = Yorg(MV_trafoBusTo,MV_trafoBusTo) +y1+y3;
        Y(MV_trafoBusFrom,MV_trafoBusTo) = -y1;
        Y(MV_trafoBusTo,MV_trafoBusFrom) = -y1;
        
        % LV grid
        [y1,y2,y3] = LV_trafo.sample(abs(vOut(i-1,LV_trafoBusTo)),LV_trafo_vRef,LV_trafo_uRef,i,day);
        Y(LV_trafoBusFrom,LV_trafoBusFrom) = Yorg(LV_trafoBusFrom,LV_trafoBusFrom) + y1+y2;
        Y(LV_trafoBusTo,LV_trafoBusTo) = Yorg(LV_trafoBusTo,LV_trafoBusTo) +y1+y3;
        Y(LV_trafoBusFrom,LV_trafoBusTo) = -y1;
        Y(LV_trafoBusTo,LV_trafoBusFrom) = -y1;
        
        % Medium Voltage
        % Wind power plant
        [pWppOut(i),qWppOut(i),pWpp(i)] = MV_wpp.sample(wIn,abs(vOut(i-1,MV_wpp_Bus)),pRef,MV_wpp_dPlim,MV_wpp_qRef,MV_wpp_vRef,i);
        Pin(i,MV_wpp_Bus) = Pin(i,MV_wpp_Bus) + pWppOut(i);
        Qin(i,MV_wpp_Bus) = Qin(i,MV_wpp_Bus) + qWppOut(i);
        
        % Low voltage 
        % PV systems
        % On bus 10
        [p,q] = LV_PV10.sample(ccIn,abs(vOut(i-1,LV_PV10_Bus)),LV_spp_dP,LV_spp_dPlim,LV_spp_qRef,LV_spp_vRef,i,day);
        Pin(i,LV_PV10_Bus) = Pin(i,LV_PV10_Bus) + p;
        Qin(i,LV_PV10_Bus) = Qin(i,LV_PV10_Bus) + q;
        % On bus 11
        [p,q] = LV_PV11.sample(ccIn,abs(vOut(i-1,LV_PV11_Bus)),LV_spp_dP,LV_spp_dPlim,LV_spp_qRef,LV_spp_vRef,i,day);
        Pin(i,LV_PV11_Bus) = Pin(i,LV_PV11_Bus) + p;
        Qin(i,LV_PV11_Bus) = Qin(i,LV_PV11_Bus) + q;
        % On bus 14
        [p,q] = LV_PV14.sample(ccIn,abs(vOut(i-1,LV_PV14_Bus)),LV_spp_dP,LV_spp_dPlim,LV_spp_qRef,LV_spp_vRef,i,day);
        Pin(i,LV_PV14_Bus) = Pin(i,LV_PV14_Bus) + p;
        Qin(i,LV_PV14_Bus) = Qin(i,LV_PV14_Bus) + q;
        % On bus 15
        [p,q] = LV_PV15.sample(ccIn,abs(vOut(i-1,LV_PV15_Bus)),LV_spp_dP,LV_spp_dPlim,LV_spp_qRef,LV_spp_vRef,i,day);
        Pin(i,LV_PV15_Bus) = Pin(i,LV_PV15_Bus) + p;
        Qin(i,LV_PV15_Bus) = Qin(i,LV_PV15_Bus) + q;    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Electrical Grid  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [vOut(i,:),pSlack(i),qSlack(i),~] = pFlow.nrLoadFlow(Y,Pin(i,:)',Qin(i,:)');
    if pFlow.nIte >= 100
        error('N-R solver not converging')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Communication from assets to controller
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vMeas(i) = clinkCtrlToWpp.sampleIn(i,vOut(i,5));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Grid Control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if incomming measurement is NaN (Indicates if data is lost).
    if isnan(vMeas(i))
        v = v;
    else
        v = vMeas(i);
    end
    % Apply PI control
    err = (abs(v)-Ctrl_vRef);               % Error
    errI(i+1) = errI(i) + err;              % Integral of error
    MV_wpp_ref(i) = -Kp*err + Ki*errI(i);   % Apply control
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Communication from controller to assets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    MV_wpp_dP(i) = clinkCtrlToWpp.sampleOut(i,MV_wpp_ref(i));
    % Check if outgoing is NaN (indicates if data is lost)
    if isnan(MV_wpp_dP(i))
        pRef = pRef;
    else
        pRef = MV_wpp_dP(i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Progress bar
    %%%%~%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~mod(i,20)
        waitbar(i/N,h,sprintf('Simulation is %.0f%% finished...',100*i/N))
    end
end
toc
close(h);

%% Display data
disp('HV/MV transformer tap position:')
MV_trafo.tapPos
disp('MV/LV transformer tap position:')
LV_trafo.tapPos


%% Plotting
close all;
t = (0:N-1)/60*60/Ts;

% P and Q slack
figure
plot(t,pSlack/1e6,t,qSlack/1e6)
xlabel('Time [hrs]')
ylabel('Pslack / Qslack [MW / MVar]')
legend('Active Power','Reactive Power')
title('P/Q slack')

% Power output of wind power plant
figure
plot(t,pWpp/1e6,t,pWppOut/1e6,t,qWppOut/1e6)
xlabel('Time [hrs]')
ylabel('Power [MW]')
legend('Available Power','P Output','Q Output')
title('Wind Power Plant')

% Voltages
figure
subplot(2,1,1)
plot(t,abs(vOut(1:N,2:8))/20e3,t,ones(1,N)*0.9,'--r',t,ones(1,N)*1.1,'--r')
title('Voltages MV')
xlabel('Time [hrs]')
ylabel('Volatges [PU]')
ylim([0.8 1.2])
subplot(2,1,2)
plot(t,abs(vOut(1:N,9:end))/400,t,ones(1,N)*0.9,'--r',t,ones(1,N)*1.1,'--r')
title('Voltages LV')
xlabel('Time [hrs]')
ylabel('Volatges [PU]')
ylim([0.8 1.2])

