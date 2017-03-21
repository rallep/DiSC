% 
clc; clear; close all;
Ts = 30;
N = 24*60*60/Ts;
rng(1);
onPU = false;
wMean = 10;
day = 1;
%% Setup the electrical grid
% Grid impedance data (inf indicates that a transformer will be placed between the two busses, this is done later in this script)
% format:   [from   to  Resistance(R[ohm])  Reactance(X[ohm])   length(l[km])]
Z =         [1      2   inf                  inf                  1;
             2      3   0.1                  0.1                 1;
             2      6   0.13                 0.09                5;
             3      4   0.13                 0.09                10;
             3      6   0.32                 0.15                10;
             4      5   0.13                 0.09                 7;
             5      8   0.32                 0.15                 15;
             6      7   0.1                  0.1                   5;
             7      8   0.05                 0.5                 10];
         
 %% Setup power flow module
% Number of busses
numBus = max(max(Z(:,1)),max(Z(:,2)));
% Set type of each bus (bus 1 is slack bus)
param.type = [0 1 1 1 1 1 1 1];  % 0 = slack bus, 1 = PQ bus and 2 = PV bus
% Set base voltage of each bus
param.vBase = [60e3 20e3*ones(1,numBus-1)]; 
% Create power flow module
pFlow = powerFlow(param);
% Get bus Admittance matrix and take a copy of the original. This is used
% for the tap-changing transformers
Y = pFlow.lDataToY(Z);
Yorg = Y;
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

% Power at busses
PF = 0.9;   % Power Factor
P1_sch = zeros(N,1);
Q1_sch = zeros(N,1);
P2_sch = -10*ones(N,1);
Q2_sch = P2_sch*tan(acos(PF));
P3_sch = -10*ones(N,1);
Q3_sch = P3_sch*tan(acos(PF));
P4_sch = -10*ones(N,1);
Q4_sch = P4_sch*tan(acos(PF));
P5_sch = -10*ones(N,1);
Q5_sch = P5_sch*tan(acos(PF));
P6_sch = -10*ones(N,1);
Q6_sch = P6_sch*tan(acos(PF));
P7_sch = -10*ones(N,1);
Q7_sch = P7_sch*tan(acos(PF));
P8_sch = -1e6*ones(N,1);
Q8_sch = P8_sch*tan(acos(PF));

Pin = [P1_sch P2_sch P3_sch P4_sch P5_sch P6_sch P7_sch P8_sch];
Qin = [Q1_sch Q2_sch Q3_sch Q4_sch Q5_sch Q6_sch Q7_sch Q8_sch];

%% Setup assets
% On the MV grid
% Wind power plant (Wind farm consisting of 10 2MW wind turbines)
MV_wpp_Bus = 6;         % Placement of wind power plant
param.numWt = 4;       % Number of wind turbines in the wind power plant
param.onPU = onPU;      % Indicate if simulation is on per unit
param.vBase = 20e3;     % Base voltage of the connection point
param.pRated = 1e6;     % Rated power of each wind turbine
param.sMax = 1e6;       % Maximum apparent power of each wind turbine
param.wMin = 3;         % Cut-in wind speed
param.wMax = 25;        % Cut-out wind speed
param.wRated = 12;      % Rated wind speed
param.Ts = Ts;          % Simulation time
param.z = 35;           % Height of each wind turbine from ground
% Create wind power plant object
MV_wpp = wppAsset(param);

%% Input to assets
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
MV_wpp_qRef = 0;         % Reactive power reference
MV_wpp.setQmode(0);  % Set reactive control mode (0 = constant power factor, 1 = follow reference, 2 = voltage droop control)

% Allocate memory
vOut = ones(N+1,numBus);    % Voltage output from power flow solution    
pSlack = zeros(N,1);        % Active power at slack bus
qSlack = zeros(N,1);        % Reactive power at slack bus
pWpp = zeros(N,1);          % Available wind power plant active power
MV_wpp_dP = zeros(N,1);     % Reference for the wind power plant
vMeas = zeros(N,1);
MV_wpp_ref = zeros(N,1);
pWppOut = zeros(N,1);
qWppOut = zeros(N,1);

for i=1:N
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set invironmental data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wIn = wMean+3*sin(i*0.005/(60/Ts));     % Mean wind changes during the day
    
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
        
        % Medium Voltage
        % Wind power plant
        [pWppOut(i),qWppOut(i),pWpp(i)] = MV_wpp.sample(wIn,abs(vOut(i-1,MV_wpp_Bus)),MV_wpp_dP(i),MV_wpp_dPlim,MV_wpp_qRef,MV_wpp_vRef,i);
        Pin(i,MV_wpp_Bus) = Pin(i,MV_wpp_Bus) + pWppOut(i);
        Qin(i,MV_wpp_Bus) = Qin(i,MV_wpp_Bus) + qWppOut(i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Electrical Grid  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [vOut(i,:),pSlack(i),qSlack(i),~] = pFlow.nrLoadFlow(Y,Pin(i,:)',Qin(i,:)');
end

%% Plotting
close all;
t = 0:N-1;
figure
plot(t,pSlack,t,qSlack)

figure
plot(t,abs(vOut(1:N,2:end))/20e3)

figure
plot(t,pWppOut/1e6,t,qWppOut/1e6)
legend('P','Q')
ylabel('Power [MW]')
xlabel('Time [min]')

figure
plot(t,abs(vOut(1:N,MV_wpp_Bus))/20e3)
