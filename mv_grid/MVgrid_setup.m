% Setup and data for the MV benchmark grid
%% Setup Base for convertion to p.u.
Sbase = 15e3;               % [VA], 15 MVA
Vbase = 20e3;               % [V], 20 kV
Zbase = Vbase^2/Sbase;      % [ohm]

%% Setup the grid

% High voltage grid parameters
SCP = 249e6;                % [VA], Short Circuit Power 249
SCR = SCP/Sbase;            % Short Circuit Ratio
XtoR = 0.1;                 % Resistance to reactance relationship (R/X)

% High voltage grid impedance (Not sure if this is needed)
Rgrid = (1/SCR) * cos(atan(XtoR));
Xgrid = (1/SCR) * sin(atan(XtoR));
Rgrid_pu = Rgrid/Zbase;
Xgrid_pu = Xgrid/Zbase;

% Line impedance matrix
Z = benchmarkImpedanceMatrix('MV');
Z(:,3:4) = Z(:,3:4)./Zbase;

% Calculate bus admitance matrix
Y = ztoybus(Z);

%% Setup busses
type = [0 ones(1,length(Y)-1)]; % 0=slack, 1=PQ
% Power factors
pfIndu = 0.9;       % Power factor industry
pfAgri = 0.9;       % Power factor agriculture
pfResi = 0.97;      % Power factor Residential
pfComm = 0.95;      % Power factor commercial

% P and Q values for 24 hours sampled every hour (Summer period)
% Industry
PinduSummer = [41 42 45.72 76.6 58.2 97.96 107.44 117.84 111.64 115 111.84 115.76 117.68 121.2 112.2 105.04 91.84 ...
                60.2 58 39.88 38.24 38.12 37.88 41.08]*1000;
QinduSummer = PinduSummer*tan(acos(pfIndu));
% Agriculture
PagriSummer = [4.18 3.55 3.59 3.66 4.02 4.06 4.37 4.09 5.22 6.05 5.67 5.03 3.47 5.47 3.41 3.6 5.18 7.1 6.7 5.61 ...
                5.28 5.41 6.09 4.05]*1000;
QagriSummer = PagriSummer*tan(acos(pfAgri));
% Residential
PresiSummer = [61.29 57.14 55.07 54.85 63.21 80.43 101.23 104.58 104.49 86.58 85.61 77.35 81.45 83.26 85.52 98.8 ...
                142.74 187.92 168.89 136.94 117.25 114.58 100.53 86.58]*1000;
QresiSummer = PresiSummer*tan(acos(pfResi)); 
% Commercial
PcommSummer = [13.62 12.88 12.73 12.95 12.59 12.88 13.77 14.74 21.92 20.51 20.13 18.18 24.17 24.34 25.76 26.27 31.17 ...
                22.31 16.22 17.12 15.52 15.88 15.16 15.44]*1000;
QcommSummer = PcommSummer*tan(acos(pfComm));

% Wind power plant
PWP = 8e3+3e3.*randn(1,length(PinduSummer))*0;
QWP = PWP*tan(acos(0.9));

% PV power plant
PPV = 8e3+4e3.*randn(1,length(PinduSummer))*0;
QPV = PPV*tan(acos(0.9));

% P and Q on each bus in p.u.
Pbus1 = zeros(1,length(PinduSummer));
Qbus1 = zeros(1,length(PinduSummer));
Pbus2 = zeros(1,length(PinduSummer));
Qbus2 = zeros(1,length(PinduSummer));
Pbus3 = zeros(1,length(PinduSummer));
Qbus3 = zeros(1,length(PinduSummer));
Pbus4 = (-PinduSummer./Sbase);
Qbus4 = (-QinduSummer./Sbase);
Pbus5 = PPV./Sbase;
Qbus5 = QPV./Sbase;
Pbus6 = (-PcommSummer./Sbase);
Qbus6 = (-QcommSummer./Sbase);
Pbus7 = zeros(1,length(PinduSummer));
Qbus7 = zeros(1,length(PinduSummer));
Pbus8 = -PagriSummer./Sbase;
Qbus8 = -QagriSummer./Sbase;
Pbus9 = (-PcommSummer./Sbase);
Qbus9 = (-QcommSummer./Sbase);
Pbus10 = (-PresiSummer./Sbase);
Qbus10 = (-QresiSummer./Sbase);
Pbus11 = zeros(1,length(PinduSummer));
Qbus11 = zeros(1,length(PinduSummer));
Pbus12 = PWP./Sbase;
Qbus12 = QWP./Sbase;

Pin = [Pbus1' Pbus2' Pbus3' Pbus4' Pbus5' Pbus6' Pbus7' Pbus8' Pbus9' Pbus10' Pbus11' Pbus12'];
Qin = [Qbus1' Qbus2' Qbus3' Qbus4' Qbus5' Qbus6' Qbus7' Qbus8' Qbus9' Qbus10' Qbus11' Qbus12'];
% Voltage at each bus in p.u.
Vbus1 = 1.0;
Vbus2 = 0;
Vbus3 = 0;
Vbus4 = 0;
Vbus5 = 0;
Vbus6 = 0;
Vbus7 = 0;
Vbus8 = 0;
Vbus9 = 0;
Vbus10 = 0;
Vbus11 = 0;
Vbus12 = 0;

Vin = [Vbus1 Vbus2 Vbus3 Vbus4 Vbus5 Vbus6 Vbus7 Vbus8 Vbus9 Vbus10 Vbus11 Vbus12];

% Setup tolerence and max iterations
tol = 0.0005;
maxIte = 100;