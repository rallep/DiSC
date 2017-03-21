% Test of tap changing transformer
clc; close all; clear;

% Simulation parameters
Ts = 60;
day = 1;
N=24*60*60/Ts;
onPU = false;

% Setup simple bus systen
Sbase = 50e6;
Vbase = 20e3;
Zbase = Vbase^2/Sbase;


% Format:
%     [from       to      R                           X                             l]
Z =   [1          2       0.3                         0.13                          1;
       2          3       0.13                        0.1                           1; 
       3          4       0.4                         0.4                           1;
       4          5       0.13                        0.1                           1;
       4          6       0.13                        0.1                           1;
       5          6       0.13                        0.1                           1;
       6          7       0.13                        0.1                           1;
       ];

% Setup power flow object
param.numBus = max(max(Z(:,1)),max(Z(:,2)));
param.type = [0 ones(1,param.numBus-1)];
param.vBase = 400*ones(param.numBus,1);

pFlow = powerFlow(param);

Y =pFlow.lDataToY(Z);
Yorg = Y;
%%

% sw Objects
param.zBase = Zbase;
param.Ts = Ts;
param.Z = 0.02-1i*0.03;
param.onPU = onPU;
sw1 = swAsset(param);
uRef = 0;

% Power factor
PF = 0.97;



Pbus1 = zeros(N,1);
Qbus1 = zeros(N,1);
Pbus2 = -5e3*ones(N,1);
Qbus2 = Pbus2*tan(acos(PF));
Pbus3 = zeros(N,1);
Qbus3 = zeros(N,1);
Pbus4 = zeros(N,1);
Qbus4 = zeros(N,1);
Pbus5 = -5e3*ones(N,1);
Qbus5 = Pbus5*tan(acos(PF));
Pbus6 = -5e3*ones(N,1);
Qbus6 = Pbus6*tan(acos(PF));
Pbus7 = -5e3*ones(N,1);
Qbus7 = Pbus7*tan(acos(PF));

Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7];
Qin = [Qbus1 Qbus2 Qbus3 Qbus4 Qbus5 Qbus6 Qbus7];


%% Setup Newton-Raphson load flow solver
Vout = ones(N,param.numBus);
Pslack = zeros(N,1);
Qslack = zeros(N,1);
for i=1:N
    if ~mod(i,24*60*60/Ts)
       if day == 365
           day = 1;
       else
           day = day +1;
       end
    end
    
    if i == 500
        uRef = 1;
    end
    if i == 700
        uRef = 0;
    end
    % Check switch state
    [state,y] = sw1.sample(uRef,i,day);
    if state
        Y(2,2) = Yorg(2,2)+y;
        Y(7,7) = Yorg(7,7)+y;
        Y(2,7) = -y;
        Y(7,2) = -y;
    else
        Y = Yorg;
    end
    
    [Vout(i,:),Pslack(i),Qslack(i),~]=pFlow.nrLoadFlow(Y,Pin(1,:)',Qin(i,:)');
    if pFlow.nIte == 100
        error('N-R not converging');
    end
end
%% Plotting
t = (0:N-1)/(60*60/Ts);
figure
plot(t,abs(Vout)/400)
title('All Voltages')
ylabel('Voltages [p.u.]')
xlabel('Time [hrs]')

figure
plot(t,abs(Vout(:,7))/400)
title('Voltage at Bus 7')
ylabel('Voltage [p.u.]')
xlabel('Time [hrs]')