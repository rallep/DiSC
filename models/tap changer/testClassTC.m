% Test of tap changing transformer
clc; clear; close all;
Ts = 60;
day = 1;
N=(24/24)*24*60*60/Ts;

% Setup simple bus system
VbaseMV = 20e3;

Z = [1      2       inf         inf         1;
     1      3       inf         inf         1;
     2      4       0.13        0.1         10;
     3      5       0.13        0.1         10];

% Power flow object
% Setup busses
numBus = max(Z(:,2));
% Type
type = [0 ones(1,numBus-1)]; % 0=slack, 1=PQ 
% Power factor
PF = 0.97;
% Voltages in p.u.
param.type = type;
param.vBase = [60e3 20e3*ones(1,numBus-1)]; 
nrpf = powerFlow(param);

Y = nrpf.lDataToY(Z);
Yorg = Y;
%%

% tc Objects
param.Z = 3+1i*13;
param.Ts = Ts;
param.onPU = false;
param.nomTapRatio = 3;
param.vBase = 20e3;

MV_tc1 = tcAsset(param);
MV_tc1.setMode(1);
MV_tc1_uRef = 0;
MV_tc1_vRef = 20e3;

MV_tc2 = tcAsset(param);
MV_tc2.setMode(1);
MV_tc2_vRef = 20e3;
MV_tc2_uRef = 0;


Pbus1 = zeros(N,1);
Qbus1 = zeros(N,1);
Pbus2 = zeros(N,1);
Qbus2 = zeros(N,1);
Pbus3 = zeros(N,1);
Qbus3 = zeros(N,1);
Pbus4 = -40e5*ones(N,1);
Qbus4 = Pbus4*tan(acos(PF));
Pbus5 = -70e5*ones(N,1);
Qbus5 = Pbus5*tan(acos(PF));

Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5];
Qin = [Qbus1 Qbus2 Qbus3 Qbus4 Qbus5];


%% Setup Newton-Raphson load flow solver
maxIte = 100;
tol = 1e-6;
Vout = ones(N,numBus);
Pslack = zeros(N,1);
Qslack = zeros(N,1);
MV_tc1Tap = zeros(N,1);
MV_tc2Tap = zeros(N,1);

[y1,y2,y3] = MV_tc1.sample(VbaseMV,VbaseMV,0,1,1);
Y(1,1) = Yorg(1,1) + y1+y2;
Y(2,2) = Yorg(2,2) +y1+y3;
Y(1,2) = Yorg(1,2)-y1;
Y(2,1) = Yorg(2,1)-y1;

[y1,y2,y3] = MV_tc2.sample(VbaseMV,VbaseMV,0,1,1);
Y(1,1) = Y(1,1) + y1+y2;
Y(3,3) = Yorg(3,3) + y1+y3;
Y(1,3) = Yorg(1,3)-y1;
Y(3,1) = Yorg(3,1)-y1;

for i=1:N
    if ~mod(i,24*60*60/Ts)
       if day == 365
           day = 1;
       else
           day = day +1;
       end
    end
    
    %%%% Insert disturbance
    if i>500
        Pin(i,4) = 10e5;
        %Pin(i,4) = 60e3;
    end
    
    if i>1000
        Pin(i,4) = -60e5;
    end

    %%%& Insert tap changers
    if i == 1
        
    else
        [y1,y2,y3] = MV_tc1.sample(abs(Vout(i-1,2)),MV_tc1_vRef,MV_tc1_uRef,i,day);
        Y(1,1) = Yorg(1,1) + y1+y2;
        Y(2,2) = Yorg(2,2) +y1+y3;
        Y(1,2) = -y1;
        Y(2,1) = -y1;
        MV_tc1Tap(i) = MV_tc1.tapPos;
        
        [y1,y2,y3] = MV_tc2.sample(abs(Vout(i-1,3)),MV_tc2_vRef,MV_tc2_uRef,i,day);
        Y(1,1) = Y(1,1) + y1+y2;
        Y(3,3) = Yorg(3,3) +y1+y3;
        Y(1,3) = -y1;
        Y(3,1) = -y1;
        MV_tc2Tap(i) = MV_tc2.tapPos;
    end
    
    [Vout(i,:),Pslack(i),Qslack(i),~]=nrpf.nrLoadFlow(Y,Pin(i,:)',Qin(i,:)');
end
Vout;
Pslack;
Qslack;
%% Plotting
t = 0:N-1;

figure
subplot(2,1,1)
plot(t,abs(Vout(:,2:5))/VbaseMV)
ylabel('Voltage MV [PU]')
legend('B2','B3','B4','B5')
subplot(2,1,2)
plot(t,MV_tc1Tap,t,MV_tc2Tap)
ylabel('Tap Position [-]')
xlabel('Time [min]')
legend('TC MV1','TC MV2')
