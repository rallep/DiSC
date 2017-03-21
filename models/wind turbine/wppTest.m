clc
clear
close all

%% Inputs
pRated = 3e5;   % [W]
sMax = 3e6;     % [VA]
wMin = 3;       % [m/s]
wMax = 25;      % [m/s]
wRated = 12;    % [m/s]
tauEnv = 1;     % [s]
tauWT = 1;      % [s]
Ts = 1;         % [min]
% Gain from wEnv^3 to p
kWT = pRated/(wRated^3);

sBase = 3e6;    % [VA]
vBase = 20e3;   % [V]

par.sBase = sBase;
par.vBase = vBase;
par.pRated = pRated;
par.sMax = sMax;
par.wMin = wMin;
par.wMax = wMax;
par.wRated = wRated;
par.kWT = kWT;
N = 50;
v = 0.1*sin(0.3*(0:1:N-1)'+pi)+1;
wEnv = 12*sin(0.3*(0:1:N-1)')*0+24;
dP = -0.1;
dPlim = 0.2;
Qmode = 1;
PF = 0.95;
qRef = 0.2;
vRef = 1;

WT = wtAsset(par);
WT.setPF(PF);
WT.setQmode(Qmode);


%% Test wtAsset
for i=1:N
   [p(i),q(i)] = wtAssetFun(sBase,vBase,v(i),wEnv(i),par,dP,dPlim,Qmode,PF,qRef,vRef); 
   [pp(i),qq(i)] = WT.sample(wEnv(i),v(i),dP,dPlim,qRef,vRef);
end

%% Plotting
t = 0:N-1;
figure
plot(t,p,t,pp+1,t,q,t,qq+1)

figure
plot(t,wEnv)

figure
plot(t,v)

%%

% % Gain from wEnv^3 to p
% kWT = pRated/(wRated^3);
% 
% 
% 
% t0 = 0;
% tf = 40;
% x0 = [2 0]';
% 
% %% Control Signal
% tu = [t0 tf/2 tf]';
% 
% % Reference for change in active power of the WT
% dP = [0 1e6 1e6]';
% 
% % Reference for derating the WT (this in principle changes the rated power)
% dPlim = [1e6 1e6 1e6]';
% 
% % Reference to reactive power
% qref = [1 0 2]';
% 
% u = [dP dPlim qref];
% 
% %%
% % Environmental disturbance
% wEnv = 12*sin(0.3*(0:1:tf)')*0+11;
% twEnv = (0:1:tf)';
% 
% options = '';
% 
% if(sum(dPlim<0)>0)
%     error('The reference for derated power can only contain nonnegative numbers.');
% end
% 
% % [t,x] = ode45(@wtVectorField,[t0 tf],x0,options,[tu u],[twEnv wEnv],tauEnv,tauWT,kWT,pRated,wMin,wMax);
% % 
% % figure
% % plot(t,x(:,2))
% % figure
% % plot(t,x(:,1))
% 
% %% Discrete time
% clear p dPlim dP
% dP = 0 % [W]
% dPlim = 0 % [W]
% PF = 1;
% v = 1;
% Qmode = 1;
% qRef = 1;
% wEnv = 11 % [m/s]
% if(Ts<0.5 || Ts>5)
%     error('This model is applicable for sampling times in the interval [0.5 5] min')
% end
% if(dPlim<0)
%     error('dPlim must be nonnegative.')
% end
% 
% %% Calculate Active Power Output
% dP = min(dP,dPlim); % dP is limited by dPlim
% pW = kWT*wEnv^3;    % [W] Available power
% vRef = 1;
% % Calculate gain from wind disturbance to power output
% if(wEnv<wMin || wEnv>wMax) % Out of normal operating range [wMin wMax] (in reality these should be 10 min mean values)
%     p = 0;
% elseif(pW>=pRated+dP-dPlim) % In full load operation
%     p = pRated+dP-dPlim;
% else   
%     dP = min(dP,0); % In partial load operation up-regulation is not possible in partial load operation
%     p = max(0,pW+dP); % The output power must be positive
% end
% %%
% qRef = qControl(Qmode,PF,qRef,vRef,v,alpha)
% %% Calculate Reactive Power Output
% qMax = sqrt(sMax^2-p^2); % [Var] Upper limit on reactive power
% 
% if(abs(qRef)>qMax)
%     q = sign(qRef)*qMax;
% else
%     q = qRef;
% end