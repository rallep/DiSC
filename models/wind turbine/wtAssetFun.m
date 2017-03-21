function [p,q] = wtAsset(sBase,vBase,v,w,param,dP,dPlim,Qmode,PF,qRef,vRef,qFun)
% wtAsset is a discrete time model of a wind turbine.
%
%   Input:
%       - sBase is the complex power base [VA].
%       - vBase is the voltage base [V].
%       - v is the bus voltage [PU].
%       - w is the wind speed [m/s].
%       - param is a struct with wind turbine parameters:
%           - param.pRated is the rated power [W].
%           - param.sMax is the maximum apperent power [VA].
%           - param.wMin is the cut-in wind speed [m/s].
%           - param.wMax is the cut-out wind speed [m/s].
%           - param.wRated is the rated wind speed [m/s].
%           - param.kWT is the maximum gain from w^3 to p.
%       - dP is the reference to active power change [PU].
%       - dPlim is the reference to derated power [PU].
%       - Qmode determines the reactive power control:
%           - Qmode=0: Power factor is kept constantly equal to PF.
%           - Qmode=1: Reactive power is to qRef [PU].
%           - Qmode=2: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
%       - vRef is the voltage reference [PU].
%
%   Output:
%       - p is the active power produced by the wind turbine [PU].
%       - q is the reactive power produced by the wind turbine [PU].

%% Test Validity of Input Arguments
if(dPlim<0)
    error('dPlim must be nonnegative.')
end

%% Wind Turbine Parameters
pRated = param.pRated;
sMax = param.sMax;
wMin = param.wMin;
wMax = param.wMax;
kWT = param.kWT;

%% Convert References to Physical Quantities
dP = dP*sBase;
dPlim = dPlim*sBase;
qRef = qRef*sBase;
v = v*vBase;
vRef = vRef*vBase;

%% Calculate Active Power Output
pW = kWT*w^3;    % [W] Available power

dP = min(dP,dPlim); % [W] dP is limited by dPlim

% Calculate gain from wind disturbance to power output
if(w<wMin || w>wMax) % Out of normal operating range [wMin wMax] (in reality these should be 10 min mean values)
    p = 0;
elseif(pW>=pRated+dP-dPlim) % In full load operation
    p = pRated+dP-dPlim;
else   
    dP = min(dP,0); % In partial load operation up-regulation is not possible in partial load operation
    p = max(0,pW+dP); % The output power must be positive
end

%% Reactive Power Control
if(exist('qFun')==0)
    qRef = qControl(vBase,v,p,Qmode,PF,qRef,vRef,sMax);
else
    qRef = qControl(vBase,v,p,Qmode,PF,qRef,vRef,sMax,qFun);
end
%% Calculate Reactive Power Output
qMax = sqrt(sMax^2-p^2); % [Var] Upper limit on reactive power
if(abs(qRef)>qMax)
    q = sign(qRef)*qMax;
else
    q = qRef;
end

%% Normalize Output to Per Unit (PU)
p = p/sBase;
q = q/sBase;