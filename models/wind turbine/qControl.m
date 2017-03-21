function qRef = qControl(vBase,v,p,mode,PF,qRef,vRef,sMax,qFun)
% qControl is a reactive power control.
%
%   Input:
%       - vBase is the voltage base [V].
%       - v is the bus voltage [PU].
%       - p is the active power [W].
%       - mode determines the reactive power control:
%           - mode=0: Power factor is kept constantly equal to PF.
%           - mode=1: Reactive power is to qRef [VAR].
%           - mode=2: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
%       - vRef is the voltage reference [PU].
%       - sMax is the maximum apperent power [VA].
%
%   Output:
%       - qRef is the reactive power reference [VAR].

% Calculate Reactive Power Reference
if(mode == 0)
    qRef = p*tan(acos(PF));
elseif(mode == 1)
    qRef = qRef;
elseif(mode == 2)
    if(exist('qFun')==0)
        qFun = @defaultQDroop;
    end
    qRef = feval(qFun,vBase,v,vRef,sMax);
else
    error('Invalid mode for reactive power control.');
end
