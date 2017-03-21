function qRef = defaultQDroop(vBase,v,vRef,sMax)
% defaultQDroop is a reactive power droop control.
%
%   Input:
%       - vBase is the voltage base [V].
%       - v is the bus voltage [PU].
%       - vRef is the voltage reference [PU].
%       - sMax is the maximum apperent power [VA].
%
%   Output:
%       - qRef is the reactive power reference [VAR].
vMax = 1.1;
vMin = 0.9;

k = 2/(vMax-vMin);
qRef = sMax*k*(v-vRef)/vBase;