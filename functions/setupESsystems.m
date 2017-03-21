function [ esSys, busses ] = setupESsystems( Ts, randomPlacement, num, sBase, vBase, busInt, sMaxInt, pRatedIntMax, pRatedIntMin, eRatedInt )
%setupESsystems returns a set of ES system structs. 
%   Setup a portfolio of energy storage systems. The ES systems could be placed in
%   e.g., a low voltage distribution grid. This function relies on the
%   class esAsset. The class esAsset implements a single P controllable
%   ES system, with voltage droop control on active power.
%
%   The systems are randomly constructed according to the intervals given
%   by the inputs described below. 
%
% Inputs:
%   - randomPlacement [-], is a flag describing if the placement of ES
%     systems on the busses, should be random based on the busInt input or
%     if the busInt describes exactly which busses the ES systems should be
%     placed.
%   - num [-], is the number of ES systems wanted.
%   - sBase [VA], is the complex power base of the grid studiet.
%   - vBase [V], is the voltage base of the grid studiet.
%   - busInt [-], is the interval of busses the ES systems can be placed
%     on, e.g., [2 30] means that the systems can be placed on bus 2 to 30.
%     If the systems are not to be placed randomly this input describes the
%     exact placement of PV systems, i.e., busInt = [2 3 15 34 ...]. Thus,
%     length(busInt) = num.
%   - sMaxInt [VA], is the apparent power interval, e.g., [6e3 11e3].
%   - pRatedIntMax [W], is the rated power interval of consumption.
%   - pRatedIntMin [W], is the rated power output.
%   - eRated [Wh], is the rated energy storage level interval.
%
% Outputs:
%   - esSys [struct of structs], object containing the ES systems structs.
%   - busses [-], a ordered list describing which busses the ES systems are
%     placed. That is: pvSys(1) is placed on busses(1) and so forth. Note:
%     busses(1) might be any number in the busInt.
%
% R. Pedersen, 25-6-2014, Aalborg University

% Error checking
if randomPlacement == false && length(busInt)~= num
    error('You have chosen to place the ES systems manually, but the number of busses specified in busInt does not match the number of busses specified in num.')
end

% Set parameter struct
param.sBase = sBase;
param.vBase = vBase;
param.Ts = Ts;

for i=1:num
    % Create parameter struct
    param.pRatedMax = ceil(pRatedIntMax(1)+(pRatedIntMax(2)-pRatedIntMax(1))*rand(1));
    param.pRatedMin = ceil(pRatedIntMin(1)+(pRatedIntMin(2)-pRatedIntMin(1))*rand(1));
    param.sMax = ceil(sMaxInt(1)+(sMaxInt(2)-sMaxInt(1))*rand(1));
    param.eRated = ceil(eRatedInt(1)+(eRatedInt(2)-eRatedInt(1))*rand(1));
    if max(param.pRatedMax,abs(param.pRatedMin)) > param.sMax
        param.sMax = max(param.pRatedMax,abs(param.pRatedMin));
    end
    
    % Create PV system
    esSys(i) = esAsset(param);
end

% Place on busses
if randomPlacement == true
    busses = ceil(busInt(1)+(busInt(2)-busInt(1))*rand(num,1));
else
    busses = busInt;
end

end

