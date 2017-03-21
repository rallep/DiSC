function [ siSys ] = setupSolarIrradiance( Ts, num, latitude, pAirInt, transInt )
%setupSolarIrradiance returns a set containing the specified number of
% solar irradiance objects.  
%   The set of solar irradiance objects can be used together with a fleet
%   of PV systems. The set of solar irradiance objects could match the set
%   of PV objects or be smaller, allowing some PV systems to experiance the
%   same irradiance.
%
% Inputs:
%   - Ts [sec], is the systems sampling time.
%   - num [-], is the number of solar irradiance objects.
%   - latitude [deg], is the latitude of the stidiet grid.
%   - pAirInt [kPa], is the air pressure at the studiet site.
%   - transInt [-], is the transmittance level at the site.
%
% Outputs:
%   - siSys [-], is the set of solar irradiance objects.
%
% R. Pedersen, 23-6-2014, Aalborg University

% Error checking

% set Parameter Struct
param.Ts = Ts;
param.lat = latitude;

% Create set of solar irradiance objects
for i=1:num
    % Create parameter struct
    param.p = ceil(pAirInt(1)+(pAirInt(2)-pAirInt(1))*rand(1));
    param.t = transInt(1)+(transInt(2)-transInt(1))*rand(1);
    % Create object
    siSys(i) = solarIrradiance(param);
end


end

