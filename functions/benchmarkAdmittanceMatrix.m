function [ Y ] = benchmarkAdmittanceMatrix( grid, oltcMV, oltcLV , zBaseMV, zBaseLV,onPU)
%Forms the grid admittance matrix.
%   Function is used to form the grid admittance matrix used for load flow
%   analysis. The admittance matrix is put on a per unit base, based on the
%   base impedances of the medium and lov voltage grids. 
%
%   Inputs
%       - grid, string decides which grid impedance matric is outputted.
%         'MV' = Only MV grid
%         'LV' = Only LV grid
%         'both' = Both
%       - oltcMV, indicates if there should be included an on-load tap
%         changing transformer on the medium voltage level. Note, that
%         there are two transformers on the MV level.
%         false = no OLTC
%         true = with OLTC
%       - oltcLV, the same as above but for the low voltage level.
%       - zBaseMV, is the impedance base of the medium voltage grid.
%       - zBaseLV, is the impedance base of the LV grid.
%       - onPU, indicates if the system is to be simulated on the per unit
%         basis. true = no, false = yes
%
%   Output
%       - Y, is the benchmark grids admittance matrix.
%
% Revision:
% 10-09-2014, R. Pedersen, Aalborg University. Notes:

% Error checking
if nargin ~= 6
    error('myApp:argChk', 'Wrong number of input arguments')
end

% temp = strcmp(grid,'MV')+strcmp(grid,'LV')+strcmp(grid,'both');
% if temp ~= 1
%     error('myApp:valChk', 'Input is wrong. Accepted strings are: ''MV'',''LV'', ''both''')
% end

% Form admittance matrix
% MV grid
if strcmp(grid,'MV')
    Z = MVlineImp(oltcMV);
    if onPU == true
        Z(:,3:4) = Z(:,3:4)./zBaseMV;
    end
    Y = ztoybus(Z);
% LV grid
elseif strcmp(grid,'LV')
    Z = LVlineImp(oltcLV);
    if onPU == true
        Z(:,3:4) = Z(:,3:4)./zBaseLV;
    end
    Y = ztoybus(Z);
% Both collected into one
elseif strcmp(grid,'both')
    % Load impedance matrices
    Zmv = MVlineImp(oltcMV);
    Zlv = LVlineImp(oltcLV);
    % cat them
    Zlv(1,1) = Zlv(1,1)+9;
    Zlv(1,2) = Zlv(1,2)+11;
    Zlv(2:end,1:2) = Zlv(2:end,1:2)+11; 
    Z = [Zmv; Zlv];
    if onPU == true
        Z(1:12,3:4) = Z(1:12,3:4)./zBaseMV;
        Z(13:end,3:4) = Z(13:end,3:4)./zBaseMV;
    end
    Y = ztoybus(Z);
end

    
        

end

