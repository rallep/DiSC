function Z = MVlineImp(oltc)
% Function for setting up the medium voltage grid impedance matrix
%   The impedance matrix, can be set to run with an on-load tap changing
%   transformer. If this is the case, the resistance and reactance of the
%   impedance matrix is set to infinity. Then when the OLTC is added to the
%   simulation, the impedances are changed to the correct value.
%
% Input:
%   - oltc, is a boolean input where: false = no OLTC, and true = with OLTC
%
% Output:
%   - Z, is the impedance matrix
%
% Revision:
% 10-09-2014, R. Pedersen, Aalborg University. Notes:

% Line parameters
Rl1 = 0.1;  Xl1 = 0.1;  ll1 = 1;
Rl2 = 0.1;  Xl2 = 0.1;  ll2 = 1;
Rl3 = 0.13; Xl3 = 0.09; ll3 = 5;
Rl4 = 0.13; Xl4 = 0.09; ll4 = 10;
Rl5 = 0.13; Xl5 = 0.09; ll5 = 5;
Rl6 = 0.32; Xl6 = 0.15; ll6 = 5;
Rl7 = 0.13; Xl7 = 0.09; ll7 = 5;
Rl8 = 0.32; Xl8 = 0.15; ll8 = 5;
Rl9 = 0.13; Xl9 = 0.09; ll9 = 10;



% Trafo parameters
if oltc == true
    Rtrafo1 = inf;
    Xtrafo1 = inf;
    Rtrafo2 = inf;
    Xtrafo2 = inf;
else
    Rtrafo1 = 0.3;
    Xtrafo1 = 0.13;
    Rtrafo2 = 0.3;
    Xtrafo2 = 0.13;
end

% Format:
%     [from       to      R       X       l]
Z =   [1          2       Rtrafo1 Xtrafo1  1;       
       1          3       Rtrafo2 Xtrafo2  1;
       2          4       Rl1     Xl1     ll1;
       4          6       Rl2     Xl2     ll2;
       6          5       Rl4     Xl4     ll4;
       6          7       Rl3     Xl3     ll3;
       7          8       Rl6     Xl6     ll6;
       7          9       Rl5     Xl5     ll5;
       9          11      Rl7     Xl7     ll7;
       11         10      Rl8     Xl8     ll8;
       11         12      Rl9     Xl9     ll9];


end

