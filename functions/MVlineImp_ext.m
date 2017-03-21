function Z = MVlineImp_ext(oltc)
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
% 23-02-2015, R. Pedersen, Aalborg University. Notes: Extended the medium
% voltage grid line impedance function (MVlineImp.m) with an additional
% feeder. 

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
Rl10 = Rl1; Xl10 = Xl1; ll10 = ll1;
Rl11 = Rl2; Xl11 = Xl2; ll11 = 112;
Rl12 = Rl3; Xl12 = Xl3; ll12 = ll3;
Rl13 = Rl4; Xl13 = Xl4; ll13 = ll4;
Rl14 = Rl6; Xl14 = Xl6; ll14 = ll6;
Rl15 = Rl7; Xl15 = Xl7; ll15 = ll7;
Rl16 = Rl8; Xl16 = Xl8; ll16 = ll8;


% Trafo parameters
if oltc == true
    Rtrafo1 = inf;
    Xtrafo1 = inf;
else
    Rtrafo1 = 3;
    Xtrafo1 = 13;
    Rtrafo2 = 3;
    Xtrafo2 = 13;
end

% Format:
%     [from       to      R       X       l]
Z =   [1          2       Rtrafo1 Xtrafo1  1;       
       2          12      Rl1     Xl1      5;
       2          3       Rl1     Xl1     ll1;
       3          5       Rl2     Xl2     ll2;
       5          4       Rl4     Xl4     ll4;
       5          6       Rl3     Xl3     ll3;
       6          7       Rl6     Xl6     ll6;
       6          8       Rl5     Xl5     ll5;
       8          10      Rl7     Xl7     ll7;
       10         9       Rl8     Xl8     ll8;
       10         11      Rl9     Xl9     ll9;
       12         13      Rl10    Xl10    ll10;
       13         14      Rl11    Xl11    ll11;
       14         15      Rl12    Xl12    ll12;
       15         16      Rl13    Xl13    ll12;
       15         17      Rl14    Xl14    ll14;
       17         18      Rl15    Xl15    ll15;
       18         19      Rl16    Xl16    ll16];
end

