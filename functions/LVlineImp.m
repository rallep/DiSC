function Z = LVlineImp(oltc)
% Function for setting up the low voltage grid impedance matrix
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

% Generic parameters
R = 0.208;
X = 0.052;
l = 0.5;

% Trafo parameters
if oltc == true
    Rtrafo = inf;
    Xtrafo = inf;
else
    Rtrafo = 0.004;
    Xtrafo = 0.04;
end

% Format:
%     [from       to      R       X       l]
Z =   [1          2       Rtrafo  Xtrafo  1;       
       2          3       R       X       l;
       2          18      R       X       l;
       2          33      R       X       l;
       3          4       R       X       l;
       3          5       R       X       l;
       4          6       R       X       l;
       4          7       R       X       l;
       5          8       R       X       l;
       6          9       R       X       l;
       6          10      R       X       l;
       8          11      R       X       l;    
       9          12      R       X       l;
       10         13      R       X       l;
       11         14      R       X       l;
       12         15      R       X       l;
       13         16      R       X       l;
       14         17      R       X       l;
       18         19      R       X       l;
       18         20      R       X       l;
       18         33      R       X       l;       
       19         21      R       X       l;
       20         22      R       X       l;
       22         23      R       X       l;
       22         24      R       X       l;
       23         25      R       X       l;
       24         26      R       X       l;
       25         27      R       X       l;
       26         28      R       X       l;
       27         29      R       X       l;
       28         30      R       X       l;
       29         31      R       X       l;
       30         32      R       X       l;
       33         34      R       X       l;
       33         35      R       X       l;
       34         36      R       X       l;
       35         37      R       X       l;
       36         38      R       X       l;
       38         39      R       X       l;
       39         40      R       X       l;
       39         41      R       X       l;
       40         42      R       X       l;];


end

