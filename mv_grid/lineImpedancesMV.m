% Script for constructiong the MV grid line impedance matrix
% Also the trafo impedance is included.

% Line parameters
Rl1 = 0.1; Xl1 = 0.1; ll1 = 1;
Rl2 = 0.1; Xl2 = 0.1; ll2 = 1;
Rl3 = 0.13; Xl3 = 0.09; ll3 = 5;
Rl4 = 0.13; Xl4 = 0.09; ll4 = 10;
Rl5 = 0.13; Xl5 = 0.09; ll5 = 5;
Rl6 = 0.32; Xl6 = 0.15; ll6 = 5;
Rl7 = 0.13; Xl7 = 0.09; ll7 = 5;
Rl8 = 0.32; Xl8 = 0.15; ll8 = 5;
Rl9 = 0.13; Xl9 = 0.09; ll9 = 10;


% Trafo parameters
Rtrafo = 3;
Xtrafo = 13;

% Format:
%     [from       to      R       X       l]
Zmv = [1          2       Rtrafo  Xtrafo  1;       
       1          3       Rtrafo  Xtrafo  1;
       2          4       Rl1     Xl1     ll1;
       4          6       Rl2     Xl2     ll2;
       6          5       Rl4     Xl4     ll4;
       6          7       Rl3     Xl3     ll3;
       7          8       Rl6     Xl6     ll6;
       7          9       Rl5     Xl5     ll5;
       9          11      Rl7     Xl7     ll7;
       11         10      Rl8     Xl8     ll8;
       11         12      Rl9     Xl9     ll9];
 
 clear Rtrafo Xtrafo