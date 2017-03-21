clc
clear 
close all

Svdh = [0.5 0.8 1.1 1.5 2.5 4 4.8 2.8 1.5 1.9 1.6 0.9 0.7...
                0.45 0.4 0.35 0.3 0.25 0.22 0.2];
x = logspace(-3,0,length(Svdh));   % Converted to [rad/s] from [cycles/hour] 

Svdh = Svdh./x

i=1;
A1 = 2/pi*sqrt(1/2*(Svdh(i)+Svdh(i+1))*(x(i+1)-x(i))) 