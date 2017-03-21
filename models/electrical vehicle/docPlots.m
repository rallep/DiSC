close all
%% Trips in progress (car driver trips)
% Time
time = 0:1:23;

% Week day
TiP_WD = [4 2 2  2  5  21  54  160  273  174  146  152  151  141  155  205  213  242  188  122  74  54  37  19];
% Saturday
TiP_Sat = [9  4  3  3  4  11  22  48  93  149  202  215  207  182  178  159  150  136 115 87  57  37  37  29];
% Sunday
TiP_Sun = [ 12  4  2  2  3  6  16  29  48  96  154  172  173  158  146  142  129  109  95  75  52  33  26  16 ];

figure
plot(time,TiP_WD,'r')
hold on
plot(time,TiP_Sat,'color',[0 0.5 0])
plot(time,TiP_Sun,'b')
xlim([0 23])
title('Trips in Progress')
%%
CommuteStart = [0 0 0 0 0 3 7 16 15 4 2 2 3 3 3 4 9 15 7 2 2 1 1 1]/100;
figure
plot(time,CommuteStart)
xlim([0 23])
title('Probability of Commute Start')

ToCommute =   [0 0 0 0 0 3 7 16 15 4 2 2 1 0 0 0 0 0 0 0 0 0 0 0]/50;
FromCommute = [0 0 0 0 0 0 0 0 0 0 0 0 2 3 3 4 9 15 7 2 2 1 1 1]/50;

figure
plot(time,ToCommute)
xlim([0 23])
title('Probability of Driving to Work')

figure
plot(time,FromCommute)
xlim([0 23])
title('Probability of Driving from Work')

%%
ShoppingStart = [0 0 0 0 0 0 0 1 3 8 12 13 11 9.5 9.5 8 7 6 5 4 2 1 0 0]/100;
figure
plot(time,ShoppingStart)
xlim([0 23])
title('Probability of Shopping Start')

ToShopping =    [0 0 0 0 0 0 0 1 2 5 7 7 6 4.5 4.5 3 3 2 2 2 1 0 0 0]/50;
FromShopping =  [0 0 0 0 0 0 0 0 1 3 5 6 5 5 5 5 4 4 3 2 1 1 0 0]/50;

figure
plot(time,ToShopping)
xlim([0 23])
title('Probability of Driving to Shop')

figure
plot(time,FromShopping)
xlim([0 23])
title('Probability of Driving from Shop')