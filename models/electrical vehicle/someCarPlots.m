%%
time = 0:1:23;

WD = [4 2 2 2 5 21 54 160 273 174 146 152 151 141 155 205 213 242 188 122 74 54 37 19];

plot(time,WD)

%% Based on Table NTS0503
CommuteStart = [0 0 0 0 0 3 7 16 15 4 2 2 3 3 3 4 9 15 7 2 2 1 1 1];
BusinessStart = [0.1 0.1 0 0.1 0.3 1 3 7 10 8 7 7 7 7 7 8 9 8 4 2 2 1 1 0.5];
AllStart = [0 0 0 0 0 1 2 5 12 6 6 6 6 6 6 11 8 8 6 4 3 2 1 1];

EducationStart = [0 0 0 0 0 0 0 7 40 2 1 1 2 2 3 32 5 2 1 0 0 0 0 0];
EscortEducationStart = [0 0 0 0 0 0 0 3 36 7 0 2 2 1 9 33 4 2 1 0 0 0 0 0];
ShoppingStart = [0 0 0 0 0 0 0 1 3 8 12 13 11 9.5 9.5 8 7 6 5 4 2 1 0 0];
OtherPersonalStart = [0 0 0 0 0 0 1 4 8 8 8 8 8 7 7 8 9 8 6 4 2 2 1 1];
VisitStart = [0 0 0 0 0 0 1 1 2 4 5 6 6 6 6 8 8 9 10 10 6 5 4 3];
HolidayStart = [0 0 0 0 0 0 2 4 5 7 9 8 6 7 8 9 9 7 7 5 3 2 1 1];
figure
plot(time,CommuteStart)
hold on
figure
plot(time,BusinessStart)
plot(time,AllStart,'r')

% figure
% plot(time,EducationStart)

figure
plot(time,ShoppingStart)
close all
%% Modification of data
ToEducation = [0 0 0 0 0 0 0 7 40 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
FromEducation = [0 0 0 0 0 0 0 0 0 0 0 1 2 2 3 32 5 2 1 0 0 0 0 0];

figure
stairs(time,ToEducation, 'g')
hold on
stairs(time,FromEducation,'r')
title('Education')

%%
ToShopping =    [0 0 0 0 0 0 0 1 2 5 7 7 6 4.5 4.5 3 3 2 2 2 1 0 0 0];
FromShopping =  [0 0 0 0 0 0 0 0 1 3 5 6 5 5 5 5 4 4 3 2 1 1 0 0];
figure
stairs(time,ToShopping, 'g')
hold on
stairs(time,FromShopping,'r')
title('Shopping')
%%

ToCommute =   [0 0 0 0 0 3 7 16 15 4 2 2 1 0 0 0 0 0 0 0 0 0 0 0];
FromCommute = [0 0 0 0 0 0 0 0 0 0 0 0 2 3 3 4 9 15 7 2 2 1 1 1];

figure
stairs(time,ToCommute, 'g')
hold on
stairs(time,FromCommute,'r')
title('Commute')

%%

BusinessStart = [0 0 0 0 0.5 1 3 7 10 8 7 7 7 7 7 8 9 8 4 2 2 1 1 0.5];
ToBusiness =   [0 0 0 0 0.5 1 3 7 10 8 6 5 3.5 2 2 1 1 0 0 0 0 0 0 0];
FromBusiness = [0 0 0 0 0   0 0 0 0  0 1 2 3.5 5 5 7 8 8 4 2 2 1 1 0.5];
figure
stairs(time,ToBusiness, 'g')
hold on
stairs(time,FromBusiness,'r')
title('Business')

%%
ToEscortEducation =   [0 0 0 0 0 0 0 3 36 7 0 2 2 0 0 0 0 0 0 0 0 0 0 0];
FromEscortEducation = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 9 33 4 2 1 0 0 0 0 0];

figure
stairs(time,ToEscortEducation, 'g')
hold on
stairs(time,FromEscortEducation,'r')
title('Escort Education')
%%
ToOtherPersonal =   [0 0 0 0 0 0 1 4 8 8 7 7 6 4 2 1 1 1 0 0 0 0 0 0];
FromOtherPersonal = [0 0 0 0 0 0 0 0 0 0 1 1 2 3 5 7 8 7 6 4 2 2 1 1];
figure
stairs(time,ToOtherPersonal, 'g')
hold on
stairs(time,FromOtherPersonal,'r')
title('Other Personal')
%%
ToVisit =   [0 0 0 0 0 0 1 1 2 4 5 6 6 5 4 4 3 3 3 2 1 0 0 0];
FromVisit = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 5 6 7 8 5 5 4 3];

figure
stairs(time,ToVisit, 'g')
hold on
stairs(time,FromVisit,'r')
title('Visit')
%%

ToHoliday =    [0 0 0 0 0 0 2 4 5 7 8 6 3 3 3 2 2 2 2 1 0 0 0 0];
FromHoliday =  [0 0 0 0 0 0 0 0 0 0 1 3 3 4 5 7 7 5 5 4 3 2 1 1];

figure
stairs(time,ToHoliday, 'g')
hold on
stairs(time,FromHoliday,'r')
title('Holiday')