clc
clear
close all

%%
%someCarPlots

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

%% Modified data used in Simulation
ToEducation = [0 0 0 0 0 0 0 7 40 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
FromEducation = [0 0 0 0 0 0 0 0 0 0 0 1 2 2 3 32 5 2 1 0 0 0 0 0];

ToShopping =    [0 0 0 0 0 0 0 1 2 5 7 7 6 4.5 4.5 3 3 2 2 2 1 0 0 0];
FromShopping =  [0 0 0 0 0 0 0 0 1 3 5 6 5 5 5 5 4 4 3 2 1 1 0 0];

ToCommute =   [0 0 0 0 0 3 7 16 15 4 2 2 1 0 0 0 0 0 0 0 0 0 0 0];
FromCommute = [0 0 0 0 0 0 0 0 0 0 0 0 2 3 3 4 9 15 7 2 2 1 1 1];

BusinessStart = [0 0 0 0 0.5 1 3 7 10 8 7 7 7 7 7 8 9 8 4 2 2 1 1 0.5];
ToBusiness =   [0 0 0 0 0.5 1 3 7 10 8 6 5 3.5 2 2 1 1 0 0 0 0 0 0 0];
FromBusiness = [0 0 0 0 0   0 0 0 0  0 1 2 3.5 5 5 7 8 8 4 2 2 1 1 0.5];

ToEscortEducation =   [0 0 0 0 0 0 0 3 36 7 0 2 2 0 0 0 0 0 0 0 0 0 0 0];
FromEscortEducation = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 9 33 4 2 1 0 0 0 0 0];

ToOtherPersonal =   [0 0 0 0 0 0 1 4 8 8 7 7 6 4 2 1 1 1 0 0 0 0 0 0];
FromOtherPersonal = [0 0 0 0 0 0 0 0 0 0 1 1 2 3 5 7 8 7 6 4 2 2 1 1];

ToVisit =   [0 0 0 0 0 0 1 1 2 4 5 6 6 5 4 4 3 3 3 2 1 0 0 0];
FromVisit = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 5 6 7 8 5 5 4 3];

ToHoliday =    [0 0 0 0 0 0 2 4 5 7 8 6 3 3 3 2 2 2 2 1 0 0 0 0];
FromHoliday =  [0 0 0 0 0 0 0 0 0 0 1 3 3 4 5 7 7 5 5 4 3 2 1 1];

%%
[ToCommutingCDF, FromCommutingCDF] = getToFromCDFs(ToCommute,FromCommute);
[ToBusinessCDF, FromBusinessCDF] = getToFromCDFs(ToBusiness,FromBusiness);
[ToEducationCDF, FromEducationCDF] = getToFromCDFs(ToEducation,FromEducation);
[ToEscortEducationCDF, FromEscortEducationCDF] = getToFromCDFs(ToEscortEducation,FromEscortEducation);
[ToShoppingCDF, FromShoppingCDF] = getToFromCDFs(ToShopping,FromShopping);
[ToPersonalBusinessOtherEscortCDF, FromPersonalBusinessOtherEscortCDF] = getToFromCDFs(ToOtherPersonal,FromOtherPersonal);
[ToVisitFriendsSportCDF, FromVisitFriendsSportCDF] = getToFromCDFs(ToVisit,FromVisit);
[ToHolidayCDF, FromHolidayCDF] = getToFromCDFs(ToHoliday,FromHoliday);

%% 
daysPerYear = 365;
individualDayOfWeekPerYear = daysPerYear/7;
%weekendDaysPerYear = daysPerYear-weekDaysPerYear;

%% Round-trip rates as car driver by day of weak

% Probability of being driver of car/van for trip (based on Table NTS0409)
CommutingCarDriverProb = 84/146;
BusinessCarDriverProb = 22/31;
EducationOrEduEscortCarDriverProb = 26/116;
ShoppingCarDriverProb = 85/189;
OtherEscortCarDriverProb = 50/87;
PersonalBusinessCarDriverProb = 41/94;
LeisureCarDriverProb = 93/248;

% Trip rates by day of weak from Table NTS0504 (We only consider round-trips; hence, the number of trips from the data is divided by two)
CommutingPerYear = [26  28 28 28 26 9 5]/2; 
BusinessPerYear = [5  6 6 6 5 2 1]/2; 
EducationPerYear = [12  12 13 12 11 0 0]/2; 
EscortEducationPerYear = [9 9 9 9 9 0 0]/2; 
ShoppingPerYear = [25 24 25 26 29 42 23]/2; 
OtherEscortPerYear = [13 14 14 13 14 14 10]/2; 
PersonalBusinessPerYear = [15 16 16 16 15 9 11]/2; 
VisitFriendsPrivateHomePerYear = [12 12 13 13 14 20 21]/2; 
VisitFriendsElsewherePerYear = [4 4 5 6 8 10 9]/2; 
SportPerYear = [8 9 10 9 8 12 8]/2; 
HolidayPerYear = [5 5 6 6 6 7 8]/2; 

% Round-trip rates as car driver (note that probability of car by purpose is only available for groups of purposes)
CommutingPerYear = CommutingPerYear*CommutingCarDriverProb;
BusinessPerYear = BusinessPerYear*BusinessCarDriverProb; 
EducationPerYear = EducationPerYear*EducationOrEduEscortCarDriverProb; 
EscortEducationPerYear = EscortEducationPerYear*EducationOrEduEscortCarDriverProb; 
ShoppingPerYear = ShoppingPerYear*ShoppingCarDriverProb; 
OtherEscortPerYear = OtherEscortPerYear*OtherEscortCarDriverProb; 
PersonalBusinessPerYear = PersonalBusinessPerYear*PersonalBusinessCarDriverProb; 
VisitFriendsPrivateHomePerYear = VisitFriendsPrivateHomePerYear*LeisureCarDriverProb; 
VisitFriendsElsewherePerYear = VisitFriendsElsewherePerYear*LeisureCarDriverProb; 
SportPerYear = SportPerYear*LeisureCarDriverProb; 
HolidayPerYear = HolidayPerYear*LeisureCarDriverProb; 
%% Average trip length in miles Table NTS0405

% We only consider round-trips; hence, the length from the data is multiplied by two
CommutingLength = 9*2;
BusinessLength = 19.6*2;
EducationLength = 3.4*2;
EscortEducationLength = 2.4*2;
ShoppingLength = 4.4*2;
OtherEscortLength = 5.5*2;
PersonalBusinessLength = 5.2*2;
VisitFriendsPrivateHomeLength = 10.2*2; 
VisitFriendsElsewhereLength = 6.4*2;
% Sport and Entertainment separated in NTS0405; hence, average is used
SportLength = (7.6+6.4)/2*2;
% Holiday and Day trip separated in NTS0405; hence, average is used
HolidayLength = (43.3+13.1)/2*2;


%% Inputs
dayNumber = 2;

%% Parameters
avgMileage = 2; % pct battery/mile

%% Simulation
dSOC = zeros(1,length(FromBusiness));
CommutingAway = zeros(1,length(FromBusiness));
if(rand<CommutingPerYear(dayNumber)/individualDayOfWeekPerYear)
    [CommutingAway,CommutingStart] = getAwayPeriod(ToCommutingCDF,FromCommutingCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(CommutingStart) = dSOC(CommutingStart)-CommutingLength*avgMileage;
end
BusinessAway = zeros(1,length(FromBusiness));
if(rand<BusinessPerYear(dayNumber)/individualDayOfWeekPerYear)
    [BusinessAway,BusinessStart] = getAwayPeriod(ToBusinessCDF,FromBusinessCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(BusinessStart) = dSOC(BusinessStart)-BusinessLength*avgMileage;
end
EducationAway = zeros(1,length(FromBusiness));
if(rand<EducationPerYear(dayNumber)/individualDayOfWeekPerYear)
    [EducationAway,EducationStart] = getAwayPeriod(ToEducationCDF,FromEducationCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(EducationStart) = dSOC(EducationStart)-EducationLength*avgMileage;
end
EscortEducationAway = zeros(1,length(FromBusiness));
if(rand<EscortEducationPerYear(dayNumber)/individualDayOfWeekPerYear)
    [EscortEducationAway,EscortEducationStart] = getAwayPeriod(ToEscortEducationCDF,FromEscortEducationCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(EscortEducationStart) = dSOC(EscortEducationStart)-EscortEducationLength*avgMileage;
end
ShoppingAway = zeros(1,length(FromBusiness));
if(rand<ShoppingPerYear(dayNumber)/individualDayOfWeekPerYear)
    [ShoppingAway,ShoppingStart] = getAwayPeriod(ToShoppingCDF,FromShoppingCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(ShoppingStart) = dSOC(ShoppingStart)-ShoppingLength*avgMileage;
end
OtherEscortAway = zeros(1,length(FromBusiness));
if(rand<OtherEscortPerYear(dayNumber)/individualDayOfWeekPerYear)
    [OtherEscortAway,OtherEscortStart] = getAwayPeriod(ToPersonalBusinessOtherEscortCDF,FromPersonalBusinessOtherEscortCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(OtherEscortStart) = dSOC(OtherEscortStart)-OtherEscortLength*avgMileage;
end
PersonalBusinessAway = zeros(1,length(FromBusiness));
if(rand<PersonalBusinessPerYear(dayNumber)/individualDayOfWeekPerYear)
    [PersonalBusinessAway,PersonalBusinessStart] = getAwayPeriod(ToPersonalBusinessOtherEscortCDF,FromPersonalBusinessOtherEscortCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(PersonalBusinessStart) = dSOC(PersonalBusinessStart)-PersonalBusinessLength*avgMileage;
end
VisitFriendsPrivateHomeAway = zeros(1,length(FromBusiness));
if(rand<VisitFriendsPrivateHomePerYear(dayNumber)/individualDayOfWeekPerYear)
    [VisitFriendsPrivateHomeAway,VisitFriendsPrivateHomeStart] = getAwayPeriod(ToVisitFriendsSportCDF,FromVisitFriendsSportCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(VisitFriendsPrivateHomeStart) = dSOC(VisitFriendsPrivateHomeStart)-VisitFriendsPrivateHomeLength*avgMileage;
end
VisitFriendsElsewhereAway = zeros(1,length(FromBusiness));
if(rand<VisitFriendsElsewherePerYear(dayNumber)/individualDayOfWeekPerYear)
    [VisitFriendsElsewhereAway,VisitFriendsElsewhereStart] = getAwayPeriod(ToVisitFriendsSportCDF,FromVisitFriendsSportCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(VisitFriendsElsewhereStart) = dSOC(VisitFriendsElsewhereStart)-VisitFriendsElsewhereLength*avgMileage;
end
SportAway = zeros(1,length(FromBusiness));
if(rand<SportPerYear(dayNumber)/individualDayOfWeekPerYear)
    [SportAway,SportStart] = getAwayPeriod(ToVisitFriendsSportCDF,FromVisitFriendsSportCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(SportStart) = dSOC(SportStart)-SportLength*avgMileage;
end
HolidayAway = zeros(1,length(FromBusiness));
if(rand<HolidayPerYear(dayNumber)/individualDayOfWeekPerYear)
    [HolidayAway,HolidayStart] = getAwayPeriod(ToHolidayCDF,FromHolidayCDF);
    % Calculate change in state of charge (SoC), where negative means discharge
    dSOC(HolidayStart) = dSOC(HolidayStart)-HolidayLength*avgMileage;
end
Away = CommutingAway|BusinessAway|EducationAway|EscortEducationAway|ShoppingAway|OtherEscortAway|PersonalBusinessAway|VisitFriendsPrivateHomeAway|VisitFriendsElsewhereAway|SportAway|HolidayAway;

figure
subplot(2,1,1)
stairs(Away)
hold on
stairs(CommutingAway,'--k')
subplot(2,1,2)
stairs(dSOC,'r')
