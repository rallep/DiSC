classdef evAsset < handle
    %evAsset implements a simple electrical vehical mobility- and battery model.
    %   The active power consumption of the battery charging cycle can be
    %   controlled when ever the electrical vehical is available at 
    %   a household charging spot.
    %   To simulate when the electrical vehical is available data from the
    %   National Travel Survey from Great Britain is used.
    %   
    %   To model the battery a simple first order model of an energy
    %   storages is used. The batterys state of charges, when arriving at a
    %   charging station/spot, is dependent on the trip the car has conducted. 
    %
    %   
    % Revision:
    % 16-07-2014, C.E. Sloth and R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 3e6;                % [VA]. Complex power base
        vBase = 20e3;               % [V]. Voltage base
        pRated = 6e6;               % [W]. Rated output power
        sMax = 3e6;                 % [VA]. Maximum apparent power
        PF = 1;                     % [-]. Power factor 
        Ts = 60;                    % [s]. Sampling time
        onPU = false;               % [-]. Indicates if the electrical vehicle is used to simulate a system on per unit basis. true = yes, false = no
        
        % Control parameters
        pMode = 0;                  % [-]. Indicates which control mode the charging of the EV is in
        qMode = 0;                  % [-]. Reactive power control mode
        
        % Parameters battery
        avgMileage = 2;             % [% battery/mile]. Discharge parameter
        eRated = 65e3*60*60;        % [J]. Battery capacity
        eOld = 0;                   % [J]. Place holder for battery state
        dSOC = 0;                   % [%]. Change in battery state of charge
        a = 1;                      % [-]. Drain rate of battery
        eta = 1;                    % [-]. Efficiency of battery
        pOld = 0;                   % [W]. Power consumption last sample
        pRate = 20;                 % [W/s]. Constraint on power rate of change
        
        % Mobility model
        away                        % [-]. Indicates if vehicle is home and ready for charging, 0 = at home and 1 = away
        dayOld = -1;                % [-]. Used to check if day has changed
        h = 2;                      % [-]. Sample of the day
        temp = 0;                   % [-]. Flag used for the mobility model
        
        % Modified data based on Table NTS0503
        toEducation = [0 0 0 0 0 0 0 7 40 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
        fromEducation = [0 0 0 0 0 0 0 0 0 0 0 1 2 2 3 32 5 2 1 0 0 0 0 0];
        toShopping =    [0 0 0 0 0 0 0 1 2 5 7 7 6 4.5 4.5 3 3 2 2 2 1 0 0 0];
        fromShopping =  [0 0 0 0 0 0 0 0 1 3 5 6 5 5 5 5 4 4 3 2 1 1 0 0];
        toCommute =   [0 0 0 0 0 3 7 16 15 4 2 2 1 0 0 0 0 0 0 0 0 0 0 0];
        fromCommute = [0 0 0 0 0 0 0 0 0 0 0 0 2 3 3 4 9 15 7 2 2 1 1 1];
        businessStart = [0 0 0 0 0.5 1 3 7 10 8 7 7 7 7 7 8 9 8 4 2 2 1 1 0.5];
        toBusiness =   [0 0 0 0 0.5 1 3 7 10 8 6 5 3.5 2 2 1 1 0 0 0 0 0 0 0];
        fromBusiness = [0 0 0 0 0   0 0 0 0  0 1 2 3.5 5 5 7 8 8 4 2 2 1 1 0.5];
        toEscortEducation =   [0 0 0 0 0 0 0 3 36 7 0 2 2 0 0 0 0 0 0 0 0 0 0 0];
        fromEscortEducation = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 9 33 4 2 1 0 0 0 0 0];
        toOtherPersonal =   [0 0 0 0 0 0 1 4 8 8 7 7 6 4 2 1 1 1 0 0 0 0 0 0];
        fromOtherPersonal = [0 0 0 0 0 0 0 0 0 0 1 1 2 3 5 7 8 7 6 4 2 2 1 1];
        toVisit =   [0 0 0 0 0 0 1 1 2 4 5 6 6 5 4 4 3 3 3 2 1 0 0 0];
        fromVisit = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 4 5 6 7 8 5 5 4 3];
        toHoliday =    [0 0 0 0 0 0 2 4 5 7 8 6 3 3 3 2 2 2 2 1 0 0 0 0];
        fromHoliday =  [0 0 0 0 0 0 0 0 0 0 1 3 3 4 5 7 7 5 5 4 3 2 1 1];
        
        % CDFs of mobility model
        toCommutingCDF
        fromCommutingCDF
        toBusinessCDF
        fromBusinessCDF
        toEducationCDF
        fromEducationCDF
        toEscortEducationCDF
        fromEscortEducationCDF
        toShoppingCDF
        fromShoppingCDF
        toPersonalBusinessOtherEscortCDF
        fromPersonalBusinessOtherEscortCDF
        toVisitFriendsSportCDF
        fromVisitFriendsSportCDF
        toHolidayCDF
        fromHolidayCDF
        
        % Average trip length in miles Table NTS0405
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
        
        % parameters
        daysPerYear = 365;
        individualDayOfWeekPerYear = 1;
        
        % Round-trip rates as car driver by day of weak
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
    end
    
    methods
        % Constructor
        function obj = evAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.vBase         [V]. Voltage base
            % - param.sMax          [VA]. Is the maximum apparent power of the EV charging
            % - param.pRated        [W]. Is the rated power of the charging spot
            % - param.pRate         [W/s]. Constraint on power rate of change
            % - param.eRated        [J]. Is the rated capacity of the EV battery    
            % - param.Ts            [s]. Sampling time
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.sBase = param.sBase;
            obj.vBase = param.vBase;
            obj.sMax = param.sMax;
            obj.pRated = param.pRated;
            obj.pRate = param.pRate*param.Ts; 
            obj.eRated = obj.eRated;
            obj.Ts = param.Ts;
            obj.onPU = param.onPU;
            
            % Compute parameters
            obj.eOld = obj.eRated;
            obj.individualDayOfWeekPerYear = obj.daysPerYear/7;
            
            % Setup CDFs of mobility model
            [obj.toCommutingCDF, obj.fromCommutingCDF] = getToFromCDFs(obj.toCommute,obj.fromCommute);
            [obj.toBusinessCDF, obj.fromBusinessCDF] = getToFromCDFs(obj.toBusiness,obj.fromBusiness);
            [obj.toEducationCDF, obj.fromEducationCDF] = getToFromCDFs(obj.toEducation,obj.fromEducation);
            [obj.toEscortEducationCDF, obj.fromEscortEducationCDF] = getToFromCDFs(obj.toEscortEducation,obj.fromEscortEducation);
            [obj.toShoppingCDF, obj.fromShoppingCDF] = getToFromCDFs(obj.toShopping,obj.fromShopping);
            [obj.toPersonalBusinessOtherEscortCDF, obj.fromPersonalBusinessOtherEscortCDF] = getToFromCDFs(obj.toOtherPersonal,obj.fromOtherPersonal);
            [obj.toVisitFriendsSportCDF, obj.fromVisitFriendsSportCDF] = getToFromCDFs(obj.toVisit,obj.fromVisit);
            [obj.toHolidayCDF, obj.fromHolidayCDF] = getToFromCDFs(obj.toHoliday,obj.fromHoliday);
            
            % Round-trip rates as car driver (note that probability of car by purpose is only available for groups of purposes)
            obj.CommutingPerYear = obj.CommutingPerYear*obj.CommutingCarDriverProb;
            obj.BusinessPerYear = obj.BusinessPerYear*obj.BusinessCarDriverProb; 
            obj.EducationPerYear = obj.EducationPerYear*obj.EducationOrEduEscortCarDriverProb; 
            obj.EscortEducationPerYear = obj.EscortEducationPerYear*obj.EducationOrEduEscortCarDriverProb; 
            obj.ShoppingPerYear = obj.ShoppingPerYear*obj.ShoppingCarDriverProb; 
            obj.OtherEscortPerYear = obj.OtherEscortPerYear*obj.OtherEscortCarDriverProb; 
            obj.PersonalBusinessPerYear = obj.PersonalBusinessPerYear*obj.PersonalBusinessCarDriverProb; 
            obj.VisitFriendsPrivateHomePerYear = obj.VisitFriendsPrivateHomePerYear*obj.LeisureCarDriverProb; 
            obj.VisitFriendsElsewherePerYear = obj.VisitFriendsElsewherePerYear*obj.LeisureCarDriverProb; 
            obj.SportPerYear = obj.SportPerYear*obj.LeisureCarDriverProb; 
            obj.HolidayPerYear = obj.HolidayPerYear*obj.LeisureCarDriverProb; 
        end
        
        
        % sample electrical vehicle
        function [p,q,e,away] = sample(obj,k,day,v,pRef,qRef,vRef)
            % Input:
            %   - k [-], is the sample number.
            %   - day [day number], is the day of the year 1-365.
            %   - v [-], is the voltage at the connection point.
            %   - pRef [W], is the reference to active power.
            %   - qRef [VAR], is the reference to reactive power.
            %   - vRef [V], is the voltage reference.
            %
            % Output:
            %   - p [W], is the active power output.
            %   - q [VAR], is the reactive power output.
            %   - e [J], is the EV batery energy level.
            %   - away [-]. Indicates if the EV is available for charging.
            
            % If day changes calculate when EV is away from home and how
            % much the state of charge has changed.
            if day ~= obj.dayOld
                [obj.away, obj.dSOC] = availability(obj,day);
                obj.dayOld = day;
                obj.h = 2;
            end
                
            % Active power control
            [p,e] = pCtrl(obj,pRef);         
            
            % Reactive power control
            q = qCtrl(obj,p,v,qRef,vRef);
            
            % Normalize Output to Per Unit (PU)
            if obj.onPU == true
                p = p/obj.sBase;
                q = q/obj.sBase;
                e = e/obj.eRated;
            end
            
            % Iterate hour of the day.
            if ~mod(k,ceil(60*60/obj.Ts))
                obj.h = obj.h + 1;
                obj.temp = 0;
            end
            p = -p;                 % Consumption is negative power flow
            q = -q;
            away = obj.away(obj.h); % Flag can be used by external controller

        end
        
        % Active power control
        function [p,e] = pCtrl(obj,pRef)
            if obj.away(obj.h) == 1 && obj.away(obj.h-1) == 0
                if obj.temp == 0
                    obj.eOld = max(0,obj.eOld+obj.eOld*sum(obj.dSOC)/100);
                end
                p = 0;
                e = 0;
                obj.temp = 1;
            elseif obj.away(obj.h) == 1
                p = 0;
                e = 0;
            elseif obj.away(obj.h) == 0
                if obj.pMode == 0
                    pRef = obj.pRated;
                elseif obj.pMode == 1
                    pRef = pRef;
                else
                    error('Invalid mode for active power control.');
                end
                [p,e] = charging(obj,pRef);
                
            else
                p = 0;
                e = 0;
            end
        end
        
        % Battery charging
        function [p,e] = charging(obj,pRef)
            
            % Handle constraint on charge rate
            if abs(pRef - obj.pOld) > obj.pRate
                pRef = obj.pOld + obj.pRate*sign(pRef - obj.pOld);
            end
            
            if obj.eOld >= obj.eRated && pRef > 0
                e = obj.eRated;
                p = 0;
            elseif obj.eOld <= 0 && pRef < 0
                e = 0;
                p = 0;
            else
                if pRef <= 0
                    p = 0;
                elseif pRef >= obj.pRated
                    p = obj.pRated;
                else
                    p = pRef;
                end
                e = obj.a*obj.eOld + obj.Ts*obj.eta*p;
            end
            obj.eOld = e;
            obj.pOld = p;
        end
        
        % Reactive power control
        function q = qCtrl(obj,p,v,qRef,vRef)
            % Calculate Reactive Power Reference
            if(obj.qMode == 0)
                qRef = p*tan(acos(obj.PF));
            elseif(obj.qMode == 1)
                qRef = qRef;
            elseif(obj.qMode == 2)
                qRef = droop(obj,v,vRef,'');
            else
                error('Invalid mode for reactive power control.');
            end
            % Calculate Reactive Power Output
            qMax = sqrt(obj.sMax^2-p^2); % [Var] Upper limit on reactive power
            if(abs(qRef)>qMax)
                q = sign(qRef)*qMax;
            else
                q = qRef;
            end
        end
        
        % Voltage reactive power droop control function
        function qRef = droop(obj,v,vRef,qFun)
            % Choose between functions
            switch qFun
                % case 'qFun'
                % Input user defined droop control functions and add to
                % class
                otherwise
                    % Choose the default droop control function
                    qRef = defaultDroop(obj,v,vRef);
            end
        
        end
        
        % Default droop control function
        function y = defaultDroop(obj,v,vRef)
            % defaultDroop is a active power droop control.
            %
            %   Input:
            %       - v is the bus voltage [PU or V].
            %       - vRef is the voltage reference [PU or V].
            %
            %   Output:
            %       - pRef is the active power reference [VAR].
            
            if obj.onPU == true
                vMax = 1.1;
                vMin = 0.9;
                k = -2/(vMax-vMin);
                y = obj.sMax*k*(vRef-v)/obj.vBase;
            else
                vMax = obj.vBase * 1.1;
                vMin = obj.vBase * 0.9;
                k = -2/(vMax-vMin);
                y = obj.sMax*k*(vRef-v);
            end 
        end
        
        
        % Calculate if car is available at home or away and the change its
        % state of charge
        function [away, dSOC]  = availability(obj,dayNumber)
            if ~mod(dayNumber,7)
                dayNumber = 7;
            else
                dayNumber = mod(dayNumber,7);
            end
            
            dSOC = zeros(1,length(obj.fromBusiness));
            CommutingAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.CommutingPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [CommutingAway,CommutingStart] = getAwayPeriod(obj.toCommutingCDF,obj.fromCommutingCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(CommutingStart) = dSOC(CommutingStart)-obj.CommutingLength*obj.avgMileage;
            end
            BusinessAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.BusinessPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [BusinessAway,BusinessStart] = getAwayPeriod(obj.toBusinessCDF,obj.fromBusinessCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(BusinessStart) = dSOC(BusinessStart)-obj.BusinessLength*obj.avgMileage;
            end
            EducationAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.EducationPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [EducationAway,EducationStart] = getAwayPeriod(obj.toEducationCDF,obj.fromEducationCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(EducationStart) = dSOC(EducationStart)-obj.EducationLength*obj.avgMileage;
            end
            EscortEducationAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.EscortEducationPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [EscortEducationAway,EscortEducationStart] = getAwayPeriod(obj.toEscortEducationCDF,obj.fromEscortEducationCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(EscortEducationStart) = dSOC(EscortEducationStart)-obj.EscortEducationLength*obj.avgMileage;
            end
            ShoppingAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.ShoppingPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [ShoppingAway,ShoppingStart] = getAwayPeriod(obj.toShoppingCDF,obj.fromShoppingCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(ShoppingStart) = dSOC(ShoppingStart)-obj.ShoppingLength*obj.avgMileage;
            end
            OtherEscortAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.OtherEscortPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [OtherEscortAway,OtherEscortStart] = getAwayPeriod(obj.toPersonalBusinessOtherEscortCDF,obj.fromPersonalBusinessOtherEscortCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(OtherEscortStart) = dSOC(OtherEscortStart)-obj.OtherEscortLength*obj.avgMileage;
            end
            PersonalBusinessAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.PersonalBusinessPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [PersonalBusinessAway,PersonalBusinessStart] = getAwayPeriod(obj.toPersonalBusinessOtherEscortCDF,obj.fromPersonalBusinessOtherEscortCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(PersonalBusinessStart) = dSOC(PersonalBusinessStart)-obj.PersonalBusinessLength*obj.avgMileage;
            end
            VisitFriendsPrivateHomeAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.VisitFriendsPrivateHomePerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [VisitFriendsPrivateHomeAway,VisitFriendsPrivateHomeStart] = getAwayPeriod(obj.toVisitFriendsSportCDF,obj.fromVisitFriendsSportCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(VisitFriendsPrivateHomeStart) = dSOC(VisitFriendsPrivateHomeStart)-obj.VisitFriendsPrivateHomeLength*obj.avgMileage;
            end
            VisitFriendsElsewhereAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.VisitFriendsElsewherePerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [VisitFriendsElsewhereAway,VisitFriendsElsewhereStart] = getAwayPeriod(obj.toVisitFriendsSportCDF,obj.fromVisitFriendsSportCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(VisitFriendsElsewhereStart) = dSOC(VisitFriendsElsewhereStart)-obj.VisitFriendsElsewhereLength*obj.avgMileage;
            end
            SportAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.SportPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [SportAway,SportStart] = getAwayPeriod(obj.toVisitFriendsSportCDF,obj.fromVisitFriendsSportCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(SportStart) = dSOC(SportStart)-obj.SportLength*obj.avgMileage;
            end
            HolidayAway = zeros(1,length(obj.fromBusiness));
            if(rand<obj.HolidayPerYear(dayNumber)/obj.individualDayOfWeekPerYear)
                [HolidayAway,HolidayStart] = getAwayPeriod(obj.toHolidayCDF,obj.fromHolidayCDF);
                % Calculate change in state of charge (SoC), where negative means discharge
                dSOC(HolidayStart) = dSOC(HolidayStart)-obj.HolidayLength*obj.avgMileage;
            end
            away = CommutingAway|BusinessAway|EducationAway|EscortEducationAway|ShoppingAway|...
                    OtherEscortAway|PersonalBusinessAway|VisitFriendsPrivateHomeAway|...
                        VisitFriendsElsewhereAway|SportAway|HolidayAway;
            away = [away 0 0]; 
        end
        
        
        %% Set functions
        % Set power factor
        function setPF(obj,PF)
            if PF<0
                PF = 0;
            elseif PF>1;
                PF = 1;
            end
            obj.PF = PF;
        end
        
        % Set pMode 
        function setPmode(obj,pMode)
            % Input:
            %   - pMode=0: Charge as fast as possible, i.e., pRef = pRated.
            %   - pMode=1: Follow power reference
            
            % Set mode
            obj.pMode = pMode;
        end
        
        % Set qMode 
        function setQmode(obj,qMode)
            % Input:
            %   - pMode=0: Follow active power reference.
            %   - pMode=1: Active power is controlled as p = pFun(v,vRef) [W]
            
            % Set mode
            obj.qMode = qMode;
        end

    end
    
end

%% Internal functions
% Mobility model functions
% Away period
function [away,awayStart] = getAwayPeriod(toCDF,fromCDF)
    away = zeros(1,length(toCDF));
    % Pick from distribution
    eventTime1 = find(rand<toCDF,1,'first');
    eventTime2 = find(rand<fromCDF,1,'first');
    awayStart = min(eventTime1,eventTime2);
    % Away hours
    for i = awayStart:1:max(eventTime1,eventTime2)
        away(i) = true;
    end
end
        
% To/from CDF's
function [toCDF,fromCDF] = getToFromCDFs(toPDF,fromPDF)
    toCDF = zeros(1,length(toPDF));
    fromCDF = zeros(1,length(fromPDF));
    toNormalize = sum(toPDF);
    fromNormalize = sum(fromPDF);
    toPDF = toPDF/toNormalize;
    fromPDF = fromPDF/fromNormalize;

    toCDF(1) = toPDF(1);
    fromCDF(1) = fromPDF(1);
    for i=2:1:length(toPDF)
        toCDF(i) = toCDF(i-1)+toPDF(i);
        fromCDF(i) = fromCDF(i-1)+fromPDF(i);
    end
end

