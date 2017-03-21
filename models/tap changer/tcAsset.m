classdef tcAsset < handle
    %tcAsset implements a controllable tap changing transformer.
    %   The tcAsset can change the tap position of a transformer and 
    %   thereby effectively change the transformers ratio.
    %   The tap changer can be used to simulate e.g., an on-load tap changing
    %   transformer. The tap changer works by altering the grid admittance
    %   matrix, when ever the tap position is changed. To ensure realistic
    %   simulation scenarios a minimum time between tap changes has been
    %   incorporated, along with a specified maximum number of tap changes
    %   per day.
    %   
    %   The model is intended for studing voltage control problems and can
    %   be set to run in two different modes; follow tap reference and
    %   apply voltage control. The model is based on the one found in "Power Systems
    %   Stability and Control, P. Kundur".
    %
    % Revision:
    % 03-09-2014, R. Pedersen, Aalborg University. Notes:

    properties
        zBase = 1;                  % [Ohm]. Base impedance on the primary side of the transformer
        vBase = 20e3;               % [V]. Voltage base
        Z = 0.03-1i*0.13;           % [Ohm]. Trafo impedance
        y1 = 1;                     % [S]. Admittance of trafo between busses
        y2 = 1;                     % [S]. Admittance from ground to trafo, primary side
        y3 = 1;                     % [S]. Admittance from ground to trafo, secondary side
        mode = 0;                   % [-]. Mode of the tap position control (0=follow tap reference, 1=voltage control
        Ts = 1;                     % [s]. Sampling time
        onPU = false;               % [-]. Indicates if the transformer is used to simulate a system on per unit basis. true = yes, false = no
        nomTapRatio = 3;            % [-]. Nominal transformer ratio
        nTapUp = 5;                 % [-]. Number of tap position upwards
        nTapDown = -5;              % [-]. Number of tap positions downwards
        nTapPrDay = 10;             % [-]. Maximum number of tap changes per day
        chPrTap = 0.0125;           % [-]. Ratio of change from nominal voltage per tap position on a per unit basis 
        tMin = 5*60;                % [s]. Minimum time between tap changes
        kMin                        % [-]. Minimum number of samples between tap changes
        tapChOk = 0;                % [-]. Flag indicating if tap change is allowed (1=yes, 0=no) 
        tapPos = 0;                 % [-]. Tap position
        tapChDay = 0;               % [-]. Number of tap changes that has occured in a day
        dayOld = -1;                % [-]. Used to check if day has changed
        tapTimer = 0;               % [-]. Timer indicating if tap changing is allowed, i.e., the timer must be greater than kMin to do tap changing
        vHyst = 0.02;               % [PU]. Hysterisis limit for voltage control
        avgFilterBuf                % [-]. Buffer for the moving average filter
        avgFilterTime = 10*60;      % [s]. Time window of the moving average filter
        
        % Placement parameters
        busFrom = nan;                % [-]. Primary side bus connection
        busTo = nan;                  % [-]. Secondary side bus connection
    end
    
    methods
        % Constructor
        function obj = tcAsset(param)
            % Param format
            % - param.vBase         [V]. Voltage base 
            % - param.Ts            [s]. Sampling time
            % - param.nomTapRatio   [-]. Nominal tap ratio e.g., 60/20 will give a nominal tap ratio of 3
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            % - param.Z             [Ohm]. Transformer impedance.
            
            % Set parameters
            obj.Ts = param.Ts;
            obj.vBase = param.vBase;
            obj.nomTapRatio = param.nomTapRatio;
            obj.onPU = param.onPU;
            if obj.onPU == true
                obj.zBase = param.zBase;
                obj.Z = param.Z/param.zBase;
                obj.vHyst = 0.03;
                vAvgFilterTemp = 1;
            else
                obj.Z = param.Z;
                obj.vHyst = 0.05*param.vBase;
                vAvgFilterTemp = param.vBase;
            end
            
            % Compute minimum samples between tap changes
            obj.kMin = ceil(obj.tMin/obj.Ts);
            
            % Setup the buffer for the moving average filter
            obj.avgFilterBuf = vAvgFilterTemp*ones(ceil(obj.avgFilterTime/obj.Ts),1);
            
            % Change the admittances of the trafo
            changeTrafoAdm(obj);      
        end
        
        
        % sample tap changing transformer
        function [y1,y2,y3] = sample(obj,v,vRef,uRef,k,day)
            % Input:
            %   - v, is the voltage at the connection point.
            %   - vRef, is the voltage reference for the voltage
            %     control.
            %   - uRef, is the reference to tap position.
            %   - k, is the sample.
            %   - day, is the day of the simulation (1-365).
            % 
            % Output:
            %   - y1, y2, y3, is the transformer admittances between busses,
            %     from ground to primary side bus and from ground to
            %     secondary side bus, respectively.
          
            
            % If day changes reset the number of tap changes a day
            if day ~= obj.dayOld
                obj.tapChDay = 0;
                obj.dayOld = day;
            end
            
            % Tap position control
            [y1,y2,y3] = tapCtrl(obj,v,vRef,uRef,k);
        end
        
        % Tap position control
        function [y1,y2,y3] = tapCtrl(obj,v,vRef,uRef,k)
            % Compute the tap reference
            if obj.mode == 0
                uRef = uRef;
            elseif obj.mode == 1
                uRef = voltageControl(obj,v,vRef);
            else
                error('Invalid mode for tap changing transformer.');
            end
            
            % Check if tap changing is allowed
            tapChAllowed(obj,k);
            
            % Check if tap has changed
            if uRef ~= obj.tapPos
                [y1,y2,y3] = tapChangingCtrl(obj,uRef,k);
            else
                y1 = obj.y1;
                y2 = obj.y2;
                y3 = obj.y3;
            end         
        end
        
        % Voltage control
        function uRef = voltageControl(obj,v,vRef)
            % Moving average filter on voltage
            v = avgFilter(obj,v);
            if v > vRef + obj.vHyst && obj.tapPos > obj.nTapDown 
                uRef = obj.tapPos - 1;
            elseif v < vRef - obj.vHyst && obj.tapPos < obj.nTapUp
                uRef = obj.tapPos + 1; 
            else
                uRef = obj.tapPos;
            end
        end
        
        % Tap changing control
        function [y1,y2,y3] = tapChangingCtrl(obj,uRef,k)
            if uRef>obj.tapPos && obj.tapPos<obj.nTapUp && obj.tapChOk == 1
                obj.tapPos = obj.tapPos + 1;
                [y1,y2,y3] = changeTrafoAdm(obj);
                obj.tapChDay = obj.tapChDay + 1;
                obj.tapTimer = k+obj.kMin;
            elseif uRef<obj.tapPos && obj.tapPos>obj.nTapDown && obj.tapChOk == 1
                obj.tapPos = obj.tapPos - 1;
                [y1,y2,y3] = changeTrafoAdm(obj);
                obj.tapChDay = obj.tapChDay + 1; 
                obj.tapTimer = k+obj.kMin;
            else
                obj.tapPos = obj.tapPos;
                [y1,y2,y3] = changeTrafoAdm(obj);
            end
        end
        
        % Change admittance matrix
        function [y1,y2,y3] = changeTrafoAdm(obj)
            % Setup trafo admittances
            Ye = 1/obj.Z;
            if obj.onPU == true
                n = 1-obj.chPrTap*obj.tapPos;
            else
                n = obj.nomTapRatio - obj.chPrTap*obj.nomTapRatio*obj.tapPos;
            end
            c=1/n;
            y1 = (c)*Ye;
            y2 = (c^2-c)*Ye;
            y3 = (1-c)*Ye;
            obj.y1 = y1;
            obj.y2 = y2;
            obj.y3 = y3;
        end
        
        % Check if tap change is allowed
        function tapChAllowed(obj,k)
            if obj.tapChDay >= obj.nTapPrDay
                obj.tapChOk = 0;
            elseif k < obj.tapTimer
                obj.tapChOk = 0;
            elseif k > obj.tapTimer
                obj.tapChOk = 1;
            end
        
        end
        
        % Moving average filter on voltage for the tap changing
        function vOut = avgFilter(obj,v)
            % Update data
            obj.avgFilterBuf(2:end) = obj.avgFilterBuf(1:end-1);
            obj.avgFilterBuf(1) = v;
            vOut = sum(obj.avgFilterBuf)/length(obj.avgFilterBuf);
        end
        
        %% Set functions
        % Set mode 
        function setMode(obj,mode)
            % Input:
            %   - pMode=0: Follow tap position reference.
            %   - pMode=1: Tap position is controlled according to voltage
            
            % Set mode
            obj.mode = mode;
        end
        
        % Set tap specifications
        function setTapSpec(obj,nrTapUp,nrTapDown,nrTapPrDay,tMin,chPrTap,hyst)
            % Input:
            %   - nrTapUp, is the number of tap positions upwards (Default 5)
            %   - nrTapDown, is the number of tap positions downwards (Default 5)
            %   - nrTapPrDay, is the maximum number of tap changes per day (Default 10)
            %   - tmin, is the minimum time between taps (Default 5 min, 5*60)
            %   - chPrTap, is the change in ratio from nomial position on a
            %     per unit basis (Default 0.0125, equal to 1.25%)
            %   - hyst, is the hysterisis on the voltage when applying
            %     voltage control (Default is 0.02 p.u.)
            
            obj.nTapUp = nrTapUp;
            obj.nTapDown = -nrTapDown;
            obj.nTapPrDay = nrTapPrDay;
            obj.tMin = tMin;
            obj.kMin = tMin/obj.Ts;
            obj.chPrTap = chPrTap;
            if obj.onPU == true
                obj.vHyst = hyst;
            else
                obj.vHyst = hyst*obj.vBase;
            end
        end


    end
    
end

