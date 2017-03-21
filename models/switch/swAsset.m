classdef swAsset < handle
    %swAsset implements a controllable switch.
    %   The swAsset can change the switch position between open and closed.
    %   The switch is used to change the electrical grid topology.
    %   To ensure realistic simulation scenarios a minimum time between
    %   switch state change has been incorporated, along with a maximum
    %   number of state changes per day.
    %   
    %   The model is intended for studing voltage control problems and
    %   power loss minimization. The switch works by inserting an
    %   admittance into the admittance matrix whenever the switch state is
    %   closed.
    %   
    % R. Pedersen 05-9-2014, Aalborg University
    
    properties
        zBase = 1;                  % [Ohm]. Base impedance of the grid
        Z = 0.03-1i*0.13;           % [Ohm]. Switch impedance
        Y = 1;                      % [S]. Admittance of switch
        Ts = 1;                     % [s]. Sampling time
        onPU = false;               % [-]. Indicates if the switch is used to simulate a system on per unit basis. true = yes, false = no
        nSwPrDay = 10;              % [-]. Maximum number of switch state changes per day
        tMin = 10*60;               % [s]. Minimum time between switch position changes
        kMin                        % [-]. Minimum number of samples between switch position changes
        swChOk = 1;                 % [-]. Flag indicating if switch state can be changed (1=yes, 0=no) 
        swState = 0;                % [-]. Switch state (0=open, 1=closed)
        swChDay = 0;                % [-]. Number of switch state changes that has occured in a day
        dayOld = -1;                % [-]. Used to check if day has changed
        swTimer = 0;                % [-]. Timer indicating if switch state change is allowed, i.e., the timer must be greater than kMin to apply change
    end
    
    methods
        % Constructor
        function obj = swAsset(param)
            % - param format:
            % - param.zBase           [Ohm]. Base impedance used for per unit simulation
            % - param.Ts              [s]. Sampling time
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.Ts = param.Ts;
            obj.onPU = param.onPU;
            % Check if simulation is on PU
            if obj.onPU == true
                obj.zBase = param.zBase;
            	obj.Z = param.Z/param.zBase;
            else
                obj.Z = param.Z;
            end
            obj.Y = 1/obj.Z;
            
            % Compute minimum samples until switch state change
            obj.kMin = ceil(obj.tMin/obj.Ts);    
        end
        
        
        % sample switch
        function [state,y] = sample(obj,uRef,k,day)
            % Input:
            %   - uRef [-], is the reference to switch state (0=open, 1=closed).
            %   - k [-], is the sample.
            %   - day [-], is the day of the simulation (1-365).
            %
            % Output:
            %   - state [-], is the state of the switch, either open or closed.
            %   - y [-], is the admittance to be put into the grid admittance matrix.
            
            % If day changes reset the number of switches a day
            if day ~= obj.dayOld
                obj.swChDay = 0;
                obj.dayOld = day;
            end
            
            % Switch state control
            state = swCtrl(obj,uRef,k);
            y = obj.Y;
        end
        
        % Switch state control
        function state = swCtrl(obj,uRef,k)
            
            % Check if switch state changing is allowed
            swChAllowed(obj,k);
            
            % Check if switch state has changed
            if uRef ~= obj.swState
                state = swChangingCtrl(obj,uRef,k);
            else
                state = obj.swState;
            end         
        end
        
        % Switch state changing control
        function state = swChangingCtrl(obj,uRef,k)
            if obj.swChOk == 1
                state = uRef;
                obj.swState = state;
                obj.swChDay = obj.swChDay + 1;
                obj.swTimer = k+obj.kMin;
            else
                state = obj.swState;
            end
        end
        
        % Check if switch state change is allowed
        function swChAllowed(obj,k)
            if obj.swChDay >= obj.nSwPrDay
                obj.swChOk = 0;
            elseif k < obj.swTimer
                obj.swChOk = 0;
            elseif k > obj.swTimer
                obj.swChOk = 1;
            end
        
        end
        
        %% Set functions
        % Set tap specifications
        function setSwSpec(obj,nrSwPrDay,tMin)
            % Input:
            %   - nrSwPrDay [-], is the number of switches per day
            %   - tMin [s], is the minimum time between switches
            obj.nSwPrDay = nrSwPrDay;
            obj.tMin = tMin;
            obj.kMin = ceil(tMin/obj.Ts);
        end


    end
    
end

