classdef cbAsset < handle
    %cbAsset implements a capacitor bank.
    %   The cbAsset can be controlled to inject reactive power into the
    %   grid in steps. It is characterized by having a number of steps and
    %   a reactive output per step (i.e., it can inject from 0 - Qmax in 
    %   N steps). 
    %   
    %   
    % Revision:
    % 15-10-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 10e3;               % [VA]. Complex power base
        vBase = 400;                % [V]. Voltage base
        qMax = 10e3;                % [VAR]. Maximum reactive power output
        Ts = 1;                     % [sec]. Sampling time
        onPU = false;               % [-]. Indicates if the energy storage is used to simulate a system on per unit basis. true = yes, false = no
        qMode = 0;                  % [-]. Mode for control (0 = follow reference, 1 = voltage control)
        
        qStep                       % [VAR], Reactive power per step.
        nStepPrDay = 20;            % [-]. Maximum number of steps per day
        tMin = 10*60;               % [s]. Minimum time between steps
        kMin                        % [-]. Minimum number of samples between steps
        stepPos = 0;                % [-]. Position of the capacitor bank (how many capacitors are coupled in)
        nSteps = 10;                % [-]. Number of steps to reach qMax.
        stepOk = 1;                 % [-]. Flag indicating if step is allowed (1=yes, 0=no) 
        stepChDay = 0;              % [-]. Number of steps that has occured in a day
        dayOld = -1;                % [-]. Used to check if day has changed
        stepTimer = 0;              % [-]. Timer indicating if a step is allowed, i.e., the timer must be greater than kMin to apply change
        vHyst = 0.05;               % [%]. Hysterisis limit for voltage control
        
        avgFilterBuf                % [-]. Buffer for the moving average filter
        avgFilterTime = 10*60;      % [sec]. Time window of the moving average filter
    end
    
    methods
        % Constructor
        function obj = cbAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.vBase         [V]. Voltage base
            % - param.qMax          [VAR]. Maximum reactive power output
            % - param.nSteps        [-]. Number of steps in the capacitor bank
            % - param.Ts            [s]. Sampling time
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)

            % Set parameters
            obj.sBase = param.sBase;
            obj.vBase = param.vBase;
            obj.qMax = param.qMax;
            obj.nSteps = param.nSteps;
            obj.Ts = param.Ts;
            obj.onPU = param.onPU;
            
            % Check if the simulation is on per unit
            if obj.onPU == true
                vAvgFilterTemp = 1;
            else
                vAvgFilterTemp = param.vBase;
            end
            
            % Compute reactive power per step
            obj.qStep = param.qMax/param.nSteps;
            
            % Compute minimum samples between capacitor bank step
            obj.kMin = ceil(obj.tMin/obj.Ts);
            
            % Setup the buffer for the moving average filter
            obj.avgFilterBuf = vAvgFilterTemp*ones(ceil(obj.avgFilterTime/obj.Ts),1);
        end
        
        
        % sample capacitor bank
        function q = sample(obj,v,vRef,uRef,k,day)
            % Input:
            %   - v [V], is the voltage at the connection point.
            %   - qRef [VAR], is the reference to active power.
            %   - vRef [V], is the voltage reference for the droop control.
            %
            % Output:
            %   - q [VAR], is the reactive power output    
            
            % If day changes reset the number of steps a day
            if day ~= obj.dayOld
                obj.stepChDay = 0;
                obj.dayOld = day;
            end
            
            % Reactive power control
            q = qCtrl(obj,v,vRef,uRef,k);
            
            % Check if simulation is on PU basis
            if obj.onPU == true
                % Normalize Output to Per Unit (PU)
                q = q/obj.sBase;
            end

        end
        
        % Reactive power control
        function q = qCtrl(obj,v,vRef,uRef,k)
            % Compute the tap reference
            if obj.qMode == 0
                uRef = uRef;
            elseif obj.qMode == 1
                uRef = voltageControl(obj,v,vRef);
            else
                error('Invalid mode for tap changing transformer.');
            end
            
            % Check if step changing is allowed
            stepChAllowed(obj,k);
            
            % Check if step has changed
            if uRef ~= obj.stepPos
                q = stepChangingCtrl(obj,uRef,k);
            else
                q = obj.stepPos*obj.qStep;
            end
            
        end
        
        % Step changing control
        function q = stepChangingCtrl(obj,uRef,k)
            if uRef>obj.stepPos && obj.stepPos<obj.nSteps && obj.stepOk == 1
                obj.stepPos = obj.stepPos + 1;
                q = obj.stepPos*obj.qStep;
                obj.stepChDay = obj.stepChDay + 1;
                obj.stepTimer = k+obj.kMin;
            elseif uRef<obj.stepPos && obj.stepPos>0 && obj.stepOk == 1
                obj.stepPos = obj.stepPos - 1;
                q = obj.stepPos*obj.qStep;
                obj.stepChDay = obj.stepChDay + 1; 
                obj.stepTimer = k+obj.kMin;
            else
                obj.stepPos = obj.stepPos;
                q = obj.stepPos*obj.qStep;
            end
        end
        
        % Voltage control
        function uRef = voltageControl(obj,v,vRef)
            % Moving average filter on voltage
            v = avgFilter(obj,v);
            % Voltage control algorithm
            if v > vRef + obj.vHyst && obj.stepPos > 0 
                uRef = obj.stepPos - 1;
            elseif v < vRef - obj.vHyst && obj.stepPos < obj.nSteps
                uRef = obj.stepPos + 1; 
            else
                uRef = obj.stepPos;
            end
        end
        
        % Check if step change is allowed
        function stepChAllowed(obj,k)
            if obj.stepChDay >= obj.nStepPrDay
                obj.stepOk = 0;
            elseif k < obj.stepTimer
                obj.stepOk = 0;
            elseif k > obj.stepTimer
                obj.stepOk = 1;
            end
        
        end
        
        % Moving average filter on voltage for the stepping
        function vOut = avgFilter(obj,v)
            % Update data
            obj.avgFilterBuf(2:end) = obj.avgFilterBuf(1:end-1);
            obj.avgFilterBuf(1) = v;
            vOut = sum(obj.avgFilterBuf)/length(obj.avgFilterBuf);
        end
        
        %% Set functions  
        % Set pMode 
        function setQmode(obj,qMode)
            % Input:
            %   - qMode=0: Follow active power reference.
            %   - qMode=1: Active power is controlled as p = pFun(v,vRef) [W]
            
            % Set mode
            obj.qMode = qMode;
        end


    end
    
end

