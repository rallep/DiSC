classdef esAsset < handle
    %esAsset implements a simple controllable energy storage.
    %   The esAsset can be set to consume active power or inject active
    %   power. The energy storage consumes reactive power with a constant
    %   power factor. The energy storage can also be set to do voltage
    %   droop control.
    %   
    %   
    % Revision:
    % 23-06-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 10e3;               % [VA]. Complex power base
        vBase = 400;                % [V]. Voltage base
        pRatedMax = 10e3;           % [W]. Rated power consumption
        pRatedMin = -10e3           % [W]. Rated power output
        pRate = 20;                 % [W/sample]. Constraint on power rate of change
        pOld = 0;                   % [W]. Placeholder for power. Used for rate constraint
        sMax = 10e3;                % [VA]. Maximum apparent power
        PF = 1;                     % [-]. Power factor
        eRated = 10e6;              % [J]. Rated energy leve
        eOld = 0;                   % [J]. Placeholder for energy level
        Ts = 1;                     % [sec]. Sampling time
        onPU = false;               % [-]. Indicates if the energy storage is used to simulate a system on per unit basis. true = yes, false = no
        a = 1;                      % [-]. Energy storage drain rate
        etaIn = 1;                  % [-]. Energy storage effeciency when charging
        etaOut = 1;                 % [-]. Energy storage efficiency when discharging
        pMode = 0;                  % [-]. Mode for the control of active power
        qMode = 0;                  % [-]. Mode for the control of reactive power
        
        avgFilterBuf                % [-]. Buffer for the moving average filter
        avgFilterTime = 10*60;      % [s]. Time window of the moving average filter
    end
    
    methods
        % Constructor
        function obj = esAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.vBase         [V]. Voltage base
            % - param.sMax          [VA]. Is the maximum apparent power of the energy storage
            % - param.pRatedMax     [W]. Is the maximum power of the energy storage (Consumption)
            % - param.pRatedMin     [W]. Is the minimum power of the energy storage (Injection)  
            % - param.pRate         [W]. Constraint om power rate of change
            % - param.eRated        [J]. Is the rated capacity of the energy storage    
            % - param.Ts            [s]. Sampling time
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.sBase = param.sBase;
            obj.vBase = param.vBase;
            obj.sMax = param.sMax;
            obj.pRatedMax = param.pRatedMax;
            obj.pRatedMin = param.pRatedMin;
            obj.pRate = param.pRate*param.Ts; 
            obj.eRated = param.eRated;
            obj.Ts = param.Ts;
            obj.onPU = param.onPU;
            
            if obj.onPU == true
                vAvgFilterTemp = 1;
            else
                vAvgFilterTemp = param.vBase;
            end
            
            % Setup the buffer for the moving average filter
            obj.avgFilterBuf = vAvgFilterTemp*ones(ceil(obj.avgFilterTime/obj.Ts),1);
        end
        
        
        % sample energy storage system
        function [p,q,e] = sample(obj,v,pRef,qRef,vRef)
            % Input:
            %   - v, is the voltage at the connection point.
            %   - pRef, is the reference to active power.
            %   - vRef, is the voltage reference for the droop control
            %
            % Output:
            %   - p [W], is the active power consumption
            %   - q [VAR], is the reactive power output
            %   - e [J], is the energy storage energy level
            
            % Active power control and energy level
            [p,e] = pCtrl(obj,v,vRef,pRef);
            
            % Reactive power control
            q = qCtrl(obj,p,v,qRef,vRef);
            
            if obj.onPU == true
                % Normalize Output to Per Unit (PU)
                p = p/obj.sBase;
                q = q/obj.sBase;
                e = e/obj.eRated;
            end
            p = p;   % Negative is consumption
            q = q;

        end
        
        % Active power control
        function [p,e] = pCtrl(obj,v,vRef,pRef)
            % Calculate Active Power Output and energy level
            if obj.pMode == 0
            elseif obj.pMode == 1
                %v = avgFilter(obj,v);
                pRef = droopP(obj,v,vRef,'A');
            else
                error('Invalid mode for active power control.');
            end
            
            % Handle constraint on charge rate
            if abs(pRef - obj.pOld) > obj.pRate
                pRef = obj.pOld + obj.pRate*sign(pRef - obj.pOld);
            end
            
            % Handle constraints on maximum and minimum levels
            if obj.eOld >= obj.eRated && pRef <= 0
                e = obj.eRated;
                p = 0;
            elseif obj.eOld <= 0 && pRef >= 0
                e = 0;
                p = 0;
            else
                if pRef <= obj.pRatedMin
                    p = obj.pRatedMin;
                elseif pRef >= obj.pRatedMax
                    p = obj.pRatedMax;
                else
                    p = pRef;
                end
                if p <= 0
                    e = obj.a*obj.eOld + obj.Ts*obj.etaIn*(-p);
                else
                    e = obj.a*obj.eOld + obj.Ts*obj.etaOut*(-p);
                end
            end
            obj.pOld = p;
            obj.eOld = e;
        end
        
        % Reactive power control
        function q = qCtrl(obj,p,v,qRef,vRef)
            % Calculate Reactive Power Reference
            if(obj.qMode == 0)
                qRef = p*tan(acos(obj.PF));
            elseif(obj.qMode == 1)
                qRef = qRef;
            elseif(obj.qMode == 2)
                %v = avgFilter(obj,v);
                qRef = droopQ(obj,v,vRef,'');
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
        
        % Voltage active power droop control function
        function pRef = droopP(obj,v,vRef,pFun)
            % Choose between functions
            switch pFun
                case 'A'
                    pRef = agressiveDroop(obj,v,vRef);
                % Input user defined droop control functions and add to
                % class
                otherwise
                    % Choose the default droop control function
                    pRef = defaultDroop(obj,v,vRef);
            end
        
        end
        
        % Voltage reactive power droop control function
        function qRef = droopQ(obj,v,vRef,qFun)
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
        
        % Agressive droop control function
        function pRef = agressiveDroop(obj,v,vRef)
            k = -0.2;
            if obj.onPU == true
                pRef = obj.sMax*k*(vRef-v)/obj.vBase;
            else
                pRef = obj.sMax*k*(vRef-v);
            end
                
        end
        
        % Moving average filter on voltage active power droop control
        function vOut = avgFilter(obj,v)
            % Update data
            obj.avgFilterBuf(2:end) = obj.avgFilterBuf(1:end-1);
            obj.avgFilterBuf(1) = v;
            vOut = sum(obj.avgFilterBuf)/length(obj.avgFilterBuf);
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
        
        % Set drain rate (For ideal storage, a = 1)
        function setDrain(obj,a)
            if a<0
                a = 0;
            elseif a>1;
                a = 1;
            end
            obj.a = a;
        end
        
        % Set pMode 
        function setPmode(obj,pMode)
            % Input:
            %   - pMode=0: Follow active power reference.
            %   - pMode=1: Active power is controlled as p = pFun(v,vRef) [W]
            
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

