classdef scAsset < handle
    %scAsset implements a simple controllable static compensator asset.
    %   The scAsset can be set to inject or consume reactive power based
    %   on a reference  or do voltage droop control.
    %   
    %   This model is based on a simple algebraic principle, as it is
    %   intended for sampling times higher than 1 second.
    %   
    % Revision:
    % 23-06-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 10e3;               % [VA]. Complex power base
        vBase = 400;                % [V]. Voltage base
        qRatedMax = 10e3;           % [VAR]. Rated reactive power output
        qRatedMin = -10e3;          % [VAR]. Rated reactive power consumption
        Qmode = 0;                  % [-]. Reactive power control mode
        onPU = false;               % [-]. Indicates if the static compensator is used to simulate a system on per unit basis. true = yes, false = no
    end
    
    methods
        % Constructor
        function obj = scAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.vBase         [V]. Voltage base
            % - param.qRatedMax     [VAR]. Maximum reactive power
            % - param.qRatedMin     [VAR]. Minimum reactive power    
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.sBase = param.sBase;
            obj.vBase = param.vBase;
            obj.qRatedMax = param.qRatedMax;
            obj.qRatedMin = param.qRatedMin;
            obj.onPU = param.onPU;
        end
        
        
        % sample capacitor bank system
        function q = sample(obj,v,qRef,vRef)
            % Input:
            %   - v [V], is the voltage magnitude on the bus.
            %   - qRef [VAR], is the reference to reactive power.
            %   - vRef [V], is the voltage reference.
            %
            % Output:
            %   - q [VAR], reactive power output.
            
            if obj.onPU == true
                % Convert References to Physical Quantities
                qRef = qRef*obj.sBase;
            end
            
            % Reactive power control
            q = qCtrl(obj,v,qRef,vRef);
            
            if obj.onPU == true
                % Normalize Output to Per Unit (PU)
                q = q/obj.sBase;
            end
        end
        
        % Reactive power control
        function q = qCtrl(obj,v,qRef,vRef)
            % Calculate Reactive Power Reference
            if(obj.Qmode == 0)
                qRef = qRef;
            elseif(obj.Qmode == 1)
                qRef = droop(obj,v,vRef,'');
            else
                error('Invalid mode for reactive power control.');
            end
            % Calculate Reactive Power Output and limit
            if qRef>obj.qRatedMax
                q = obj.qRatedMax;
            elseif qRef < obj.qRatedMin;
                q = obj.qRatedMin;
            else
                q = qRef;
            end
        end
        
        % Voltage droop control function
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
        function qRef = defaultDroop(obj,v,vRef)
            % defaultQDroop is a reactive power droop control.
            %
            %   Input:
            %       - v is the bus voltage [PU].
            %       - vRef is the voltage reference [PU].
            %
            %   Output:
            %       - qRef is the reactive power reference [VAR].
            vMax = 1.1;
            vMin = 0.9;

            k = -2/(vMax-vMin);
            qRef = max(obj.qRatedMax,abs(obj.qRatedMin))*k*(v-vRef)/obj.vBase;
        end
        
        %% Set functions
        % Set Qmode and Q control function
        function setQmode(obj,Qmode)
            % Input:
            %   - Qmode=0: Reactive power is to qRef [PU].
            %   - Qmode=1: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
            
            % Set mode
            obj.Qmode = Qmode;
        end

    end
    
end

