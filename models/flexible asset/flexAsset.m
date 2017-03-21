classdef flexAsset < handle
    %flexAsset implements a generic flexible asset.
    %   The flexAsset can be set to consume active and/or reactive power
    %   or inject active/reactive power.   
    %   
    % Revision:
    % 05-03-2015, R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 10e3;               % [VA]. Complex power base
        sMax = 10e3;                % [VA]. Maximum apparent power
        vBase = 400;                % [V]. Voltage base
        pRate = 20;                 % [W/sample]. Constraint on power rate of change (NOTE: not implemented)
        pOld = 0;                   % [W]. Placeholder for power. Used for rate constraint
        onPU = false;               % [-]. Indicates if the energy storage is used to simulate a system on per unit basis. true = yes, false = no
        pRefOld = 0;                % [W]. Used to handle communication loss.            
        qRefOld = 0;                % [Var]. Used to handle communication loss.
    end
    
    methods
        % Constructor
        function obj = flexAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.sMax          [VA]. Is the maximum apparent power of the energy storage
            % - param.vBase         [V]. Voltage base
            % - param.pRate         [W]. Constraint om power rate of change   
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.vBase = param.vBase;
            obj.sMax = param.sMax;
            obj.onPU = param.onPU;
            
            
            if obj.onPU == true
                obj.sBase = param.sBase;
            end
        end
        
        
        % sample flexible asset
        function [p,q,satUp,satDown,pMin,pMax,qMin,qMax] = sample(obj,pRef,qRef)
            % Input:
            %   - pRef, is the reference to active power.
            %   - qRef, is the reference for reactive power.
            %
            % Output:
            %   - p [W], is the active power consumption.
            %   - q [VAR], is the reactive power output.
            %   - satUp [-], indicates that the asset is in upper saturation.
            %   - satDown [-], indicates that the asset is in lower saturation.
            %   - pFlexUp [W], is upwards flexibility in active power.
            %   - pFlexDown [W], is downwards flexibility in active power.
            %   - qFlexUp [Var], is upwards flexibility in reactive power.
            %   - qFlexDown [Var], is downwards flexibility in reactive power.
            
            % Check if system is in saturation
            % Active power
            satUp = 0;
            satDown = 0;
            if isnan(pRef)
                pRef = obj.pRefOld;
            end
            if pRef >= obj.sMax
                pSet = obj.sMax;
                satUp = 1;
            elseif pRef <= -obj.sMax
                pSet = -obj.sMax;
                satDown = 1;
            else
                pSet = pRef;
            end
            obj.pRefOld = pSet;
            % Reactive power
            if isnan(qRef)
                qRef = obj.qRefOld;
            end
            qMax = sqrt(obj.sMax^2-pSet^2);
            if qRef > qMax
                qSet = qMax;
            elseif qRef < -qMax
                qSet = -qMax;
            else
                qSet = qRef;
            end
            obj.qRefOld = qSet;
            
            % Calculate flexibility
            pMax = obj.sMax;
            pMin = -obj.sMax;
            qMax = sqrt(obj.sMax^2-pSet^2);
            qMin = sqrt(obj.sMax^2-pSet^2);
            
            % Set output
            p = pSet;
            q = qSet;
            
            if obj.onPU == true
                % Normalize Output to Per Unit (PU)
                p = p/obj.sBase;
                q = q/obj.sBase;
            end
        end
        


    end
    
end

