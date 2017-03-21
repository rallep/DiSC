classdef smAsset < handle
%Supermarket with flexible defrost cycles for demand response purpose
%   This class implements a model of a supermarket with flexible
%   defrost cycles. The defrost cycles of low temperature display cases
%   are allowed to be shifted in time under certain bounds. The reason
%   for only using low temperature dispaly cases is that these use
%   electrical heating elements for actively melting the ice during
%   defrost cycles. These heating elements consume a substantial amount
%   of power.
%
%   For a formal description of supermarket defrost flexibility, see the
%   note on the homepage:
%   http://kom.aau.dk/project/SmartGridControl/DiSC/documentation.html
%
% Revision:
% 16-07-2014, R. Pedersen, Aalborg University. Notes:


    % Public properties
    properties
        % Setup supermarket 
        sBase = 3e6;            % [VA]. Complex power base
        Ts = 1;                 % [sec]. Sampling time
        onPU = false;           % [-]. Indicates if the supermarket is used to simulate a system on per unit basis. true = yes, false = no
        numDF = 1;              % [-]. Number of low temperature display cases
        pDefrost                % [W]. Vector containing defrost power consumption
        readyTimers             % [sec]. Set of ready state timers, one for each defrost
        defrostTimers           % [sec]. Set of defrost state timers
        waitTimers              % [sec]. Set of wait state timers
        states                  % [-]. Set of states (0 = ready, 1 = defrost, 2 = wait)
        PF = 0.95;              % [-]. Power factor
        readyGuards             % [sec]. Guards for the ready state
        defrostGuards           % [sec]. Guards for the defrost state
        waitGuards              % [sec]. Guards for the wait state
    end

    methods
        % Constructor
        function obj = smAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis   
            % - param.Ts            [s]. Sampling time
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            % - param.pDefrost      [W]. Is a vector containing the power consumption of each defrost cycle
            
            % Set parameters
            obj.sBase = param.sBase;
            obj.Ts = param.Ts;
            obj.onPU = param.onPU;
            obj.pDefrost = param.pDefrost;
            
            % Calculate other parameters 
            obj.numDF = length(param.pDefrost);
            
            % Initialize timers and states
            obj.readyTimers = zeros(obj.numDF,1);
            obj.defrostTimers = zeros(obj.numDF,1);
            obj.waitTimers = ones(obj.numDF,1)*60*60*3;     % First defrost from 7-9
            obj.states = ones(obj.numDF,1)*2;               % Start all defrost in wait state
            obj.readyGuards = ones(obj.numDF,1)*60*60*2;    % Ready state for 2 hours
            obj.defrostGuards = ones(obj.numDF,1)*45*60;    % Defrost state for 45 minutes
            obj.waitGuards = ones(obj.numDF,1)*60*60*10;    % Wait state for 10 hours             
            
        end

        % Run supermarket
        function [p,q,timeReady,timeDefrost] = sample(obj,exeDefrost,pRS)
            % Input:
            %   - exeDefrost [-], is the external control signal to execute defrost cycles.
            %   - pRS [W], is the power consumption of the refrigeration system, this should be provided as data.
            %
            % Output:
            %   - p [W], is the active power consumption of the supermarket.
            %   - q [VAR], is the reactive power consumption of the supermarket.
            %   - timeReady [s], is the time until defrost cycles are ready to be executed.
            %   - timeDefrost [s], is the time until the defrost cycles will be executed by the supermarket.
            
            % Power and time flexibility from defrost cycle
            [pDF,timeReady,timeDefrost] = pCtrlDefrost(obj,exeDefrost);
            
            % Normalize output to per unit (PU)
            if obj.onPU == true
                p = -(pRS + pDF)/obj.sBase;
                q = p*tan(acos(obj.PF))/obj.sBase;
            else
                p = -(pRS + pDF);
                q = p*tan(acos(obj.PF));
            end
            
        end
        
        % Defrost cycles
        function [pDF,timeReady,timeDefrost] = pCtrlDefrost(obj,exeDefrost)
            % Input:
            %   - k [-], is the sample
            %   - exeDefrost [-], is the input that can start a defrost
            %     cycle. It is a vector with as many entries as
            %     controllable defrost cycles. A 1 in a entry indicates
            %     activation
            %
            % Output:
            %   - pDF [W], is the aggregated power consumption of all
            %     defrost cycles.
            %   - timeDefrost [s], is the time until a defrost is ready for
            %     activation. It is a vector with as many entries as
            %     controllable defrost cycles.
            
            % Get time to ready state and time to defrost
            [timeReady,timeDefrost] = timeToDefrost(obj);
            
            % Check defrost state and timers
            switchState(obj,exeDefrost);
            
            % Output power
            pVec = zeros(obj.numDF,1);
            for i=1:obj.numDF
                if obj.states(i) == 1
                    pVec(i) = obj.pDefrost(i);
                else
                    pVec(i) = 0;
                end
            end
            pDF = sum(pVec);
        end
        
        % Check time to ready and time to defrost
        function [timeReady, timeDefrost]= timeToDefrost(obj)
            timeReady = zeros(obj.numDF,1);
            timeDefrost = zeros(obj.numDF,1);
            for i=1:obj.numDF
                % Ready state
                if obj.states(i) == 0
                    timeReady(i) = 0;
                    timeDefrost(i) = obj.readyGuards(i)-obj.readyTimers(i);
                % Defrost state    
                elseif obj.states(i) == 1
                    timeReady(i) = obj.defrostGuards(i) + obj.waitGuards(i) - obj.defrostTimers(i);
                    timeDefrost(i) = obj.defrostGuards(i) + obj.waitGuards(i) + obj.readyGuards(i) - obj.defrostTimers(i);
                % Wait state
                elseif obj.states(i) == 2
                    timeReady(i) = obj.waitGuards(i) - obj.waitTimers(i);
                    timeDefrost(i) = obj.waitGuards(i) + obj.readyGuards(i) - obj.waitTimers(i);
                else
                    error('Something went wrong...')
                end
            end
        end
        
        % Switch defrost state and update timers
        function switchState(obj,u)
            % Update timers
            obj.readyTimers = obj.readyTimers + obj.Ts;
            obj.defrostTimers = obj.defrostTimers + obj.Ts;
            obj.waitTimers = obj.waitTimers + obj.Ts;
            
            % Check for state change
            for i=1:obj.numDF
                % In ready state
                if (obj.states(i) == 0 && obj.readyTimers(i) >= obj.readyGuards(i)) || (obj.states(i) == 0 && u(i) == 1)
                    obj.states(i) = 1;
                    obj.defrostTimers(i) = 0;
                % In defrost state
                elseif obj.states(i) == 1 && obj.defrostTimers(i) >= obj.defrostGuards(i)
                    obj.states(i) = 2;
                    obj.waitTimers(i) = 0;
                % In wait state    
                elseif obj.states(i) == 2 && obj.waitTimers(i) >= obj.waitGuards(i)
                    obj.states(i) = 0;
                    obj.readyTimers(i) = 0;
                else
                end
            end  
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
        
        % Set guard times
        function setGuards(obj,readyGuards,defrostGuards,waitGuards)
            obj.readyGuards = readyGuards;
            obj.defrostGuards = defrostGuards;
            obj.waitGuards = waitGuards;
        end
    end

end


