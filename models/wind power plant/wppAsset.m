classdef wppAsset < handle
    %wppAsset implements a wind power plant.
    %  The wind power plant (WPP) is a collection of wind turbine asset
    %  objects (wtAsset, see class file). The wind power plant simply
    %  distributes the references on both active and reactive power between
    %  the wind turbines. Wake effekt and displacement between wind
    %  turbines are not considered.
    %   
    %   
    % Revision:
    % 08-10-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
        N                           % [-]. Number of wind turbines in the WPP
        pOut                        % [W]. Vector of active power outputs
        qOut                        % [VAR]. Vector of reactive power outputs
        pAva                        % [W]. Vector of available power outputs
        wsMod                       % [-]. Flag used to limit sampling frequency of wind speed.
        wsOld = 0;                  % [m/s]. Placeholder for old wind speed sample.
        
        % Wind speed and wind turbine objects
        ws                          % [-]. Wind speed object
        wt                          % [-]. Set of windturbine objects
        
    end
    
    methods
        % Constructor
        function obj = wppAsset(param)
            % Parameters are a collection for a wind turbine asset and 
            % wind speed object. It also has a parameter indicating the
            % number of turbines in the wind power plant. The parameter
            % structure is as follows:
            % 
            %   param.sBase         % Complex power base of the WPP
            %   param.vBase         % Voltage base of the WPP
            %   param.pRated        % Rated power of each wind turbine
            %   param.sMax          % Maximum apparent power of each wind turbine
            %   param.wMin          % Cut-in wind speed of each wind turbine
            %   param.wMax          % Cut-out wind speed of each wind turbine
            %   param.wRated        % Rated wind speed of each wind turbine
            %   param.Ts            % Sampling time
            %   param.onPU          % Indicator if the system is simulated on a per unit basis
            %   param.numWt         % Number of wind turbines in the WPP
            %   param.z             % Is the height of each wind turbine from ground
            %   param.avgFiltTime   % Is the average filter on wind speed time constant
             
            
            % Construct wind speed object and wind turbine objects and allocate
            % memory
            obj.N = param.numWt;
            obj.pOut = zeros(obj.N,1);
            obj.qOut = zeros(obj.N,1);
            obj.pAva = zeros(obj.N,1);
            
            % Wind speed object;
            obj.ws = windSpeed(param);
            obj.wt = wtAsset(param);
            for i=1:obj.N
                %obj.ws(i) = windSpeed(param);
                obj.wt(i) = wtAsset(param);
                obj.wt(i).avgFilterTime = param.avgFiltTime;
            end  
            
            % Reduce the sampling frequency for the wind speed, if Ts is to
            % low
            obj.wsMod = ceil(20/param.Ts);
        end
        
        
        % sample wind turbine
        function [p,q,pAva] = sample(obj,wMean,v,dP,dPlim,qRef,vRef,k)
            % Input:
            %   - wMean, is the mean wind speed.
            %   - v, is the voltage at the connection point.
            %   - dP, is the reference to active power change.
            %   - dPlim, is the reference to derated power.
            %   - qRef, is the reactive power reference.
            %   - vRef, is the voltage reference.
            %   - k, is the sample number
            
            % Sample wind speed and wind turbine objects
            for i=1:obj.N
                if ~mod(k,obj.wsMod)
                     wSpeed = obj.ws.sample(k,wMean);
                else
                    wSpeed = obj.wsOld;
                end
%                 wSpeed = obj.ws.sample(k,wMean);
                [obj.pOut(i),obj.qOut(i),obj.pAva(i)] = ...
                    obj.wt(i).sample(wSpeed,v,dP/obj.N,dPlim/obj.N,qRef/obj.N,vRef);
            end
            obj.wsOld = wSpeed;
            
            p = sum(obj.pOut);
            q = sum(obj.qOut);
            pAva = sum(obj.pAva);

        end
        
        %% Set and get functions
        function setPF(obj,PF)
            if PF<0
                PF = 0;
            elseif PF>1;
                PF = 1;
            end
            for i=1:obj.N
                obj.wt(i).setPF(PF);
            end
        end
        
        % Set Qmode and Q control function
        function setQmode(obj,Qmode)
            % Input:
            %   - Qmode=0: Power factor is kept constantly equal to PF.
            %   - Qmode=1: Reactive power is to qRef [PU].
            %   - Qmode=2: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
            
            % Set mode
            for i=1:obj.N
                obj.wt(i).setQmode(Qmode);
            end
        end
        
    end
    
end

