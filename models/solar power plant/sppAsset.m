classdef sppAsset < handle
    %sppAsset implements a solar power plant.
    %  The solar power plant (SPP) is a collection of PV asset
    %  objects (pvAsset, see class file), i.e., it is a solar farm.
    %  The solar power plant simply distributes the references on both
    %  active and reactive power between the PV systems.
    %   
    %   
    % Revision:
    % 15-10-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
        N                           % [-]. Number of PV systems in the SPP
        pOut                        % [W]. Vector of active power outputs
        qOut                        % [VAR]. Vector of reactive power outputs
        pAva                        % [W]. Vector of available power outputs
        
        % Wind speed and wind turbine objects
        si                          % [-]. Solar irradiance object
        pv                          % [-]. Set of PV system objects
        
    end
    
    methods
        % Constructor
        function obj = sppAsset(param)
            % Parameters are a collection for a PV asset and 
            % solar irradiance object. It also has a parameter indicating the
            % number of PV systems in the solar power plant. The parameter
            % structure is as follows:
            % 
            %   param.sBase         % Complex power base of the WPP
            %   param.vBase         % Voltage base of the WPP
            %   param.pRated        % Rated power of each wind turbine
            %   param.sMax          % Maximum apparent power of each wind turbine
            %   param.eta           % Efficiency of each PV system cells
            %   param.A             % Area of each PV system
            %   param.Ts            % Sampling time
            %   param.onPU          % Indicator if the system is simulated on a per unit basis
            %   param.numPv         % Number of PV systems in the SPP
            %   param.lat           % Latitudal location of the SPP
            %   param.t             % Transmittance of the location
            %   param.p             % Air pressure at the location
            
            
            % Construct solar irradiance object and
            % PV objects and allocate memory
            obj.N = param.numPv;
            obj.pOut = zeros(obj.N,1);
            obj.qOut = zeros(obj.N,1);
            obj.pAva = zeros(obj.N,1);
            
            % Wind speed object;
            obj.si = solarIrradiance(param);
            obj.pv = pvAsset(param);
            for i=1:obj.N
                %obj.ws(i) = windSpeed(param);
                obj.pv(i) = pvAsset(param);
            end         
        end
        
        
        % sample wind turbine
        function [p,q,pAva] = sample(obj,cc,v,dP,dPlim,qRef,vRef,k,day)
            % Input:
            %   - cc, is the cloud cover [0,1].
            %   - v, is the voltage at the connection point.
            %   - dP, is the reference to active power change.
            %   - dPlim, is the reference to derated power.
            %   - qRef, is the reactive power reference.
            %   - vRef, is the voltage reference.
            %   - k, is the sample number
            %
            % Output:
            %   - p [W], is the active power output
            %   - q [VAR], is the reactive power output
            %   - pAva [W], is the available power
            
            % Sample solar irradiance and PV systems
            for i=1:obj.N
                sIrradiance = obj.si.sample(k,day,cc);
                [obj.pOut(i),obj.qOut(i),obj.pAva(i)] = ...
                    obj.pv(i).sample(sIrradiance,v,dP/obj.N,dPlim/obj.N,qRef/obj.N,vRef);
            end
            
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
                obj.pv(i).setPF(PF);
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
                obj.pv(i).setQmode(Qmode);
            end
        end
        
        % Set Pmode and P control function
        function setPmode(obj,Pmode)
            % Input:
            %   - Qmode=0: Power factor is kept constantly equal to PF.
            %   - Qmode=1: Reactive power is to qRef [PU].
            %   - Qmode=2: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
            
            % Set mode
            for i=1:obj.N
                obj.pv(i).setPmode(Pmode);
            end
        end
        
        % Set reactive droop function
        function setDroopQFun(obj,droopFun)
            for i=1:obj.N
                obj.pv(i).setDroopQFun(droopFun);
            end
        end
        
        % Set P/Q droop gains
        function setDroopGains(obj,qGain,pGain)
            for i=1:obj.N
                obj.pv(i).setDroopGains(qGain,pGain);
            end
        end
        
        % Set voltage filter type and parameters
        function setVoltFilter(obj,filtType,LPFTau,filtBufLen)
            for i=1:obj.N
                obj.pv(i).setVoltFilter(filtType,LPFTau,filtBufLen);
            end
        end
        
    end
    
end

