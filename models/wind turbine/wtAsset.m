classdef wtAsset < handle
    %wtAsset implements a simple P,Q controllable wind turbine system.
    %   The wtAsset can be set to follow an active power reference, output
    %   maximum power available or curtail power. It can also, at the same
    %   time be set to follow a reactive power reference, output the
    %   minimum reactive power according to the power factor or try to keep
    %   a voltage reference in the coupling point by controlling reactive
    %   power output.
    %   
    %   The model is partly based on the book
    %   "Optimal Control of Wind Energy Systems 2008".
    %   
    %   
    % Revision:
    % 23-05-2014, C. E. Sloth and R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 3e6;                % [VA]. Complex power base
        vBase = 20e3;               % [V]. Voltage base
        pRated = 3e6;               % [W]. Rated output power
        sMax = 3e6;                 % [VA]. Maximum apparent power
        wMin = 3;                   % [m/s]. Cut-in wind speed
        wMax = 25;                  % [m/s]. Cut-out wind speed
        wRated = 12;                % [m/s]. Rated wind speed
        kWT = 3e6/(12^3);           % [-]. Gain from w^3 to p
        PF = 1;                     % [-]. Power factor
        Qmode = 0;                  % [-]. Reactive power control mode
        Ts = 60;                    % [s]. Sampling time
        onPU = false;               % [-]. Indicates if the wind turbine system is used to simulate a system on per unit basis. true = yes, false = no    
        
        % Parameters
        rho = 1.225;                % [kg/m^3]. Is the density of air (at 15 degrees C)
        bLimit = 0.593;             % [-]. Betz limit
        eta = 0.9;                  % [-]. Effiency of the wind power system
        r                           % [m]. Radius of swept area
        xOldsf = [0 0];             % [-]. Place holder spatial filter
        wOld = 0;                   % [m/s]. Place holder for the old wind input
        pOld = 0;                   % [w]. Place holder for the rotor initia filter
        rfTau = 10;                 % [s]. Time constant of rotor filter
        avgFilterBuf                % [-]. Buffer for the moving average filter used for start and stop
        avgFilterTime = 10*60;      % [s]. Time that the filter should look back 
    end
    
    methods
        % Constructor
        function obj = wtAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.vBase         [V]. Voltage base
            % - param.sMax          [VA]. Is the maximum apparent power
            % - param.pRated        [W]. Is the rated power
            % - param.wMin          [m/s]. Cut-in wind speed
            % - param.wMax          [m/s]. Cut-out wind speed
            % - param.wRated        [m/s]. Rated wind speed
            % - param.Ts            [s]. Sampling time
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.onPU = param.onPU;
            obj.vBase = param.vBase;
            obj.pRated = param.pRated;
            obj.sMax = param.sMax;
            obj.wMin = param.wMin;
            obj.wMax = param.wMax;
            obj.wRated = param.wRated;
            obj.Ts = param.Ts;
            % Set sBase if simulating on a per unit basis
            if obj.onPU == true
                obj.sBase = param.sBase;
            end
            
            % Power conversion constant
            obj.kWT = param.pRated/(param.wRated^3);
            
            % Calculate the radius of the swept area
            cpMax = obj.bLimit*obj.eta;
            A = obj.pRated/(0.5*obj.rho*cpMax*obj.wRated^3);
            obj.r = sqrt(A/pi);
            
            % Setup buffer for the moving average filter
            obj.avgFilterBuf = zeros(ceil(obj.avgFilterTime/obj.Ts),1);
            
        end
        
        
        % sample wind turbine
        function [p,q,pAva] = sample(obj,w,v,dP,dPlim,qRef,vRef)
            % Input:
            %   - w [m/s], is the wind speed.
            %   - v [V], is the voltage at the connection point.
            %   - dP [W], is the reference to active power change (curtailment).
            %   - dPlim [W], is the reference to derated power (derate).
            %   - qRef [VAR], is the reactive power reference.
            %   - vRef [V], is the voltage reference.
            %
            % Output:
            %   - p [W], is active power outpur.
            %   - q [VAR], is reactive power output.
            %   - pAva [W], is the available active power.
            
            % Convert References to Physical Quantities
            if obj.onPU == true
                dP = dP*obj.sBase;
                dPlim = dPlim*obj.sBase;
                qRef = qRef*obj.sBase;
                v = v*obj.vBase;
                vRef = vRef*obj.vBase;
            end
            
            % Spatial filter
            wt = spatialFilter(obj,w);
            
            % Active power control
            [p,pAva] = pCtrl(obj,wt,dP,dPlim);
            
            % Reactive power control
            q = qCtrl(obj,p,v,qRef,vRef);
            
            p = rotorFilter(obj,p);
            % Normalize Output to Per Unit (PU)
            if obj.onPU == true
                p = p/obj.sBase;
                q = q/obj.sBase;
            end

        end
        
        % Active power control
        function [p,pAva] = pCtrl(obj,w,dP,dPlim)
            % Calculate Active Power Output
            pW = obj.kWT*obj.PF*w^3;    % [W] Available power
            % Output the available power of the wind turbine
            if pW >= obj.pRated
                pAva = obj.pRated;
            elseif pW <= 0
                pAva = 0;
            else
                pAva = pW;
            end

            dP = min(dP,dPlim); % [W] dP is limited by dPlim
            
            % Get average of wind for start and stop of turbine
            w = avgFilter(obj,w);

            % Calculate gain from wind disturbance to power output
            if(w<obj.wMin || w>obj.wMax) % Out of normal operating range [wMin wMax] (in reality these should be 10 min mean values)
                p = 0;
            elseif(pW>=obj.pRated+dP-dPlim) % In full load operation
                p = max(0,obj.pRated+dP-dPlim); % The output power must be positive
            else   
                dP = min(dP,0); % In partial load operation up-regulation is not possible
                p = max(0,pW+dP); % The output power must be positive
            end

        end
        
        % Reactive power control
        function q = qCtrl(obj,p,v,qRef,vRef)
            % Calculate Reactive Power Reference
            if(obj.Qmode == 0)
                qRef = p*tan(acos(obj.PF));
            elseif(obj.Qmode == 1)
                qRef = qRef;
            elseif(obj.Qmode == 2)
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
            if obj.onPU == true
                vMax = 1.1;
                vMin = 0.9;
                k = -2/(vMax-vMin);
                qRef = obj.sMax*k*(v-vRef)/obj.vBase;
            else
                vMax = obj.vBase * 1.1;
                vMin = obj.vBase * 0.9;
                k = -2/(vMax-vMin);
                qRef = obj.sMax*k*(v-vRef);
            end
        end
        
        % Spatial low pass filter
        function wt = spatialFilter(obj,w)
            if w<=0
                w= 0.001;
            end
                
            % This filter averages the wind speed over the area swept by
            % the rotor.
            a = 30;         % [-], is an empirical factor between 0 and 55
            gamma = 1.3;    % [-], is an empirical factor
            b = gamma*(obj.r/w);
            
            % Create filter
            [A,B,C,D] = tf2ss([b sqrt(2)],[b^2 sqrt(2)*b/sqrt(a)+b*sqrt(a) sqrt(2)]);
            sys = ss(A,B,C,D); % Needs to be on ss form to use initial conditions in lsim
            
            % Run filter
            t = 0:obj.Ts:obj.Ts;
            [y,~,x] = lsim(sys,[obj.wOld w],t,obj.xOldsf);
            obj.xOldsf = x(end,:);
            wt = y(end);
            obj.wOld = w;
        end
        
        % Low pass filter for simulating initia of rotor
        function p = rotorFilter(obj,u)
            p = u*(obj.Ts/(obj.rfTau+obj.Ts)) + obj.pOld*(obj.rfTau/(obj.rfTau+obj.Ts));
            obj.pOld = p;
        end
        
        % Moving average filter for the start and stop procidure of the
        % wind turbine
        function wOut = avgFilter(obj,w)
            % Update data
            obj.avgFilterBuf(2:end) = obj.avgFilterBuf(1:end-1);
            obj.avgFilterBuf(1) = w;
            wOut = sum(obj.avgFilterBuf)/length(obj.avgFilterBuf);
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
        
        % Set Qmode and Q control function
        function setQmode(obj,Qmode)
            % Input:
            %   - Qmode=0: Power factor is kept constantly equal to PF.
            %   - Qmode=1: Reactive power is to qRef [PU].
            %   - Qmode=2: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
            
            % Set mode
            obj.Qmode = Qmode;
        end
    end
    
end

