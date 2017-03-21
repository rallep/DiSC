classdef pvAsset < handle
    %pvAsset implements a simple P,Q controllable solar PV system.
    %   The pvAsset can be set to follow an active power reference, output
    %   maximum power available or derate power. It can also, at the same
    %   time be set to follow a reactive power reference, output the
    %   minimum reactive power according to the power factor or try to keep
    %   a voltage reference in the coupling point by controlling reactive
    %   power output.
    %   
    %   
    % Revision:
    % 19-06-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
        sBase = 6e3;                % [VA]. Complex power base
        vBase = 400;                % [V]. Voltage base
        pRated = 6e3;               % [W]. Rated output power
        sMax = 6e3;                 % [VA]. Maximum apparent power
        onPU = false;               % [-]. Indicates if the PV system is used to simulate a system on per unit basis. true = yes, false = no
        Ts = 1;                     % [s]. System sampling time
        A = 10;                     % [m^2]. Area of PV cell
        eta = 0.2;                  % [-]. Efficency of cell
        PF = 1;                     % [-]. Power factor
        Qmode = 0;                  % [-]. Reactive power control mode
        Pmode = 0;                  % [-]. Active power control mode
        droopQFun = 0;              % [-]. Flag indicating which droop control function on reactive power should be used
        kPV = 0.2*10;               % [-]. Gain from solar to electrical power
        
        % For voltage droop control on reactive and active power
        avgFilterBuf                % [-]. Buffer for moving average filter on voltage measurement. Used for local voltage control
        avgFilterLength = 5;        % [-]. Length of moving average filter in samples
        droopQGain = 1;             % [-]. Gain on droop controller on reactive power
        droopPGain = 1;             % [-]. Gain on droop controller on active power
        vLPFOld = 0;                % [V]. Former value of the low-pass filter on voltage magnitude
        vLPFTau = 5;                % [s]. Time constant of low-pass filter of voltage magnitude
        vDbLow = 0.99;              % [PU]. Lower voltage deadband
        vDbHigh = 1.01;             % [PU]. High voltage deadband
        vHyst = 0.005;              % [PU]. Voltage hysteresis.
        vP = 1.05;                  % [PU]. Voltage limit for when to curtail active power (used in droop control on active power)
        vMax = 1.1;                 % [PU]. Maximum allowable voltage deviation
        actLow = 0;                 % [-]. Flag indicating if droop control is active in the low voltage region (used for voltage droop with hystwrisis reactive power)
        actHigh = 0;                % [-]. Same as above, but for high voltage region
        filtType = 0;               % [-]. Flag determining how the voltage should be filtered. 0 = no filter (Default), 1 = first order low-pass filter, 2 = moving avereging filter
    end
    
    methods
        % Constructor
        function obj = pvAsset(param)
            % Param format
            % - param.sBase         [VA]. Complex power base. Used for simulating on PU basis
            % - param.vBase         [V]. Voltage base
            % - param.sMax          [VA]. Is the maximum apparent power
            % - param.pRated        [W]. Is the rated power 
            % - param.eta           [0-1]. Solar cell efficiency
            % - param.A             [m^2]. Area of solar cells 
            % - param.onPU          [-]. Flag indicating if simulation is on a per unit basis, (true or false)
            
            % Set parameters
            obj.vBase = param.vBase;
            obj.pRated = param.pRated;
            obj.sMax = param.sMax;
            obj.Ts = param.Ts;
            obj.eta = param.eta;
            obj.A = param.A;
            % Set sBase if simulating on a per unit basis
            if obj.onPU == true
                obj.sBase = param.sBase;
            end
            
            % Solar to power conversion factor
            obj.kPV = obj.eta*obj.A;
            
            % Setup droop gains on P and Q
            obj.droopQGain = -obj.sMax/(obj.vP*obj.vBase-obj.vBase);
            obj.droopPGain = -obj.sMax/(obj.vMax*obj.vBase-obj.vP*obj.vBase);
            
            % Setup moving average filter for voltage measurement
            obj.avgFilterBuf = ones(obj.avgFilterLength,1)*obj.vBase;
            obj.vLPFOld = param.vBase;
            
        end
        
        % sample solar PV system
        function [p,q,pAva] = sample(obj,s,v,dP,dPlim,qRef,vRef)  
            % Input:
            %   - s [w/m^2], is the solar irradiance.
            %   - v [V], is the voltage at the connection point.
            %   - dP [W], is the reference to active power change.
            %   - dPlim [W], is the reference to derated power.
            %   - qRef [VAR], is the reactive power reference.
            %   - vRef [V], is the voltage reference.
            %
            % Output:
            %   - p [W], is the active power output
            %   - q [VAR], is the reactive power output
            %   - pAva [W], is the available power
            
            % Filter voltage measurement
           switch obj.filtType
                case 1
                    vOut = voltLowPassFilter(obj,v);
                case 2
                    vOut = avgFilter(obj,v);
                otherwise
                    vOut = v;
            end
            
            % Active power control
            [p,pAva] = pCtrl(obj,s,dP,dPlim,vOut,vRef);
            
            % Reactive power control
            q = qCtrl(obj,p,vOut,qRef,vRef);
            
            if obj.onPU == true
                % Normalize Output to Per Unit (PU)
                p = p/obj.sBase;
                q = q/obj.sBase;   
            end

        end
        
        % Active power control
        function [p,pAvb] = pCtrl(obj,s,dP,dPlim,v,vRef)
            % Calculate Active Power Output
            pAvb = obj.kPV*obj.PF*s;      % [W] Available power
            switch obj.Pmode
                case 1
                    kp = obj.droopPGain;
                    if v>=obj.vP*obj.vBase
                        pCur = (v-obj.vP*obj.vBase)*kp;
                    else
                        pCur = 0;
                    end
                    
                    pRef = pAvb + pCur;
                    p = max(0,pRef);
                    p = min(p,obj.pRated);
                    
                otherwise % default is external control of curtailment and derating
                    dP = min(dP,dPlim);         % [W] dP is limited by dPlim

                    % Calculate gain from solar irradinace to power output
                    if(pAvb>=obj.pRated-dPlim) % In full load operation
                        p = max(0,obj.pRated+dP-dPlim);
                        p = min(p,obj.pRated);
                    else   
                        dP = min(dP,0); % In partial load operation up-regulation is not possible
                        p = max(0,pAvb+dP); % The output power must be positive
                    end
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
                qRef = droop(obj,v,vRef,obj.droopQFun);
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
                case 1
                    qRef = userDefined1(obj,v,vRef);
                case 2
                    qRef = userDefined2(obj,v,vRef);
                case 3
                    qRef = userDefined3(obj,v,vRef);
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
            %       - v is the bus voltage [V].
            %       - vRef is the voltage reference [V].
            %
            %   Output:
            %       - qRef is the reactive power reference [VAR].
            if obj.onPU == true
                vMax_local = 1.1;
                vMin_local = 0.9;
                k = -2/(vMax_local-vMin_local);
                qRef = obj.sMax*k*(v-vRef)/obj.vBase;
            else
                vMax_local = obj.vBase * 1.1;
                vMin_local = obj.vBase * 0.9;
                k = -2/(vMax_local-vMin_local);
                qRef = obj.sMax*k*(v-vRef);
            end
        end
        
        % userDefined1 droop control function (Linear with adjustable gain)
        function qRef = userDefined1(obj,v,vRef)
            % defaultQDroop is a reactive power droop control.
            %
            %   Input:
            %       - v is the bus voltage [V].
            %       - vRef is the voltage reference [V].
            %
            %   Output:
            %       - qRef is the reactive power reference [VAR].
            k = obj.droopQGain;
            qRef = k*(v-vRef);
        end
        
        % userDefined2 droop control function (With deadband)
        function qRef = userDefined2(obj,v,vRef)
            %   Input:
            %       - v is the bus voltage [PU].
            %       - vRef is the voltage reference [PU].
            %
            %   Output:
            %       - qRef is the reactive power reference [VAR].
            
            % Set droop gain
            k = obj.droopQGain;
            % Put voltages on PU basis
            v = v/obj.vBase;
            vRef = vRef/obj.vBase;
            % Check if voltage is in deadband
            if v<=obj.vDbLow
                qRef = k*(v-vRef-(obj.vDbLow-1))*obj.vBase;
            elseif v>obj.vDbLow && v<obj.vDbHigh
                qRef = 0;
            else
                qRef = k*(v-vRef-(obj.vDbHigh-1))*obj.vBase;
            end
        end
        
        % userDefined3 droop control function (With deadband and hysteresis)
        function qRef = userDefined3(obj,v,vRef)
            % defaultQDroop is a reactive power droop control.
            %
            %   Input:
            %       - v is the bus voltage [PU].
            %       - vRef is the voltage reference [PU].
            %
            %   Output:
            %       - qRef is the reactive power reference [VAR].
            
            % Set droop gain
            k = obj.droopQGain;
            % Put voltages on er PU basis
            v = v/obj.vBase;
            vRef = vRef/obj.vBase;
            % Check if droop control should be activated
            if v <= obj.vDbLow && obj.actLow ~= 1
                obj.actLow = 1;
            elseif v > obj.vDbLow+obj.vHyst && obj.actLow == 1
                obj.actLow = 0;
            elseif v >= obj.vDbHigh && obj.actHigh ~= 1
                obj.actHigh = 1;
            elseif v < obj.vDbHigh-obj.vHyst && obj.actHigh == 1
                obj.actHigh = 0;
            end
            % Apply voltage control
            if obj.actLow == 1
                qRef = k*(v-vRef-(obj.vDbLow-1))*obj.vBase;
            elseif obj.actHigh == 1
                qRef = k*(v-vRef-(obj.vDbHigh-1))*obj.vBase;
            else
                qRef = 0;
            end
        end
        
        % Moving average filter for the voltage measurement
        function vOut = avgFilter(obj,v)
            % Update data
            obj.avgFilterBuf(2:end) = obj.avgFilterBuf(1:end-1);
            obj.avgFilterBuf(1) = v;
            vOut = sum(obj.avgFilterBuf)/length(obj.avgFilterBuf);
        end
        
        % First order low pass filter for voltage
        function vOut = voltLowPassFilter(obj,u)
            vOut = u*(obj.Ts/(obj.vLPFTau+obj.Ts)) + obj.vLPFOld*(obj.vLPFTau/(obj.vLPFTau+obj.Ts));
            obj.vLPFOld = vOut;
        end
        
        
        %% Set and get functions
        % Set power factor
        function setPF(obj,PF)
            if PF<0
                PF = 0;
            elseif PF>1;
                PF = 1;
            end
            obj.PF = PF;
        end
        
        % Set Qmode (Reactive power control mode)
        function setQmode(obj,Qmode)
            % Input:
            %   - Qmode=0: Power factor is kept constantly equal to PF.
            %   - Qmode=1: Reactive power is to qRef [PU].
            %   - Qmode=2: Reactive power is controlled as q = qFun(v,vRef,sMax) [VAR]
            
            % Set mode
            obj.Qmode = Qmode;
        end
        
        % Set Pmode (Active power control mode)
        function setPmode(obj,Pmode)
            % Input:
            %   - Pmode=0: External curtailment/derating
            %   - Qmode=1: Internal droop on active power
            
            % Set mode
            obj.Pmode = Pmode;
        end
        
        % Setup reactive droop control function
        function setDroopQFun(obj,droopFun)
            obj.droopQFun = droopFun; 
        end
        
        % Set P/Q droop gains
        function setDroopGains(obj,qGain,pGain)
            obj.droopQGain = qGain;
            obj.droopPGain = pGain;
        end
        
        % Set voltage filter type
        function setVoltFilter(obj,filtType,LPFTau,filtBufLen)
            obj.filtType = filtType;
            if obj.filtType == 1
                obj.vLPFTau = LPFTau;
                obj.vLPFOld = obj.vBase;
            elseif obj.filtType == 2
                obj.avgFilterLength = filtBufLen;
                obj.avgFilterBuf = ones(obj.avgFilterLength,1)*obj.vBase;
            end
        end
        
        % Set voltage point for curtailment of active power
        
        
    end
    
end

