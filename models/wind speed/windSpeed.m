classdef windSpeed < handle
    %windSpeed implements a non-stationary stochastic wind speed model.
    %   The wind speed model can be used to investigate the behavior
    %   of e.g., wind turbine systems. The model is based on the book
    %   "Optimal Control of Wind Energy Systems 2008". 
    %
    %   The model simulates both the slow variations and the turbulence
    %   component, according to the Van der Hoven spectrum.
    %   
    % Revision:
    % 15-06-2014, C.E. Sloth and R. Pedersen, Aalborg University. Notes:
    
    properties
        Lt = 150;       % [m], is the turbulence length
        It = 1;         % [-], turbulence intencity
        Ts = 60;        % [sec], sampling time
        z = 30;         % [m], is the height of the turbine from ground
        z0 = 0.01;      % [m], is the terrains roughness length.
        
        % Shaping filter parameters
        m1 = 0.4;       % [-], parameter for the shaping filter
        m2 = 0.25;      % [-], parameter for the shaping filter
        xOld = [0 0];   % [-], placeholder for shaping filter
        
        % Van der Hoven Spectrum (For low frequency component)
        N = 50;         % [-], number of frequency points in the spectrum
        Svdh = [0.5 0.8 1.1 1.5 2.5 3.5 4.5 4.8 2.8 1.5 1.4 1.6 1.9 1.6 1 0.5 0.4 0.35 0.3 0.2];   % [(m/s)^2], Density data of spectrum 
        logx            % [cycles/hour], frequency point on a logarithmic scale
        Amp             % [-], Amplitude of harmonic frequencies
        phi             % [rad], phase of harmonic frequencies
    end
    
    methods
        % Constructor
        function obj = windSpeed(param)
            % Param format
            % - param.z             [m]. Height of turbines from ground 
            % - param.Ts            [s]. Sampling time
            
            % Set parameters
            obj.z = param.z;
            obj.Ts = param.Ts;
            
            % Set turbulence length (According to DS 742 2007)
            if obj.z >= 30
                obj.Lt = 150;
            elseif obj.z < 30
                obj.Lt = 5*obj.z;
            end
            
            % Set turbulence intencity (According to DS 742 2007)
            obj.It = 1/log(obj.z/obj.z0);
            
            % Setup Van der Hoven spectrum
            setupVanDerHovenSpectrum(obj);

        end
        
        % Sample model
        function [v, vs] = sample(obj,k,vMean)
            % Input:
            %   k [-], is the itteration number
            %   vMean [m/s], is the mean windspeed
            %
            % Output:
            %   v [m/s], is the wind speed.
            %   vs [m/s], is the mean wind speed without turbulence.
            
            % Low frequency component
            vs = sWindSpeed(obj,k,vMean);
            
            % Turbulence component
            [Kf, Tf] = staticGain(obj,vs);
            v = shapeFilter(obj,Tf,Kf);
            
            % Wind speed output 
            v = v*obj.It*vs + vs;
            if v<0
                v=0;
            end
            
        end
        
        % Low frequency component of the Van der Hoven spectrum
        function vs = sWindSpeed(obj,k,vMean)
            
            % Eq. (3.11) in the book
            vs = vMean;
            for i=1:obj.N
                vs = vs + obj.Amp(i)*cos(k*obj.Ts*obj.logx(i)...
                        *2*pi/(60*60) + obj.phi(i)); % Converted to [rad/s]
                                                     % from [cycles/hour]   
            end
        end
        
        % Shaping filter (Eq. (3.12) in book)
        function y = shapeFilter(obj,Tf,Kf)
            % Create filter 
            [A,B,C,D] = tf2ss([obj.m1*Tf*Kf Kf],[Tf*obj.m2 Tf+obj.m2*Tf 1]);
            sys = ss(A,B,C,D); % Needs to be on ss form to use initial condition in lsim
            
            % Run filter
            t = 0:obj.Ts:obj.Ts;
            e = randn(1,length(t));
            [y,~,x] = lsim(sys,e,t,obj.xOld);
            obj.xOld = x(end,:);

            y = y(end);
        end
        
        % Static gain computation (Eq. (3.9) in book)
        function [Kf, Tf] = staticGain(obj,vs)
            if vs <= 0
                Tf = 10000;
            else
                Tf = obj.Lt/vs;
            end
            Kf = sqrt(((2*pi)/(beta(1/2,1/3)))*(Tf/obj.Ts));
        end
        
        % Setup the Van der Hoven spectrum
        function setupVanDerHovenSpectrum(obj)
            % Setup logarithmic x coordinates for initial data points
            x = logspace(-3,0,length(obj.Svdh));
            
            % Convert density data to gains
            S = obj.Svdh./x;
            
            % Interpolate between points
            obj.logx = logspace(-3,0,obj.N+1);      
            Y = interp1(x,S,obj.logx,'spline');
            
            % Amplitude of the harmonic frequencies
            A = zeros(1,obj.N);
            for i = 1:obj.N
                A(i) = 2/pi*sqrt(1/2*(Y(i)+Y(i+1))*(obj.logx(i+1)-obj.logx(i)));
            end
            obj.Amp = A;
            
            % Set phase of harmonic frequencies
            obj.phi =  -pi + (pi+pi)*rand(1,obj.N);

        end
        
        %% Set/get functions
        function setNumVdhSpec(obj,N)
            if N<20
                obj.N = 20;
                warning('Minimum number of points in van der Hoven spectrum is 20. This has been changed to 20 in the simulation')
            else
                obj.N = N;
            end
            % Setup the Van der Hoven Spectrum when number of points
            % has changed
            setupVanDerHovenSpectrum(obj);
        end

        
    end
    
end
