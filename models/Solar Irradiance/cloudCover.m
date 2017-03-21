classdef cloudCover < handle
    %cloudCover implements a cloud cover model.
    %   The cloud cover model can be used to investigate the behavior
    %   of e.g., photovoltaic systems. This model can be used together with
    %   the solar irradiance model.
    %
    %   The model is based on the work by
    %   Heinriich Morf (1998-2013)
    %   
    %
    %	R. Pedersen 6-11-2014, Aalborg University
    
    properties 
        % For cloud cover model (cc)
        mean = 0.7511;  % [-], mean of cloud cover
        var = 0.1567;   % [-], variance of cloud cover
        N = 50;         % [-], number of observation points
        Aobs = 1;       % [m^2], observed area (can be set to 1 if unknown)
        A1_m = 0;       % [m^2], mean of covered area
        A0_m = 1;       % [m^2], mean of non-covered area
        P = zeros(2,1);   % [-], transition probability matrix
        
        ccOld;          % [-], array of placeholder for cc, N long
        n = 1;
    end
    
    methods
        % Constructor
        function obj = cloudCover(param)
            % Set parameters
            obj.mean = param.mean;
            obj.var = param.var;
            obj.N = param.N;
            obj.Aobs = param.Aobs;
            
            obj.ccOld = zeros(obj.N,1);
            obj.ccOld(1:2:obj.N) = 1;
            % Setup cloud cover model
            f = @(x) obj.var - 2*(1-obj.mean)^2*obj.mean^2*x ...
                *(1-(1-obj.mean)*obj.mean*x*(1-...
                exp(-1/((1-obj.mean)*obj.mean) * 1/x)));
            x = fzero(f,1);
            obj.A1_m = obj.mean*x;
            obj.A0_m = (obj.A1_m - obj.mean*obj.A1_m)/obj.mean;
            W1 = obj.mean;
            W0 = 1-obj.mean;
            % Transition probability matrix

            obj.P(1) = W1-W1*exp((-obj.Aobs/obj.N)*((obj.A0_m+obj.A1_m)/(obj.A0_m*obj.A1_m)));
            obj.P(2) = W0-W0*exp((-obj.Aobs/obj.N)*((obj.A0_m+obj.A1_m)/(obj.A0_m*obj.A1_m)));

        end
        
        % Sample model
        function cc = sample(obj)
            if obj.n > obj.N
                obj.n = 1;
            end
            r = rand(1);
            if obj.ccOld(obj.n) == 1
                if r < obj.P(1)
                    obj.ccOld(obj.n) = 0;
                else
                    obj.ccOld(obj.n) = 1;
                end
            elseif obj.ccOld == 0
                if r < obj.P(2)
                    obj.ccOld(obj.n) = 1;
                else
                    obj.ccOld(obj.n) = 0;
                end

            end
            cc = sum(obj.ccOld)/obj.N;
            obj.n = obj.n+1;   
        end

        
        %% Set functions
        % set mean and variance of cloud cover model
        function setMeanVar(obj,mean,var)
            obj.mean = mean;
            obj.var = var;
            % Setup model
            f = @(x) obj.var - 2*(1-obj.mean)^2*obj.mean^2*x ...
                *(1-(1-obj.mean)*obj.mean*x*(1-...
                exp(-1/((1-obj.mean)*obj.mean) * 1/x)));
            x = fzero(f,1);
            obj.A1_m = obj.mean*x;
            obj.A0_m = (obj.A1_m - obj.mean*obj.A1_m)/obj.mean;
            W1 = obj.mean;
            W0 = 1-obj.mean;
            % Transition probability matrix
            obj.P(1,1) = W0+W1*exp((-obj.Aobs/obj.N)*((obj.A0_m+obj.A1_m)/(obj.A0_m*obj.A1_m)));
            obj.P(1,2) = W1-W1*exp((-obj.Aobs/obj.N)*((obj.A0_m+obj.A1_m)/(obj.A0_m*obj.A1_m)));
            obj.P(2,1) = W0-W0*exp((-obj.Aobs/obj.N)*((obj.A0_m+obj.A1_m)/(obj.A0_m*obj.A1_m)));
            obj.P(2,2) = W1+W0*exp((-obj.Aobs/obj.N)*((obj.A0_m+obj.A1_m)/(obj.A0_m*obj.A1_m)));
        end

        
    end
    
end


