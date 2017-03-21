classdef comLink < handle
    %comLink implements a communication link between assets/controllers in the grid.
    % The comLink implements a communication link between different
    % entities in the grid simulation. That could be a controller and a 
    % asset or between two controllers or even between a controller and 
    % several assets. To generat delay and loss, the comLink uses inverse
    % transform sampling. This implies the need for inverse cumulative 
    % distribution functions.
    %
    %
    % Revision:
    % 16-10-2014, R. Pedersen, C. Sloth, Aalborg University. Notes:
    
    properties
        maxDelay                        % [-]. Maximum delay time in samples
        bufLength                       % [-]. Buffer length
        Ts                              % [s]. Sampling time
        numLinksIn                      % [-]. Number of communication links in
        numLinksOut                     % [-]. Number of communication links in
        inBuf                           % [-]. Buffer containing input data (e.g., measurements)
        inTimeBuf                       % [-]. Buffer containing time delay info
        outBuf                          % [-]. Buffer containing output data (e.e., control output)
        outTimeBuf                      % [-]. Buffer containing time delay info
        
        % Parameters for the inverse CDF functions
        % Delay
        iCDFdelay = 0;                  % [-]. Inverse CDF function for delay
        iCDFdelayNormMu = 1;            % [-]. Mean of delay normal distribution
        iCDFdelayNormSigma = 1;         % [-]. Standard deviation of delay normal distribution
        % Loss
        prLoss = 0.1;                   % [-]. Probability of packet loss
        
    end
    
    methods
        % Constructor
        function obj = comLink(param)
            % Param format
            % - param.maxDelay      [s]. Maximum delay time
            % - param.Ts            [s]. Sampling time
            % - param.numLinks      [-]. Number of communication links 
            
            % Set parameters
            obj.maxDelay = ceil(param.maxDelay/param.Ts);
            obj.Ts = param.Ts;
            obj.numLinksIn = param.numLinksIn;
            obj.numLinksOut = param.numLinksOut;
            
            % Compute parameters
            obj.bufLength = ceil(param.maxDelay/param.Ts);
            obj.inBuf = ones(obj.bufLength,obj.numLinksIn);
            obj.inTimeBuf = zeros(obj.bufLength,obj.numLinksIn);
            obj.outBuf = ones(obj.bufLength,obj.numLinksOut);
            obj.outTimeBuf = zeros(obj.bufLength,obj.numLinksOut);
        end
        
        % Sample the communication link incomming
        function out = sampleIn(obj,k,in,delayData)
            % Input:
            %   - k [-], is the sample number
            %   - in [-], is the incomming data
            %   - delayData [s], used when delay profiles are provided by
            %   data
            %
            % Output:
            %   - out [-], is the data output with delay and loss
            
            % Check if delayData is provided
            if ~strcmp(obj.iCDFdelay,'dataTrace')
                delayData = 0;
            end
            
            % Change packets which are lost to NaN
            in = loss(obj,in,obj.numLinksIn);
            
            % Generate delay    
            % Put incomming data into buffer
            obj.inBuf(mod(k-1,obj.bufLength)+1,:) = in;
            % Determin delay
            obj.inTimeBuf(mod(k-1,obj.bufLength)+1,:) = k+delay(obj,obj.numLinksIn,delayData);
            
            % Output data 
            [~,IndexIn] = max(obj.inTimeBuf.*(obj.inTimeBuf<=k),[],1);
            for j = obj.numLinksIn:-1:1
                out(j) = obj.inBuf(IndexIn(j),j);
            end
        end
        
        % Sample the communication outgoing
        function out = sampleOut(obj,k,in,delayData)
            % Input:
            %   - k [-], is the sample number
            %   - in [-], is the incomming data
            %
            % Output:
            %   - out [-], is the data output with delay and loss
            
            % Check if delayData is provided
            if ~strcmp(obj.iCDFdelay,'dataTrace')
                delayData = 0;
            end
            
            % Change packets which are lost to NaN
            in = loss(obj,in,obj.numLinksOut);
            
            % Generate delay    
            % Put incomming data into buffer
            obj.outBuf(mod(k-1,obj.bufLength)+1,:) = in;
            % Determin delay
            obj.outTimeBuf(mod(k-1,obj.bufLength)+1,:) = k+delay(obj,obj.numLinksOut,delayData);
            
            % Output data
            [~,IndexOut] = max(obj.outTimeBuf.*(obj.outTimeBuf<=k),[],1);
            for j = obj.numLinksOut:-1:1
                out(j) = obj.outBuf(IndexOut(j),j);
            end
        end
        
        % Generate delay
        function out = delay(obj,numLinks,delayData)
            % Input:
            %   - numLinks [-], is the "length" or "number" of data...
            %   - delayData [s], if used is the current delay in seconds of
            %   data. This is recast to samples.
            %
            % Output:
            %   - out [-], delay of data in samples.
            
            % You can put in any desired inverse CDF or base in on data
            % traces.
            switch obj.iCDFdelay
                case 'normal'
                    out = ceil(norminv(ones(1,numLinks)*rand(),obj.iCDFdelayNormMu,obj.iCDFdelayNormSigma));
                    out = out.*(out>0);
                case 'uniform'
                    out = ceil(obj.maxDelay*ones(1,numLinks)*rand());
                case 'dataTrace'
                    delayDataInSamples = ceil(delayData/obj.Ts);
                    if delayDataInSamples >= obj.maxDelay
                        delayDataInSamples = obj.maxDelay;
                    end
                    out = delayDataInSamples*ones(1,numLinks);
                otherwise % Uniform
                    out = ceil(obj.maxDelay*ones(1,numlinks)*rand());
            end
        end
        
        % Generate loss
        function out = loss(obj,in,numLinks)
            % Replace lost data with NaN
            idx = rand(1,numLinks)<obj.prLoss;
            out = in;
            out(idx) = NaN;
        end
        
        %% Set and get functions
        % Set inverse CDF for delay
        function setInverseCDFdelay(obj,iCDFFun,param)
            obj.iCDFdelay = iCDFFun;
            obj.iCDFdelayNormMu = param.mu/obj.Ts;
            obj.iCDFdelayNormSigma = param.sigma/obj.Ts;
            obj.maxDelay = param.maxDelay/param.Ts;
        end
        
        % Set probability of loss
        function setPrLoss(obj,prLoss)
            obj.prLoss = prLoss;
        end
    end
end
