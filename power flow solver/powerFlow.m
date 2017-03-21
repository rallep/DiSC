classdef powerFlow < handle
    %powerFlow implements solving power flow and setting up the grid.
    % The powerFlow object implements a Newton-Raphson load flow solver 
    % along with tools converting the grid line data into the impedance
    % matrix.
    %
    %
    % Revision:
    % 16-10-2014, R. Pedersen, Aalborg University. Notes:
    
    properties
         maxIte = 100;              % [-]. Maximum number of iterations
         nIte = 0;                  % [-]. Number of iterations
         tol = 1e-6;                % [-]. Tolerance
         type                       % [-]. Type of each bus (0=slack, 1=PQ, 2=PV)
         vBase                      % [V]. Voltage base of all busses
         
         % Parameters for NR
         nBus                       % [-]. Number of busses
         nSlack                     % [-]. Number of slack busses
         nPQ                        % [-]. Number of PQ busses
         nPV                        % [-]. Number of PV busses
         idSlack                    % [-]. Index of slack busses
         idPQ                       % [-]. Index of PQ busses
         idPV                       % [-]. Index of PV busses
         idPQV                      % [-]. Index of both PQ and PV busses
         Vm                         % [V]. Voltage magnitude
         Va                         % [deg]. Voltage angel
         
         % Parameters for calculating line loading
         
    end
    
    methods
        % Constructor
        function obj = powerFlow(param)
            % Param format
            % - param.type          [-]. Vector containing the type of each bus
            % - param.vBase         [V]. Vector containing the base voltage of each bus

            % Set parameters
            obj.type = param.type;
            obj.vBase = param.vBase;
            
            % Compute parameters
            obj.nBus = length(param.type);                      % Number of busses
             
            % Setup index parameters
            obj.nSlack = obj.nBus - length(find(obj.type));     % Number of slack busses
            obj.nPV = length(find(obj.type>1));                 % Number of PV busses
            obj.nPQ = obj.nBus-obj.nSlack-obj.nPV;              % Number of PQ busses

            obj.idSlack = find(obj.type==0);                    % Index of slack busses
            obj.idPQ = find(obj.type==1);                       % Index of PQ busses
            obj.idPV = find(obj.type==2);                       % Index of PV busses
            obj.idPQV = find(obj.type>0);                       % Index of both PQ and PV busses

            % Initial guess on voltage magnitudes and angles                           (STEP 1)
            obj.Vm = zeros(obj.nBus,1);                         % Vector of bus voltage magnitudes  
            obj.Va = zeros(obj.nBus,1);                         % Vector of bus voltage angles

            % Set slack bus voltages in Vm and Va vector
            for k=1:obj.nSlack
                obj.Vm(obj.idSlack(k)) = abs(obj.vBase(obj.idSlack(k)));
                obj.Va(obj.idSlack(k)) = angle(obj.vBase(obj.idSlack(k)));
            end
            % Set PV bus voltages in Vm and Va vector 
            for k=1:obj.nPV
                obj.Vm(obj.idPV(k)) = abs(obj.vBase(obj.idPV(k)));
                obj.Va(obj.idPV(k)) = angle(obj.vBase(obj.idSlack(1)));
            end
            % set PQ voltages in Vm and Va vector (set to vBase)
            for k=1:obj.nPQ
                obj.Vm(obj.idPQ(k)) = abs(obj.vBase(obj.idPQ(k)));
                obj.Va(obj.idPQ(k)) = angle(obj.vBase(obj.idPQ(k)));
            end
        end
        
        
        % Solve load flow equations
        function [vOut,pSlack,qSlack,qOut] = nrLoadFlow(obj,Y,Pin,Qin)
            % Initialization
            % Error and iterations
            err = 1;
            ite = 0;

            % Admitance matrix in polar coordinates
            Ym = abs(Y);                         % Magnitude of Ybus 
            Ya = angle(Y);                       % Angle of Ybus

            %% Run the Newton-Raphson Method
            while err >= obj.tol && ite < obj.maxIte

                % Iterate
                ite = ite+1;

                % Allocate Memory
                P = zeros(obj.nBus,1);
                Q = zeros(obj.nBus,1);
                J1 = zeros(obj.nBus-obj.nSlack,obj.nBus-obj.nSlack);
                J2 = zeros(obj.nBus-obj.nSlack,obj.nPQ);
                J3 = zeros(obj.nPQ,obj.nBus-obj.nSlack);
                J4 = zeros(obj.nPQ,obj.nPQ);

                % For load busses (PQ), calculate P, Q.                                (STEP 2)
                for k=1:obj.nPQ
                    for n=1:obj.nBus
                        P(obj.idPQ(k)) = P(obj.idPQ(k)) + obj.Vm(obj.idPQ(k))*obj.Vm(n)*Ym(obj.idPQ(k),n) * ...
                                        cos(Ya(obj.idPQ(k),n)-obj.Va(obj.idPQ(k))+obj.Va(n));
                        Q(obj.idPQ(k)) = Q(obj.idPQ(k)) - (obj.Vm(obj.idPQ(k))*obj.Vm(n)*Ym(obj.idPQ(k),n) * ...
                                        sin(Ya(obj.idPQ(k),n)-obj.Va(obj.idPQ(k))+obj.Va(n)));
                    end
                end

                % For voltage-controled busses (PV), calculate P.                      (STEP 3)
                for k=1:obj.nPV
                    for n=1:obj.nBus
                        P(obj.idPV(k)) = P(obj.idPV(k)) + obj.Vm(obj.idPV(k))*obj.Vm(n)*Ym(obj.idPV(k),n) * ...
                                        cos(Ya(obj.idPV(k),n)-obj.Va(obj.idPV(k))+obj.Va(n));
                    end
                end

                % Calculate the elements J1,J2,J3 and J4 of the Jacobian matrix        (STEP 4)
                % J1 has dim(nBus-nSlack,nBus-nSlack)
                for k=1:obj.nBus-obj.nSlack
                    % Diagonal elements of J1
                    for n=1:obj.nBus
                        if obj.idPQV(k)==n
                        else
                            J1(k,k) = J1(k,k) + obj.Vm(obj.idPQV(k))*obj.Vm(n)*Ym(obj.idPQV(k),n) * ...
                                        sin(Ya(obj.idPQV(k),n)-obj.Va(obj.idPQV(k))+obj.Va(n));

                        end
                    end
                    % Off diagonal elements of J1
                    for n=1:obj.nBus-obj.nSlack
                       if k==n
                       else
                            J1(k,n) = -(obj.Vm(obj.idPQV(k))*obj.Vm(obj.idPQV(n))*Ym(obj.idPQV(k),obj.idPQV(n)) * ...
                                        sin(Ya(obj.idPQV(k),obj.idPQV(n))-obj.Va(obj.idPQV(k))+obj.Va(obj.idPQV(n))));
                       end
                   end
                end

                % J2 has dim(nBus-nSlack,nPQ)
                for k=1:obj.nPQ
                    % Diagonal elements of J2
                    for n=1:obj.nBus
                        if obj.idPQV(k)==n
                        else
                            J2(k,k) = J2(k,k) + obj.Vm(n)*Ym(obj.idPQV(k),n) * ...
                                         cos(Ya(obj.idPQV(k),n)-obj.Va(obj.idPQV(k))+obj.Va(n));
                        end
                    end
                    % Add last element to diagonal of J2
                    J2(k,k) = J2(k,k) + 2*obj.Vm(obj.idPQV(k))*Ym(obj.idPQV(k),obj.idPQV(k)) * ... 
                                cos(Ya(obj.idPQV(k),obj.idPQV(k)));
                    % Off diagonal elements of J2
                    for n=1:obj.nBus-obj.nSlack
                        if k==n
                        else
                           J2(n,k) = obj.Vm(obj.idPQV(n))*Ym(obj.idPQV(n),obj.idPQV(k)) * ...
                                        cos(Ya(obj.idPQV(n),obj.idPQV(k))-obj.Va(obj.idPQV(n))+obj.Va(obj.idPQV(k))); 
                        end
                    end
                end

                % J3 has dim(nPQ,nBus-nSlack)
                for k=1:obj.nPQ
                    % Diagonal elements of J3
                    for n=1:obj.nBus
                        if obj.idPQV(k)==n
                        else
                            J3(k,k) = J3(k,k) + obj.Vm(obj.idPQV(k))*obj.Vm(n)*Ym(obj.idPQV(k),n) * ...
                                        cos(Ya(obj.idPQV(k),n)-obj.Va(obj.idPQV(k))+obj.Va(n));
                        end 
                    end
                    % Off diagonal elements of J3
                    for n=1:obj.nBus-obj.nSlack
                        if k==n
                        else
                            J3(k,n) = -(obj.Vm(obj.idPQV(k))*obj.Vm(obj.idPQV(n))*Ym(obj.idPQV(k),obj.idPQV(n)) * ...
                                        cos(Ya(obj.idPQV(k),obj.idPQV(n))-obj.Va(obj.idPQV(k))+obj.Va(obj.idPQV(n))));
                        end
                    end
                end

                % J4 has dim(nPQ,nPQ)
                for k=1:obj.nPQ
                    % Diagonal elements of J4
                    for n=1:obj.nBus
                        if obj.idPQV(k)==n
                        else
                            J4(k,k) = J4(k,k) - (obj.Vm(n)*Ym(obj.idPQV(k),n) * ...
                                        sin(Ya(obj.idPQV(k),n)-obj.Va(obj.idPQV(k))+obj.Va(n)));
                        end
                    end
                    % Add last element to diagonal of J4
                    J4(k,k) = J4(k,k) - 2*obj.Vm(obj.idPQV(k))*Ym(obj.idPQV(k),obj.idPQV(k)) * ...
                                sin(Ya(obj.idPQV(k),obj.idPQV(k)));
                    % Off diagonal elements of J4
                    for n=1:obj.nPQ
                        if k==n
                        else
                            J4(k,n) = -(obj.Vm(obj.idPQV(k))*Ym(obj.idPQV(k),obj.idPQV(n)) * ...
                                        sin(Ya(obj.idPQV(k),obj.idPQV(n))-obj.Va(obj.idPQV(k))+obj.Va(obj.idPQV(n))));
                        end
                    end
                end

                % Form J
                J = [J1 J2; J3 J4];

                % Find residuals and solve linear simultaneous equations               (STEP 5)
                dP = Pin(obj.idPQV) - P(obj.idPQV);
                dQ = Qin(obj.idPQ) - Q(obj.idPQ);
                dPQ = [dP; dQ];

                dVam = J\dPQ;

                % Compute new voltage magnitudes and phase angles                      (STEP 6)
                obj.Va(obj.idPQV) = obj.Va(obj.idPQV) + dVam(1:length(obj.Va(obj.idPQV)));
                obj.Vm(obj.idPQ) = obj.Vm(obj.idPQ) + dVam(length(obj.Va(obj.idPQV))+1:end);

                % Calculate error                                                      (STEP 7)
                err = max(abs(dPQ));

            end
            % Indicate number of iterations
            obj.nIte = ite;

            % Set outputs
            % Voltages at each bus including slack bus
            vOut = obj.Vm.*exp(1i*obj.Va);

            % Active and reactive power at slack bus
            Pslack = zeros(obj.nSlack,1);
            Qslack = zeros(obj.nSlack,1);
            for k=1:obj.nSlack
                for n=1:obj.nBus
                   Pslack(k) = Pslack(k) + obj.Vm(k)*obj.Vm(n)*Ym(k,n) * ...
                                cos(Ya(k,n)-obj.Va(k)+obj.Va(n));
                   Qslack(k) = Qslack(k) - (obj.Vm(k)*obj.Vm(n)*Ym(k,n) * ...
                                sin(Ya(k,n)-obj.Va(k)+obj.Va(n))); 
                end

            end

            % Reactive power at PV busses
            Qout = zeros(obj.nPV,1);
            for k=1:obj.nPV
                for n=1:obj.nBus
                    Qout(k) = Qout(k) - (obj.Vm(obj.idPV(k))*obj.Vm(n)*Ym(obj.idPV(k),n) * ...
                                sin(Ya(obj.idPV(k),n)-obj.Va(obj.idPV(k))+obj.Va(n))); 
                end
            end

            % Form output
            pSlack = Pslack;
            qSlack = Qslack;
            qOut = Qout;
        end
        
        % Transform grid line data into admittance matrix
        function Y = lDataToY(obj,zData)
            % Unpack data and save in object
            from = zData(:,1);                  % From bus
            to = zData(:,2);                    % To bus
            R = zData(:,3);                     % real(z)
            X = zData(:,4);                     % img(z)
            l = zData(:,5);                     % Length of cable    

            nbr = length(zData(:,1));           % Number of branches
            nbus = max(max(from),max(to));      % Number of busses

            % Branch impedance 
            Z = complex(R,X);
            Z = Z.*l;
            % Branch admitance
            y = ones(nbr,1)./Z;
            % Admitance Matrix
            Y = zeros(nbus,nbus);
            for k = 1:nbr                       % Off diagonal elements
                if from(k)>0 && to(k)>0
                    Y(from(k),to(k)) = Y(from(k),to(k)) - y(k);
                    Y(to(k),from(k)) = Y(from(k),to(k));
                end

            end
            for n = 1:nbus                      % Diagonal elements
                for k = 1:nbr
                    if from(k)==n || to(k)==n
                        Y(n,n) = Y(n,n) + y(k);
                    end
                end
            end
        end

        % Build the grid impedance matrix. Can be used for loss
        % minimization
        function [Zbus] = zbuild(obj,linedata)
            % This program forms the complex bus impedance matrix by the method
            % of building algorithm.  Bus zero is taken as reference.
            %  Copyright (c) 1998  by H. Saadat
            %

            % Update: 18-11-2014, by Rasmus Pedersen, Aalborg University for DiSC data
            % format. 
           
            % Unpack data
            from = linedata(:,1);               % From bus
            to = linedata(:,2);                 % To bus
            R = linedata(:,3);                  % Resistance of cable
            X = linedata(:,4);                  % Reactance of cable
            l = linedata(:,5);                  % Length of cable
            
            nbr = length(linedata(:,1));        % Number of branches
            nbus = max(max(from), max(to));     % Number of busses

            % Branch impedances
            ZB = complex(R,X).*l;
            % Allocate memory
            Zbus = zeros(nbus, nbus);
            ntree = zeros(nbr,1);
            tree = 0;  %%%%new
         
            % Adding a branch from a new bus to reference bus 0
            for I = 1:nbr
                ntree(I) = 1;
                if from(I) == 0 || to(I) == 0
                    if from(I) == 0
                        n = to(I);
                    elseif to(I) == 0
                        n = from(I);
                    end
                    if abs(Zbus(n, n)) == 0
                          Zbus(n,n) = ZB(I);
                          tree = tree + 1; %%new
                    else 
                        Zbus(n,n) = Zbus(n,n)*ZB(I)/(Zbus(n,n) + ZB(I));
                    end
                    ntree(I) = 2;
                end
            end
            % Adding a branch from new bus to an existing bus
            while tree < nbus  %%% new
                for n = 1:nbus
                    nadd = 1;
                    if abs(Zbus(n,n)) == 0
                        for I = 1:nbr
                            if nadd == 1;
                                if from(I) == n || to(I) == n
                                    if from(I) == n
                                        k = to(I);
                                    elseif to(I) == n
                                        k = from(I);
                                    end
                                    for m = 1:nbus
                                        if m ~= n
                                            Zbus(m,n) = Zbus(m,k);
                                            Zbus(n,m) = Zbus(m,k);
                                        end
                                    end
                                    Zbus(n,n) = Zbus(k,k) + ZB(I); tree=tree+1; %%new
                                    nadd = 2;  ntree(I) = 2;
                                end
                            end
                        end
                    end
                end
            end  %%%%%%new
            % Adding a link between two old buses
            for n = 1:nbus
                for I = 1:nbr
                    if ntree(I) == 1
                        if from(I) == n || to(I) == n
                            if from(I) == n
                                k = to(I);
                            elseif to(I) == n 
                                k = from(I);
                            end
                            DM = Zbus(n,n) + Zbus(k,k) + ZB(I) - 2*Zbus(n,k);
                            for jj = 1:nbus
                                AP = Zbus(jj,n) - Zbus(jj,k);
                                for kk = 1:nbus
                                    AT = Zbus(n,kk) - Zbus(k, kk);
                                    DELZ(jj,kk) = AP*AT/DM;
                                end
                            end
                            Zbus = Zbus - DELZ;
                            ntree(I) = 2;
                        end
                    end
                end
            end
        end % End function
        
        % Calculate line loading and loss
        function output=lineLoading(obj,zData,trafoPos,trafoZ,trafoRatio,v)
            % Unpack data and save in object
            from = zData(:,1);                  % From bus
            to = zData(:,2);                    % To bus
            R = zData(:,3);                     % real(z)
            X = zData(:,4);                     % img(z)
            l = zData(:,5);                     % Length of cable    

            nbr = length(zData(:,1));           % Number of branches
            numBus = max(max(from),max(to));      % Number of busses
            
            % Branch impedance 
            Z = complex(R,X);
            Z = Z.*l;
            % Insert trafo impedance
            Z(trafoPos) = trafoZ;
            % Branch admitance
            y = ones(nbr,1)./Z;
            % Set vector of transformer ratios
            tR = ones(nbr,1);
            tR(trafoPos) = trafoRatio;
            output.from = from;
            output.to = to;
            
            for n=1:numBus
                for j=1:nbr
                    if from(j)==n
                        k = to(j);
                        In = (v(n)-v(k)*tR(j))*(y(j)/tR(j)^2);
                        Ik = (v(k)-v(n)/tR(j))*y(j);
                        Snk = v(n)*conj(In);
                        Skn = v(k)*conj(Ik);
                        SL = Snk+Skn;
                        output.SfromTo(j) = Snk;
                        output.In(j) = In;
                        output.Ik(j) = Ik; 
                        output.Sloss(j) = SL;
                    elseif to(j)==n
                        k = from(j);
                        In = (v(n)-v(k)/tR(j))*y(j);
                        Ik = (v(k)-v(n)*tR(j))*(y(j)/tR(j)^2);
                        Snk = v(n)*conj(In);
                        Skn = v(k)*conj(Ik);
                        SL = Snk+Skn;
                        output.SfromTo(j) = Snk;
                        output.In(j) = In;
                        output.Ik(j) = Ik; 
                        output.Sloss(j) = SL;
                    end
                end
            end
        end
        
        % Compute Jacobian of power flow equations
        function [J,J1,J2,J3,J4] = gridJacobian(obj,v,Y)
            % Input:
            %   - v [V], is a vector of complex valued voltages where power flow equations are linearized around.
            %   - Y [S], is the grid admittance matrix. NOTE: remember
            %   that transformers must be placed before calling this
            %   function. Otherwise inf numbers are present.
            %
            % Output:
            %   - J [-], is the full grid Jacobian matrix.
            %   - J1 [-], is the upper left part of the Jacobian, see H.
            %   Saadat "Power System Analysis" page 233.
            %   - J2 [-], is the upper right part of the Jacobian.
            %   - J3 [-], is the lower left part of the Jacobian.
            %   - J4 [-], is the lower right part of the Jacobian.
            
            % Get input values on polar form:
            Ym = abs(Y);                         % Magnitude of Y 
            Ya = angle(Y);                       % Angle of Y 
            vm = abs(v);                         % Magnitude of v    
            va = angle(v);                       % Angle of v
            
            % Allocate memory
            J1 = zeros(obj.nBus-obj.nSlack,obj.nBus-obj.nSlack);
            J2 = zeros(obj.nBus-obj.nSlack,obj.nPQ);
            J3 = zeros(obj.nPQ,obj.nBus-obj.nSlack);
            J4 = zeros(obj.nPQ,obj.nPQ);
            
            % Calculate the elements J1,J2,J3 and J4 of the Jacobian matrix
            % J1 has dim(nBus-nSlack,nBus-nSlack)
            for k=1:obj.nBus-obj.nSlack
                % Diagonal elements of J1
                for n=1:obj.nBus
                    if obj.idPQV(k)==n
                    else
                        J1(k,k) = J1(k,k) + vm(obj.idPQV(k))*vm(n)*Ym(obj.idPQV(k),n) * ...
                                    sin(Ya(obj.idPQV(k),n)-va(obj.idPQV(k))+va(n));

                    end
                end
                % Off diagonal elements of J1
                for n=1:obj.nBus-obj.nSlack
                   if k==n
                   else
                        J1(k,n) = -(vm(obj.idPQV(k))*vm(obj.idPQV(n))*Ym(obj.idPQV(k),obj.idPQV(n)) * ...
                                    sin(Ya(obj.idPQV(k),obj.idPQV(n))-va(obj.idPQV(k))+va(obj.idPQV(n))));
                   end
               end
            end

            % J2 has dim(nBus-nSlack,nPQ)
            for k=1:obj.nPQ
                % Diagonal elements of J2
                for n=1:obj.nBus
                    if obj.idPQV(k)==n
                    else
                        J2(k,k) = J2(k,k) + vm(n)*Ym(obj.idPQV(k),n) * ...
                                     cos(Ya(obj.idPQV(k),n)-va(obj.idPQV(k))+va(n));
                    end
                end
                % Add last element to diagonal of J2
                J2(k,k) = J2(k,k) + 2*vm(obj.idPQV(k))*Ym(obj.idPQV(k),obj.idPQV(k)) * ... 
                            cos(Ya(obj.idPQV(k),obj.idPQV(k)));
                % Off diagonal elements of J2
                for n=1:obj.nBus-obj.nSlack
                    if k==n
                    else
                       J2(n,k) = vm(obj.idPQV(n))*Ym(obj.idPQV(n),obj.idPQV(k)) * ...
                                    cos(Ya(obj.idPQV(n),obj.idPQV(k))-va(obj.idPQV(n))+va(obj.idPQV(k))); 
                    end
                end
            end

            % J3 has dim(nPQ,nBus-nSlack)
            for k=1:obj.nPQ
                % Diagonal elements of J3
                for n=1:obj.nBus
                    if obj.idPQV(k)==n
                    else
                        J3(k,k) = J3(k,k) + vm(obj.idPQV(k))*vm(n)*Ym(obj.idPQV(k),n) * ...
                                    cos(Ya(obj.idPQV(k),n)-va(obj.idPQV(k))+va(n));
                    end 
                end
                % Off diagonal elements of J3
                for n=1:obj.nBus-obj.nSlack
                    if k==n
                    else
                        J3(k,n) = -(vm(obj.idPQV(k))*vm(obj.idPQV(n))*Ym(obj.idPQV(k),obj.idPQV(n)) * ...
                                    cos(Ya(obj.idPQV(k),obj.idPQV(n))-va(obj.idPQV(k))+va(obj.idPQV(n))));
                    end
                end
            end

            % J4 has dim(nPQ,nPQ)
            for k=1:obj.nPQ
                % Diagonal elements of J4
                for n=1:obj.nBus
                    if obj.idPQV(k)==n
                    else
                        J4(k,k) = J4(k,k) - (vm(n)*Ym(obj.idPQV(k),n) * ...
                                    sin(Ya(obj.idPQV(k),n)-va(obj.idPQV(k))+va(n)));
                    end
                end
                % Add last element to diagonal of J4
                J4(k,k) = J4(k,k) - 2*vm(obj.idPQV(k))*Ym(obj.idPQV(k),obj.idPQV(k)) * ...
                            sin(Ya(obj.idPQV(k),obj.idPQV(k)));
                % Off diagonal elements of J4
                for n=1:obj.nPQ
                    if k==n
                    else
                        J4(k,n) = -(vm(obj.idPQV(k))*Ym(obj.idPQV(k),obj.idPQV(n)) * ...
                                    sin(Ya(obj.idPQV(k),obj.idPQV(n))-va(obj.idPQV(k))+va(obj.idPQV(n))));
                    end
                end
            end

            % Form J
            J = [J1 J2; J3 J4]; 
        end

        
        %% Set and get functions
        function setParam(obj,maxIte,tol)
            obj.maxIte = maxIte;
            obj.tol = tol;
        end
    end
end
