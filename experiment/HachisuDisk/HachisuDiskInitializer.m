classdef HachisuDiskInitializer < Initializer
% Uses equilibrium routines to create an equilibrium disk configuration and then formats them for
% use in Imogen based on the Kojima model. The disk is created in the polytropic unit systems under
% the constraints that the gravitational constant, G, the mass of the central star, M, and the 
% polytropic constant, K, are all unity.
%
% Unique properties for this initializer:
%   q                 angular velocity exponent (omega = omega_0 * r^-q.              double
%   radiusRatio       (inner disk radius) / (radius of density max).                  double
%   edgePadding       number of cells around X & Y edges to leave blank.              double
%   pointRadius       size of the softened point at the center of the grid            double
%   useStatics        specifies if static conditions should be set for the run.       logical
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]

    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        q;              % angular velocity exponent; omega^2 = omega_0^2(r^2 + d^2)^-q  double
        d;              % 
        aspectRatio;    % inner radius / radius of density max.                         double
        structType;     % either 'torus' or 'spheroid'                                  string
        edgePadding;    % distance to leave from outer disk edge to grid edge.          double

        pointRadius;    % size of the softened point at the center of the grid.         double
        starMass;       % Mass of the central object                                    double

        bgDensityCoeff; % Min density is this times max initial density                 double
        useUpperMirror; % If 1 simulates the top half of the disk only                  logical
        useRightMirror; % If 1 simulates right half of disk only                        logical

        useStatics;     % specifies if static conditions should be set for the run.     logical
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ KojimaDiskInitializer
        function obj = HachisuDiskInitializer(input)
            obj                   = obj@Initializer();            
            obj.runCode           = 'HACHISU';
            obj.info              = 'Self-gravitating rotating structure.';
            obj.mode.fluid        = true;
            obj.mode.magnet       = false;
            obj.mode.gravity      = true;
            obj.iterMax           = 300;
            obj.bcMode.x          = ENUM.BCMODE_FADE;
            obj.bcMode.y          = ENUM.BCMODE_FADE;
            obj.bcMode.z          = ENUM.BCMODE_FADE;
%            obj.bcMode.mom.z      = ENUM.BCMODE_FLIP;
%            obj.bcMode.mom.flux.z = ENUM.BCMODE_FLIP;
            
            obj.bcInfinity        = 5;
            obj.activeSlices.xy   = true;
            obj.timeUpdateMode    = ENUM.TIMEUPDATE_PER_STEP;
            obj.bgDensityCoeff    = 1e-4;
            
            obj.gravity.constant  = 1;
            obj.pointRadius       = 0.3;
            obj.gamma             = 5/3;

            obj.structType        = 'torus';
            obj.starMass         = 0;
            obj.q                 = 2;
            obj.d                 = .1;
            obj.aspectRatio       = 0.3;
            obj.edgePadding       = 0.1;

            obj.thresholdMass     = 2*obj.bgDensityCoeff;
            obj.useStatics        = false;

            obj.useUpperMirror    = 0; % Defaults to fully global simulation
            obj.useRightMirror    = 0;
            
            obj.operateOnInput(input, [64 64 32]);
            
        end
        
%___________________________________________________________________________________________________ GS: pointMass
% Dynmaic pointMass property to connect between the gravity vars structure and run files.
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = protected) %                                          P R O T E C T E D    [M]                
        
%___________________________________________________________________________________________________ calculateInitialConditions
        function [mass, mom, ener, mag, statics] = calculateInitialConditions(obj)

            %--- Ensure that the grid dimensions are even. ---%
            %       Even grid size means that the star will be correctly placed in the center cell.
            for i=1:2
                if (bitget(obj.grid(i),1) && obj.grid(i) > 1); obj.grid(i) = obj.grid(i) + 1; end
            end

            mom     = zeros([3 obj.grid]);
            mass    = zeros(obj.grid);

            % If using only half of global domain, turn gravity mirror on
            if obj.useUpperMirror == true
                diskres             = [obj.grid(1)/2 obj.grid(3)*2];
                obj.gravity.mirrorZ = true;
                obj.bcMode.z        = ENUM.BCMODE_CONST;
            end
            if obj.useRightMirror == true
                diskres             = [obj.grid(1)/2 obj.grid(2)];
                obj.gravity.mirrorZ = true;
                obj.bcMode.z        = ENUM.BCMODE_WORMHOLE;
%		obj.bcMode.flux.z = ENUM.BCMODE_CONST;
            end
            if (obj.useUpperMirror == false) && (obj.useRightMirror == false)
                diskres             = [obj.grid(1)/2 obj.grid(3)];
                obj.gravity.mirrorZ = false;
                obj.bcMode.z        = ENUM.BCMODE_FADE;
            end
            
            [rho, lMom, dGrid, R, info] = hachisuDisk(obj.aspectRatio, obj.structType, obj.q, obj.d,...
                                   1/(obj.gamma-1), diskres, obj.starMass, obj.edgePadding, 1e-4, 4*obj.bgDensityCoeff);

            if obj.useUpperMirror == true
                z0 = obj.grid(3);
                for zct = 1:obj.grid(3)
                    [mass(:,:,zct) mom(1,:,:,zct) mom(2,:,:,zct)] = cyl2rect(R(:,zct+z0), rho(:,zct+z0), lMom(:,zct+z0), obj.grid(1)/2, dGrid);
                end
            end
            if obj.useRightMirror == true
                for yct = 1:obj.grid(2)
                    [u v w] = cyl2rect(R(:,yct), rho(:,yct), lMom(:,yct), obj.grid(1)/2, dGrid);
%		size(mass(:,yct,:))
%		size(u(:,(end/2+1):end))
                    mass(:,yct,:)  = u(:,(end/2+1):end);
                    mom(1,:,yct,:) = v(:,(end/2+1):end);
                    mom(3,:,yct,:) = w(:,(end/2+1):end);
                end
            end
            if (obj.useUpperMirror == false) && (obj.useRightMirror == false)
                for zct = 1:obj.grid(3)
                    [mass(:,:,zct) mom(1,:,:,zct) mom(2,:,:,zct)] = cyl2rect(R(:,zct), rho(:,zct), lMom(:,zct), obj.grid(1)/2, dGrid);
                end
            end

            obj.appendInfo(sprintf('Automatically set dGrid uniformly to %d', dGrid));
            obj.dGrid = dGrid*ones(1,3);

            obj.thresholdMass = 2*obj.bgDensityCoeff * max(mass(:));
            obj.minMass = obj.bgDensityCoeff * max(mass(:));
                        
            mass    = max(mass, obj.minMass);

            % Purely hydro disks at the moment
            mag     = zeros([3 obj.grid]);
            
            ener    = (max(mass, obj.minMass).^obj.gamma)/(obj.gamma - 1) ...   % internal energy
                        + 0.5*squeeze(sum(mom .* mom, 1)) ./ mass ...           % kinetic energy
                        + 0.5*squeeze(sum(mag .* mag, 1));                      % magnetic energy                    
            
            if (obj.useStatics)
                statics.values(1) = 1e-5;
                statics.values(2) = 0;
                statics.values(3) = (obj.minMass .^ obj.gamma)/(obj.gamma - 1);

                index = zeros(2,3);
                index(1,:) = max( ceil(obj.grid/2)-2, 1);
                index(2,:) = min( floor(obj.grid/2)+2, obj.grid);
                indeces = cell(1,3);
                for i=1:3;  indeces{i} = index(1,i):index(2,i); end
                stencil = uint8(zeros(size(mass)));
                stencil(indeces{:}) = 1;

                statics.mass.s = stencil;
                statics.mom.s = uint8(zeros([3 size(stencil)]));
                for i=1:3; statics.mom.s(i,:,:,:) = 2*stencil; end
                statics.ener.s = 3*stencil;
            else statics = [];
            end

            if obj.starMass > 0
            
                if obj.useZMirror == 1
                    obj.gravity.fixedPotential = grav_GetPointPotential(obj.grid, tempd, ...
                    [obj.grid(1)/2 obj.grid(2)/2 0] + [.5 .5 0], obj.starMass, obj.pointRadius); % Temporary kludge
                else
                    obj.gravity.fixedPotential = grav_GetPointPotential(obj.grid, tempd, ...
                    obj.grid/2 + [.5 .5 .5], obj.starMass, obj.pointRadius); % Temporary kludge
                end
            end

        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
