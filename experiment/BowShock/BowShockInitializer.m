classdef BowShockInitializer < Initializer
% Creates initial conditions for a bow shock simulation. 
%
% Unique properties for this initializer:
%   stencil      % File name for the statics stencil (must be in data dir).                 str
%   staticType   % Enumerated specification of how to apply static values.                  str
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
        PRIMAY_MODE = 'primary'; % Statics are applied to array classes.
        FLUX_MODE   = 'flux';    % Statics are applied to fluxes.
        FLUX_LR_MODE = 'fluxlr'; % Statics are applied to the left and right TVD fluxes only. 
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        stencil;      % File name for the statics stencil (must be in data dir).            str
        staticType;   % Enumerated specification of how to apply static values.             str

        ballCells;    % 3x1, radii in all 3 dimensions                                      double
        ballCenter;   % 3x1, center of the ball                                             double

        magneticShock;
	magX;
	magY;

        bgRho;             % Background density, x velocity and thermal pressure; double scalar
        bgVx;
        bgThermalPressure;

        ballRho;           % Obstruction density, outflow velocity, radius and thermal pressure
        ballVr;            % Double scalar
        ballXRadius;
        ballThermalPressure;

    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ BowShockInitializer
        function obj = BowShockInitializer(input)           
            obj = obj@Initializer();
            obj.gamma            = 5/3;
            obj.runCode          = 'Bow';
            obj.info             = 'Bow shock trial.';
            obj.mode.fluid       = true;
            obj.mode.magnet      = false;
            obj.mode.gravity     = false;
            obj.cfl              = 0.7;
            obj.iterMax          = 100;
            obj.bcMode.x         = 'trans';
            obj.bcMode.y         = 'fade';
            if input(3) > 1
                obj.bcMode.z     = 'fade';
            else
                obj.bcMode.z     = 'circ';
            end
            obj.activeSlices.xy  = true;
            obj.ppSave.dim2      = 25;
            
            obj.staticType       = BowShockInitializer.PRIMAY_MODE;
            obj.stencil          = 'SmallSphere_800x256.mat';
           
            obj.ballCells        = [32 32 32];
            obj.ballCenter       = round(input/2);
	    obj.ballXRadius      = 1;

            obj.magX                = 0;
            obj.magY                = 0;

            obj.bgRho               = .125;
            obj.bgVx                = 1;
            obj.bgThermalPressure   = .03;

            obj.ballRho             = 1;
            obj.ballVr              = 1;
            obj.ballThermalPressure = .25;
         
            obj.operateOnInput(input, [800, 256, 1]);
        end
               
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = protected) %                                          P R O T E C T E D    [M]
        
%___________________________________________________________________________________________________ calculateInitialConditions
        function [mass, mom, ener, mag, statics] = calculateInitialConditions(obj)
        % Returns the initial conditions for a bow shock simulation
        % USAGE: [mass, mom, ener, mag, statics, run] = getInitialConditions();
 
            %--- Background Values ---%
            mass        = zeros(obj.grid);
            mom         = zeros([3 obj.grid]);
            mag         = zeros([3 obj.grid]);
            ener        = zeros(obj.grid);

            %--- Static Values ---%
            statics = StaticsInitializer();

            [X Y Z] = ndgrid(1:obj.grid(1), 1:obj.grid(2), 1:obj.grid(3));
            Ledge = (X < 8); % Left edge - we establish plane flow here

            % The obstacle is a spheroid 
            X = X - obj.ballCenter(1);
            Y = Y - obj.ballCenter(2);
            Z = Z - obj.ballCenter(3);
            norm = sqrt((X/obj.ballCells(1)).^2 + (Y/obj.ballCells(2)).^2 + (Z/obj.ballCells(3)).^2);
            ball = (norm <= 1.0);

            obj.minMass = min(obj.ballRho, obj.bgRho) / 100;

            % set background values for momentum
            fastedge = max(round(obj.ballCenter(1) - obj.ballCells(1)-100), 16);
            mom(1, 1:fastedge,:,:) = obj.bgVx * obj.bgRho;
            mom(1, fastedge:end,1:round(obj.ballCenter(2)-obj.ballCells(1)-20),:) =  obj.bgVx * obj.bgRho;
            mom(1, fastedge:end,round(obj.ballCenter(2)+obj.ballCells(1)+20):end,:) =  obj.bgVx * obj.bgRho;
            
            mass(:) = obj.bgRho;
            mass(ball) = obj.ballRho;

            ener(:) = mass.^obj.gamma / (obj.gamma-1) + .5*squeeze(mom(1,:,:,:).^2)./mass;

	    obj.dGrid = obj.ballXRadius / obj.ballCells(1);

            statics.indexSet{1} = indexSet_fromLogical(ball); % ball
            statics.indexSet{2} = indexSet(obj.grid, 1:2, 1:(obj.grid(2)-0), 1:obj.grid(3)); % incoming blast
            statics.indexSet{3} = indexSet(obj.grid, (obj.grid(1)-4):(obj.grid(1)-0), 1:obj.grid(2), 1:obj.grid(3)); % right edge
            statics.indexSet{4} = indexSet(obj.grid, 1:(obj.grid(1)-2), 1:2, 1:obj.grid(3)); % top edge

            xhat = X/obj.ballCells(1);
            yhat = Y/obj.ballCells(2);
            zhat = Z/obj.ballCells(3); 
            ballMomRadial = obj.ballVr*obj.ballRho;
            ballEner      = (obj.ballRho^obj.gamma)/(obj.gamma-1) + .5*ballMomRadial^2.*norm(ball)/obj.ballRho;

            statics.valueSet = {0, obj.bgRho, obj.bgRho*obj.bgVx, ener(1), ...
                obj.ballRho, ballMomRadial*xhat(ball), ballMomRadial*yhat(ball), ballMomRadial*zhat(ball), ballEner, obj.magX, obj.magY, obj.bgRho.^obj.gamma / (obj.gamma-1) };

            % Force a left-edge plane flow
            statics.setFluid_allConstantBC(mass, ener, mom, 1);

        
            % Lock ball in place
            statics.associateStatics(ENUM.MASS, ENUM.SCALAR,    statics.CELLVAR, 1, 5);
            statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(1), statics.CELLVAR, 1, 6);
            statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(2), statics.CELLVAR, 1, 7);
%            if obj.mode.magnet == true;
%                statics.associateStatics(ENUM.MAG,  ENUM.VECTOR(1), statics.CELLVAR, 1, 10);
%                statics.associateStatics(ENUM.MAG,  ENUM.VECTOR(2), statics.CELLVAR, 1, 11);
%            end
            if obj.grid(3) > 1    
                statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(3), statics.CELLVAR, 1, 8);
%                if obj.mode.magnet == true; statics.associateStatics(ENUM.MAG,  ENUM.VECTOR(3), statics.CELLVAR, 1, 1); end;
            end
            statics.associateStatics(ENUM.ENER, ENUM.SCALAR,    statics.CELLVAR, 1, 9);

            % Force constant on right edge (crude approximation of fade)
            %statics.associateStatics(ENUM.MASS, ENUM.SCALAR,    statics.CELLVAR, 3, 2);
            %statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(1), statics.CELLVAR, 3, 1);
            %statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(2), statics.CELLVAR, 3, 1);
%            if obj.mode.magnet == true;
%                statics.associateStatics(ENUM.MAG,  ENUM.VECTOR(1), statics.CELLVAR, 2, 10);
%                statics.associateStatics(ENUM.MAG,  ENUM.VECTOR(2), statics.CELLVAR, 2, 11);
%            end
            %if obj.grid(3) > 1
            %    statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(3), statics.CELLVAR, 3, 1);
%                if obj.mode.magnet == true; statics.associateStatics(ENUM.MAG,  ENUM.VECTOR(3), statics.CELLVAR, 2, 1); end
            %end
            %statics.associateStatics(ENUM.ENER, ENUM.SCALAR,    statics.CELLVAR, 3, 12);

            % Force constant on top edge
            statics.setFluid_allConstantBC(mass, ener, mom, 4);
%            statics.associateStatics(ENUM.MASS, ENUM.SCALAR,    statics.CELLVAR, 4, 2);
%            statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(1), statics.CELLVAR, 4, 3);
%            statics.associateStatics(ENUM.MOM,  ENUM.VECTOR(2), statics.CELLVAR, 4, 1);
%            statics.associateStatics(ENUM.ENER, ENUM.SCALAR,    statics.CELLVAR, 4, 4);
            
        end

    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
