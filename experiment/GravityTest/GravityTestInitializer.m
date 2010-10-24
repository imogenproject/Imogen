classdef GravityTestInitializer < Initializer
% Creates a Bonner-Ebert self-gravitating hydrostatic gas sphere
%
% Unique properties for this initializer:
% rho0        Gas density at the center of the cloud
%         -- By extension, run.gravity.G will be important

%   useStatics        specifies if static conditions should be set for the run.       logical
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
    rhoSphere;     % Gas density of spheres
    rhoBG;         % Gas density of background
    enerSphere;    % Define energy/temperature of spheres

    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]

    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]

    end %PROTECTED

%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ GravityTestInitializer
        function obj = GravityTestInitializer(input)            
            obj                     = obj@Initializer();
            obj.grid = input;
            obj.runCode             = 'GRAV_TOY';
            obj.info                = 'Self-gravitating toy problem';
            obj.mode.fluid          = true;
            obj.mode.magnet         = false;
            obj.mode.gravity        = true;
            obj.iterMax             = 300;
            obj.gamma               = 2;
            obj.bcMode.x            = 'fade';
            obj.bcMode.y            = 'fade';
            obj.bcMode.z            = 'fade';

            obj.bcInfinity          = 4;

            obj.activeSlices.xy     = true;
            obj.timeUpdateMode      = ENUM.TIMEUPDATE_PER_STEP;
            
            obj.gravity.constant    = 1;
            obj.gravity.solver      = ENUM.GRAV_SOLVER_BICONJ;

            obj.rhoSphere           = 1;
            obj.enerSphere          = 0.05;
            obj.rhoBG               = 0.05;
            obj.thresholdMass       = obj.rhoBG * 1.5;
            
            obj.operateOnInput(input, [200, 100, 100]);
        end
        
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
            for i=1:3
                if (bitget(obj.grid(i),1) && obj.grid(i) > 1);  obj.grid(i) = obj.grid(i) + 1; end
            end
            
            % Numerically compute the hydrostatic density balance
            % Fine resolution appears important - 
            [mass, ener] = generateTestSystem(obj.grid, obj.rhoSphere, obj.rhoBG, ...
                                                            obj.gravity.constant, obj.enerSphere);
            
            obj.dGrid   = [2 2 2] / obj.grid(1);
            mom         = zeros([3 obj.grid]);
            obj.minMass = .5*obj.thresholdMass;
            mag         = zeros([3 obj.grid]);
            statics     = [];
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
