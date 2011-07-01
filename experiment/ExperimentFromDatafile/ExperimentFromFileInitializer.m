classdef ExperimentFromFileInitializer < Initializer
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
        bgDensityCoeff;

        fname;
        latheAboutZ;
        
        momFormat;
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
        function obj = ExperimentFromFileInitializer(input)
            obj                     = obj@Initializer();            
            obj.runCode             = 'LD_FROM_FILE';
            obj.info                = 'running data loaded from file';
            obj.mode.fluid          = true;
            obj.mode.magnet         = false;
            obj.mode.gravity        = true;
            obj.iterMax             = 300;
            obj.bcMode.x            = ENUM.BCMODE_FADE;
            obj.bcMode.y            = ENUM.BCMODE_FADE;
            obj.bcMode.z            = ENUM.BCMODE_FADE;
            obj.bcInfinity          = 5;
            obj.activeSlices.xy     = true;
            obj.timeUpdateMode      = ENUM.TIMEUPDATE_PER_STEP;
            obj.bgDensityCoeff      = 1e-4;
            
            obj.gravity.constant        = 1;
            obj.gravity.solver      = ENUM.GRAV_SOLVER_BICONJ;
            obj.gravity.tolerance = 1e-10;


            obj.gamma                   = 5/3;

            obj.thresholdMass           = 0;
            obj.grid = input;
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

            load(obj.fname);
            isangmom = 0;

            % Only cases 1 and 4 make sense for non-axial applications.
            switch obj.momFormat
                case 1; % linear momentum density, no operation needed
                case 2; isangmom = 1; % angular momentum density, divide by r
                case 3; data.mom = data.mom .* data.mass; isangmom = 1; % angular velocity, mult. by rho / r.
                case 4; data.mom = data.mom .* data.mass; % linear velocity, mult by rho
            end
            
            mass    = zeros(obj.grid);
            mom     = zeros([3 obj.grid]);
            
            if obj.latheAboutZ
                [mass, mom(1,:,:,:), mom(2,:,:,:)]  = latheMassMom(data.mass, data.mom, obj.grid, data.h, isangmom);
            else
                mass = data.mass;
                mom = data.mom;
            end

            obj.dGrid = [data.h data.h data.h];
            obj.appendInfo(sprintf('Automatically set dGrid uniformly to %d', data.h));
           
            tempd = cell(1,3);
            tempd{1} = obj.dGrid(1); tempd{2} = obj.dGrid(2); tempd{3} = obj.dGrid(3);
 
            obj.minMass = maxFinderND(mass) * obj.bgDensityCoeff;
            obj.thresholdMass = 10*obj.minMass;

            mass    = max(mass, obj.minMass);
            mag     = zeros([3 obj.grid]);
                        
            ener    = (mass.^obj.gamma)/(obj.gamma - 1) ...   % internal energy
                        + 0.5*squeeze(sum(mom .* mom, 1)) ./ mass ...           % kinetic energy
                        + 0.5*squeeze(sum(mag .* mag, 1));                      % magnetic energy                    

            statics = [];
            
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
