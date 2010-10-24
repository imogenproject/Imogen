classdef BonnerEbertInitializer < Initializer
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
    rho0;           % Gas density at the cloud's center.                            ???
    sphereGAMMA;    % ???                                                           ???
    sphereK;        % ???                                                           ???
    sphereRmax;     % ???                                                           ???
    edgePadding     % Determines how much space to keep at the edge of the          ???
                    %   grid at background density.
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]

    BgDensity;
    BgDensityCoeff;

    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]

    pBgDensityCoeff;
    

    end %PROTECTED

%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ BonnerEbertInitializer
        function obj = BonnerEbertInitializer(input)
            obj = obj@Initializer();
            obj.runCode             = 'BON_EBRT';
            obj.info                = 'Bonner-Ebert hydrostatic gas sphere simulation';
            obj.mode.fluid          = true;
            obj.mode.magnet         = false;
            obj.mode.gravity        = true;
            obj.iterMax             = 300;
            obj.bcMode.x            = 'fade';
            obj.bcMode.y            = 'fade';
            obj.bcMode.z            = 'fade';
            obj.bcInfinity          = 5;

            obj.activeSlices.xy     = true;
            obj.timeUpdateMode      = ENUM.TIMEUPDATE_PER_STEP;

            obj.pBgDensityCoeff     = 1e-3;
            
            obj.gravity.solver      = 'biconj';

            obj.gamma               = 2; % Imogen gamma

            obj.rho0                = .1;
            obj.sphereGAMMA         = 1;
            obj.sphereK             = 1;
            obj.sphereRmax          = 15;

            obj.dGrid               = [1 1 1];
            obj.thresholdMass       = obj.rho0 * obj.pBgDensityCoeff * 2;
            
            obj.operateOnInput(input, [64 64 64]);
        end
        
%___________________________________________________________________________________________________ GS: BgDensity
        function value = get.BgDensity(obj)
            value = obj.rho0*obj.pBgDensityCoeff;
        end
        %\ @private \%
        function set.BgDensity(obj, value)
            obj.pBgDensityCoeff = value/obj.rho0;
            obj.thresholdMass = obj.rho0 * obj.pBgDensityCoeff * 2;
        end

%___________________________________________________________________________________________________ GS: BgDensityCoeff
        function value = get.BgDensityCoeff(obj)
            value = obj.pBgDensityCoeff;
        end
        %\ @private \%
        function set.BgDensityCoeff(obj, value)
            obj.pBgDensityCoeff = value;
            obj.thresholdMass = obj.rho0 * obj.pBgDensityCoeff;
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
            [R RHO E]       = computeBonnerEbert(obj.rho0, obj.sphereGAMMA, ...
                                            0.5*obj.sphereRmax / min(obj.grid), ...
                                            obj.sphereRmax, obj.sphereK, ...
                                            obj.rho0 * obj.pBgDensityCoeff);

            [X Y Z]         = ndgrid(1:obj.grid(1), 1:obj.grid(2), 1:obj.grid(3));

            obj.dGrid       = 2.4*max(R) ./ obj.grid;
            X               = (X - obj.grid(1)/2) * obj.dGrid(1);
            Y               = (Y - obj.grid(2)/2) * obj.dGrid(2);
            Z               = (Z - obj.grid(3)/2) * obj.dGrid(3);
            sphereRadii     = sqrt(X.^2 + Y.^2 + Z.^2);
            
            mass            = pchip(R, RHO, sphereRadii);
            mom             = zeros([3 obj.grid]);
            ener            = pchip(R, E/(obj.gamma - 1), sphereRadii);
            mag             = zeros([3 obj.grid]);

            obj.minMass     = obj.rho0 * obj.pBgDensityCoeff;
            
            ener            = ener + ... % internal already computed by initializer
                                + 0.5*squeeze(sum(mom .* mom, 1)) ./ mass ...   % kinetic energy
                                + 0.5*squeeze(sum(mag .* mag, 1));              % magnetic energy                    
                        
            statics = [];
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
