classdef SedovTaylorBlastWaveInitializer < Initializer
% Creates initial conditions for a Sedov-Taylor blast wave simulation. This is a high energy point
% located centrally in the grid that shocks outward in a spherical manner. A good test of how well
% the Cartesian grid handles non-aligned propagation. This problem is essentially the simulation of
% an explosion, and is developed based on the original nuclear simulation work of Sedov and Taylor.
%
% Unique properties for this initializer:
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ SedovTaylorBlastWaveInitializer
        function obj = SedovTaylorBlastWaveInitializer(input)
            obj                  = obj@Initializer();
            obj.gamma            = 1.4;
            obj.runCode          = 'BLAST';
            obj.info             = 'Sedov-Taylor blast wave trial.';
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = false;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.7;
            obj.iterMax          = 100;
            obj.bcMode.x		 = ENUM.BCMODE_FADE;
            obj.bcMode.y		 = ENUM.BCMODE_FADE;
            obj.bcMode.z         = ENUM.BCMODE_FADE;
            obj.activeSlices.xy  = true;
            obj.activeSlices.xyz = true;
            obj.ppSave.dim2      = 5;
            obj.ppSave.dim3      = 20;
            
            obj.operateOnInput(input, [65, 65, 65]);
        end
               
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]
        
    end%PUBLIC
    
%===================================================================================================    
    methods (Access = protected) %                                          P R O T E C T E D    [M]

%___________________________________________________________________________________________________ calculateInitialConditions
        function [mass, mom, ener, mag, statics] = calculateInitialConditions(obj)
        
            %--- Initialization ---%
            statics         = [];
            obj.dGrid       = 1./obj.grid;
            mass            = ones(obj.grid);
            mom             = zeros([3, obj.grid]);
            mag             = zeros([3, obj.grid]);
            
            %--- Calculate Radial Distance ---%
            half      = floor(obj.grid/2);
            [X, Y, Z] = ndgrid(-half(1):half(1), -half(2):half(2), -half(3):half(3));
            distance  = sqrt(X.*X + Y.*Y + Z.*Z);
            clear X Y Z;
            
            %--- Determine Energy Distribution ---%
            radius          = max(5, round(0.06*min(obj.grid)));
            pressure        = 3*(obj.gamma-1)/(4*pi*(obj.dGrid(1)*radius)^3)*ones(obj.grid);
            pressure(distance > radius) = 1e-5;
            ener            = pressure./(obj.gamma-1);
        end
    
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
