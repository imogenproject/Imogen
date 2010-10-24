classdef OrszagTangVortexInitializer < Initializer
% Creates initial conditions for an Orszag-Tang vortex test. This is THE canonical test for high
% fidelity magnetic fluxing, as it simulates the transition from sub-sonic to super-sonic mangetic
% flows for non-axis-aligned propagation vectors. There is no analytical solution to this problem,
% but it has been tested thoroughly by a number of different MHD codes to the point where a
% consensus has formed regarding its behavior. This particular problem has been setup to conform to
% the results of the Ryu paper on magnetohydrodynamics codes.
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
        
%___________________________________________________________________________________________________ OrszagTangVortexInitializer
        function obj = OrszagTangVortexInitializer(input)
            obj                  = obj@Initializer();
            obj.gamma            = 5/3;
            obj.runCode          = 'OTVortex';
            obj.info             = 'Orszag-Tang vortex trial.';
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = true;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.35;
            obj.iterMax          = 1500;
            obj.timeMax          = 0.48;
            obj.bcMode.x		 = 'circ';
            obj.bcMode.y		 = 'circ';
            obj.bcMode.z         = 'circ';
            obj.activeSlices.xy  = true;
            obj.ppSave.dim2      = 25;
            
            obj.operateOnInput(input, [256, 256, 1]);
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
            mass            = 25/(36*pi)*ones(obj.grid);

            [y,x]           = meshgrid(0:1:(obj.grid(2)-1),0:1:(obj.grid(1)-1));
            mom             = zeros([3 obj.grid]);
            mom(1,:,:)      = - mass .* sin( 2*pi*y/(obj.grid(2)-1) );
            mom(2,:,:)      =   mass .* sin( 2*pi*x/(obj.grid(1)-1) );

            mag0            = 1/sqrt(4*pi);
            mag             = zeros([3 obj.grid]);
            mag(1,:,:)      = - mag0*sin( 2*pi*y/(obj.grid(2)-1) );
            mag(2,:,:)      =   mag0*sin( 4*pi*x/(obj.grid(1)-1) );

            ener        = 5/(12*pi)/(obj.gamma - 1) ...                % internal
                            + 0.5*squeeze(sum(mom.*mom,1)) ./ mass ...  % kinetic
                            + 0.5*squeeze(sum(mag.*mag));               % magnetic
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
