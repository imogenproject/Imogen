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
        stencil;      % File name for the statics stencil (must be in data dir).        str
        staticType;   % Enumerated specification of how to apply static values.         str
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
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = true;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.7;
            obj.iterMax          = 10;
            obj.bcMode.x		 = 'trans';
            obj.bcMode.y		 = 'trans';
            obj.bcMode.z         = 'circ';
            obj.activeSlices.xy  = true;
            obj.ppSave.dim2      = 25;
            
            obj.staticType       = BowShockInitializer.PRIMAY_MODE;
            obj.stencil          = 'SmallSphere_800x256.mat';
            
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
        % Returns the initial conditions for a magnetic shock tube according to the settings for
        % the initializer.
        % USAGE: [mass, mom, ener, mag, statics, run] = getInitialConditions();
        
            %--- Initialization ---%
            pathStr = fileparts(mfilename);
            load([pathStr filesep 'data' filesep obj.stencil]);
            if ~exist('staticArray','var'); 
                error('Imogen:Initializer',['Static stencil array not found in loaded stencil ' ...
                      'file. Unable to continue']);
            end
            
            N = size(staticArray);
            if (obj.make3D(N) ~= obj.grid); 
                warning('Imogen:Initializer',['Current grid size property doesn''t match ' ...
                        'stencil size. Adjusting grid to match.']);
                obj.grid = N;
            end

            %--- Background Values ---%
            mass	= 0.125 * ones(obj.grid);
            mom		= zeros([3 obj.grid]);
            mag		= zeros([3 obj.grid]);
            ener	= (mass.^obj.gamma )./ (obj.gamma - 1);

            %--- Static Values ---%
            normalMass		= 0.125;
            normalEner      = (normalMass^obj.gamma)/(obj.gamma-1);
            statMass		= 1;
            statXMom		= 0.75;
            statEner		= (statMass^obj.gamma)/(obj.gamma-1) + 0.5*(statXMom*statXMom)/statMass;
            statics.values  = [0, normalMass, statMass, statXMom, normalEner, statEner];

            %--- Static Arrays ---%
            L = 5; U = 10;
            switch obj.staticType
                case BowShockInitializer.PRIMAY_MODE
                    statics.mass.s = uint8(2*staticArray);		statics.mass.s(L:U,:,:) = 3;

                    statics.mom.s = uint8(zeros([3 N]));
                    statics.mom.s(1,:,:,:) = staticArray;		statics.mom.s(1,L:U,:,:) = 4;
                    statics.mom.s(2,:,:,:) = staticArray;
                    statics.mom.s(3,:,:,:) = staticArray;

                    statics.ener.s =uint8(5*staticArray);       statics.ener.s(L:U,:,:)  = 6;

                case BowShockInitializer.FLUX_MODE
                    statics.mass.flux.s = staticArray;

                    statics.mom.flux.s			= uint8(zeros([3 N]));
                    statics.mom.flux.s(1,:,:,:) = staticArray;
                    statics.mom.flux.s(2,:,:,:) = staticArray;
                    statics.mom.flux.s(3,:,:,:) = staticArray;

                    statics.ener.flux.s = staticArray;

                    statics.mag.flux.s			= uint8(zeros([3 N]));
                    statics.mag.flux.s(1,:,:,:) = staticArray;
                    statics.mag.flux.s(2,:,:,:) = staticArray;
                    statics.mag.flux.s(3,:,:,:) = staticArray;
            
                case BowShockInitializer.FLUX_LR_MODE;

                    statics.mass.fluxl.s = staticArray;
                    statics.mom.fluxl.s			 = uint8(zeros([3 N]));
                    statics.mom.fluxl.s(1,:,:,:) = staticArray;
                    statics.mom.fluxl.s(2,:,:,:) = staticArray;
                    statics.mom.fluxl.s(3,:,:,:) = staticArray;
                    statics.ener.fluxl.s = staticArray;

                    statics.mass.fluxr.s = staticArray;
                    statics.mom.fluxr.s			 = uint8(zeros([3 N]));
                    statics.mom.fluxr.s(1,:,:,:) = staticArray;
                    statics.mom.fluxr.s(2,:,:,:) = staticArray;
                    statics.mom.fluxr.s(3,:,:,:) = staticArray;
                    statics.ener.fluxr.s = staticArray;
            end
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
