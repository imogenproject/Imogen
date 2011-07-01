classdef KelvinHelmholtzInitializer < Initializer
% Creates initial conditions for a Kelvin-Helmholtz instability, where two anti-parallel shearing 
% flows are seeded with noise along the shearing barrier. The perturbation of the barrier causes
% the well known Kelvin-Helmholtz instabiliity phenomena to form. It should be noted that to 
% prevent a large shock wave from forming across the barrier the two regions, while having differing
% mass density values, have equal pressures. This makes sense given the usual observation of this
% phenomena in clouds, where such conditions exist.
%
% Unique properties for this initializer:
%     direction      % enumerated orientation of the baseline flow.              str
%     massRatio      % ratio of (low mass)/(high mass) for the flow regions.     double
%     mach           % differential mach number between flow regions.            double
%     perturb        % specifies if flow should be perturbed.                    logical
    
        
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
        X = 'x';
        Y = 'y';
        Z = 'z';
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %                           P U B L I C  [P]
        direction;      % enumerated orientation of the baseline flow.              str
        massRatio;      % ratio of (low mass)/(high mass) for the flow regions.     double
        mach;           % differential mach number between flow regions.            double
        perturb;        % specifies if flow should be perturbed.                    logical
    end %PUBLIC

%===================================================================================================
    properties (Dependent = true) %                                            D E P E N D E N T [P]
    end %DEPENDENT
    
%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %                P R O T E C T E D [P]
    end %PROTECTED
    
    
    
    
    
%===================================================================================================
    methods %                                                                     G E T / S E T  [M]
        
%___________________________________________________________________________________________________ KelvinHelmholtzInitializer
        function obj = KelvinHelmholtzInitializer(input)
            obj = obj@Initializer();
            obj.gamma            = 5/3;
            obj.runCode          = 'KelHelm';
            obj.info             = 'Kelvin-Helmholtz instability trial.';
            obj.mode.fluid		 = true;
            obj.mode.magnet		 = false;
            obj.mode.gravity	 = false;
            obj.cfl				 = 0.7;
            obj.iterMax          = 1500;
            obj.activeSlices.xy  = true;
            obj.ppSave.dim2      = 25;
            
            obj.direction        = KelvinHelmholtzInitializer.X;
            obj.massRatio        = 8;
            obj.mach             = 0.25;
            obj.perturb          = true;
            
            obj.operateOnInput(input);
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
            statics = [];
            indeces = cell(1,3);
            for i=1:3;  indeces{i} = 1:obj.grid(i); end
                       
            mass	= ones(obj.grid);
            mom		= zeros([3 obj.grid]);
            mag     = zeros([3 obj.grid]);
            speed	= speedFromMach(obj.mach, obj.gamma, 1, 1/(obj.gamma-1), 0);

            half    = ceil(obj.grid/2);
            fields  = {obj.X, obj.Y, obj.Z};
            for i=1:3
                if strcmpi(obj.direction,fields{i})
                    index            = indeces;
                    if (i ~= 1);    index{1} = half(1):obj.grid(1);
                    else            index{2} = half(2):obj.grid(2);
                    end
                    
                    mass(index{:}) = 1/obj.massRatio*mass(index{:});
                    mom(i,:,:,:)     = speed*mass;
                    mom(i,index{:}) = -squeeze(mom(i,index{:}));
                    obj.bcMode.(fields{i}) = 'circ';
                else obj.bcMode.(fields{i}) = 'const';
                end
            end
            
            if obj.perturb 
                halfPlus    = min(half+3,obj.grid);
                halfMinus   = max(half-3,1);
                for i=1:3;  index{i} = halfMinus(i):halfPlus(i); end
                mass(index{:}) = maxFinderND(mass);
            end
        
            ener	= (maxFinderND(mass)^obj.gamma)/(obj.gamma - 1) ...     % internal
                      + 0.5*squeeze(sum(mom.*mom,1))./mass ...              % kinetic
                      + 0.5*squeeze(sum(mag.*mag,1));                       % magnetic
        end
        
    end%PROTECTED
        
%===================================================================================================    
    methods (Static = true) %                                                     S T A T I C    [M]
    end
end%CLASS
