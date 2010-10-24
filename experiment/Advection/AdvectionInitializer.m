classdef AdvectionInitializer < Initializer
    
%___________________________________________________________________________________________________ 
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
    end %PROTECTED
	
	
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]
        
        function obj = AdvectionInitializer(input)
            obj                 = obj@Initializer();
            obj.gamma           = 5/3;
            obj.runCode         = 'ADVEC';
            obj.info            = 'Advection test.';
            obj.mode.fluid      = true;
            obj.mode.magnet     = false;
            obj.mode.gravity    = false;
            obj.cfl             = 0.5;
            obj.iterMax         = 1000;
            obj.ppSave.dim1     = 10;
            obj.ppSave.dim3     = 25;
            obj.activeSlices.xy = true;
            obj.activeSlices.xyz= true;
            obj.activeSlices.x  = true;
            
            obj.bcMode          = ENUM.BCMODE_CONST;
            
            obj.operateOnInput(input);
        end
        
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
        
        function [mass, mom, ener, mag, statics] = calculateInitialConditions(obj)
            statics  = [];
            mass     = 1*ones(obj.grid);
            mom      = zeros([3, obj.grid]);
            mag      = zeros([3, obj.grid]);
            ener     = ones(obj.grid);
            
            
            massBump = 2;
            mach     = 0.75;
            g        = obj.gamma;
            speed    = mach*sqrt( (g + ( 2/((g-1) * 2) - g * (g-1)/2)) ...
                          / ( 1/(g-1) + g*(mach*mach)/2 ) );
 
            range    = floor(0.2*obj.grid(1)):ceil(0.25*obj.grid(1));
            
            mass(range,:,:)  = massBump;
            mom(1,range,:,:) = massBump*speed;
            ener(range,:,:)  = ener(range,:,:) + squeeze(mom(1,range,:,:));
        end
        
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
	end%PROTECTED
	
end%CLASS
    