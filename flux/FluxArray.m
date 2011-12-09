classdef FluxArray < InitializedArray
% Array class for flux objects connected to primary variable arrays.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
		FLUX	= 'flux';
		FLUXL	= 'fluxl';
		FLUXR	= 'fluxr';
		DFLUXLL	= 'dfluxll';
		DFLUXLR = 'dfluxlr';
		FREEZE  = 'freeze';
	end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
        setArrayLimited;
	end%PUBLIC
    
%===================================================================================================
    properties (Dependent = true) %											   D E P E N D E N T [P]
		fluxType;
    end %DEPENDENT
    
    
%===================================================================================================
    methods %																	  G E T / S E T  [M]

%___________________________________________________________________________________________________ GS: fluxType
		function result = get.fluxType(obj)
            try
                result = obj.id{2};
            catch MERR
                result = 'unspecified';
            end
		end
		
	end%GET/SET
    
    
%===================================================================================================
    methods %																		P U B L I C  [M]
        
%___________________________________________________________________________________________________ FluxArray
        function obj = FluxArray(component, id, run)
            obj = obj@InitializedArray(component, id, run);
            
            obj.setArrayLimited = @obj.setArrayLimited_Vanleer;
        end
        
	end%PUBLIC

%===================================================================================================
    methods (Access = protected) %											   P R O T E C T E D [M]

%___________________________________________________________________________________________________ setArrayLimited_Vanleer
        function setArrayLimited_Vanleer(obj, dFluxL, dFluxR)
            obj.array = vanleerLimiter(obj.array, dFluxL.array, dFluxR.array);   
		end

%___________________________________________________________________________________________________ setArrayLimited_Vanleer_Statics
        function setArrayLimited_Vanleer_Statics(obj, dFluxL, dFluxR)
            obj.array = vanleerLimiter(obj.array, dFluxL.array, dFluxR.array);
            obj.array = cleanStaticCells(obj);
		end
		
	end%PROTECTED
end%CLASS
