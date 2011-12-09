classdef InitializedArray < ImogenArray
% Abstract class for arrays that are intialized in situ during runs the first time their arrays are
% set instead of prior to starting the run using input data. A good example would be a flux, which
% isn't calculated until the Flux routine determines the flux values during the first iteration.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
    end %PROPERTIES
    
%===================================================================================================
    properties (SetAccess = private, GetAccess = public) %						   P R I V A T E [P]
    end %PROPERTIES
	
	
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]       
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
		
%___________________________________________________________________________________________________ InitializedArray
        function obj = InitializedArray(component, id, run)			
            obj                 = obj@ImogenArray(component, id, run);
			obj.pUninitialized  = true;
        end
		
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]

%___________________________________________________________________________________________________ applyInitialize
% Overrides the applyInitialize method in ImogenArray for the InitializedArray case, where 
% initialization doesn't occur during object construction.
        function applyInitialize(obj)
            obj.initializeShiftingStates();
            obj.initializeBoundingEdges();
            obj.pUninitialized = false;
        end
        
	end%PROTECTED
		
end%CLASS
