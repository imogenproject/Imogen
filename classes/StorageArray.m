classdef StorageArray < InitializedArray
% The class used for first order fluxing. It contains the flux array variables necessary for the 
% predictor step in the fluxing scheme and is called a storage array because it stores the 
% intermediate value guess supplied to the second-order corrector fluxing step.

%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
		STORE	= 'store';
    end%CONSTANT

%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
        fluxR;		% right flux array											FluxArray		H
        fluxL;		% left flux array											FluxArray		H
        wArray;		% auxiliary array for fluxing								double
		threshold;
		thresholdArray;
    end %PUBLIC
	
	
	
%===================================================================================================
    methods %																		P U B L I C  [M]
		
		function result = get.thresholdArray(obj)
			if (obj.threshold > 0)
				result = obj.pArray .* (obj.pArray > obj.threshold);
			else
				result = obj.pArray;
			end
		end
		
%___________________________________________________________________________________________________ StorageArray
        function obj = StorageArray(component, id, run)
            obj = obj@InitializedArray(component, id, run);
			
			if length(obj.id) < 2
				fullID = cell(1,2);
				fullID{1} = obj.id{1};
				fullID{2} = StorageArray.STORE;
				obj.id = fullID;
            end

			if strcmpi(id,ENUM.MASS)
				obj.threshold = run.fluid.MASS_THRESHOLD;
			else
				obj.threshold = 0;
			end
			
			obj.arrayINI(component, obj.id, run);
        end
        
%___________________________________________________________________________________________________ cleanup
        function cleanup(obj)
%             obj.flux.cleanup();
            obj.fluxL.cleanup();
            obj.fluxR.cleanup();
			obj.wArray = [];
        end
        
    end %PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
		
%___________________________________________________________________________________________________ arrayINI
        function arrayINI(obj, component, id, run)
            %Prepare flux arrays
            obj.fluxL  = FluxArray(component, {id{1}, FluxArray.FLUXL, id{2}}, run);
            obj.fluxR  = FluxArray(component, {id{1}, FluxArray.FLUXR, id{2}}, run);
		end
		
    end %PROTECTED
    
end %CLASS
