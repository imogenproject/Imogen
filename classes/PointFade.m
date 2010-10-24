classdef PointFade < handle
% This class creates a point fade array that can be applied to Imogen arrays to enforce soft static 
% conditions.
		
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
        array;
        location;
        activeList;
        fluxes;
    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
    end %PROTECTED
	
	
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]
        
%___________________________________________________________________________________________________ PointFade
        function obj = PointFade(grid, location, fadeSize)
            left        = 1 - location;
            right       = grid - location;
            [x y z]     = ndgrid(left(1):right(1), left(2):right(2), left(3):right(3));
            distSq      = max(fadeSize.*fadeSize - (x.*x + y.*y + z.*z), 0);
            distSq      = 1 - distSq/maxFinderND(distSq);
            distSq      = pchip([-1 0 1 2], [1 1 0 0], reshape(distSq, [1 numel(distSq)]) );
            obj.array   = reshape(distSq, grid);
            obj.location= num2cell(location);
            obj.fluxes  = false;
        end
        
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
        
%___________________________________________________________________________________________________ fadeArray
        function result = fadeArray(obj, inArray, fadeValue)
            if isempty(fadeValue); fadeValue = inArray(obj.location{:}); end
            result = 0.95*inArray.*obj.array + inArray.*(1-obj.array);
        end
        
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
	end%PROTECTED
	
end%CLASS
