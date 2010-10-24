classdef GravityArray < InitializedArray
% Array class that encapsulates the gravitational potential and related fuctionality for inclusion
% in the source functions of Imogen.
	
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
		
%___________________________________________________________________________________________________ GravityArray
% Creates a new GravityArray object.
        function obj = GravityArray(id, run, statics)
			obj = obj@InitializedArray(ENUM.SCALAR, id, run, statics);
			if isempty(id); return; end
			obj.array = zeros(run.gridSize);
		end   
		
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
		
	end%PROTECTED
		
end%CLASS
