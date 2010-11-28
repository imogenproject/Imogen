classdef MagnetManager < handle
% Manages magnetic field routines and related settings. This is a singleton class to be accessed 
% using the getInstance() method and not instantiated directly.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
	end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %			P U B L I C  [P]
		ACTIVE;				% specifies the magnetic fluxing solver state			logical
		limiter;            % Flux limiters to use for each flux direction.         cell(3)
	end%PUBLIC

%===================================================================================================
    properties (SetAccess = public, GetAccess = private) %						   P R I V A T E [P]
		parent;			% Parent manager											ImogenManager
    end %PRIVATE

	
	
	
	
	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
		
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = private) %												P R I V A T E    [M]
		
%___________________________________________________________________________________________________ MagnetManager
% Creates a new MagnetManager instance.
		function obj = MagnetManager() 
			obj.ACTIVE = false;
		end
		
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
		
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the MagnetManager class, or creates one if none have
% been initialized yet.
		function singleObj = getInstance(create)
			persistent instance;
			if isempty(instance) || ~isvalid(instance) 
				instance = MagnetManager(); 
			end
			singleObj = instance;
		end
		
	end%STATIC
	
	
end%CLASS