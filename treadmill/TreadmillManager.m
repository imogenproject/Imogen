classdef TreadmillManager < handle
% The manager class responsible for handling treadmill related actions and data. This is a singleton 
% class to be accessed using the getInstance() method and not instantiated directly.
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %			P U B L I C  [P]
        history;		% history of treadmill action over course of run			int(?)
        last;			% previous treadmill action									int
        DIRECTION;		% spatial direction along which to treadmill (1,2,3)		int
	end%PUBLIC

%===================================================================================================
    properties (Dependent = true) %											   D E P E N D E N T [P]
		ACTIVE;			% Specifies the enabled state for treadmilling				logical
	end%PUBLIC
	
%===================================================================================================
	properties (SetAccess = public, GetAccess = private, Transient = true) %	P R I V A T E    [P]
		pHistoryIndex;	% index for the next write into the history array			int
		parent;			% parent manager											ImogenManager
		pActive;		% stores the active state for treadmilling					logical
	end%PRIVATE
	
	
	
	
	
	
%===================================================================================================
	methods %																	G E T / S E T	 [M]
		
%___________________________________________________________________________________________________ GS: ACTIVE
% Specifies whether or not treadmilling is active for the run.
		function set.ACTIVE(obj,value)
			obj.pActive = value;
			if (value); obj.history = zeros(1000,1); end
		end
		function value = get.ACTIVE(obj);	value = obj.pActive; end
			
	end%GET/SET	

%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]
		        
%___________________________________________________________________________________________________ toStruct
% Converts the TreadmillManager to a structure for saving and non-class related actions.
% # result		The structure created by the conversion.								Struct
		function result = toStruct(obj)
			if (obj.pActive) 
				result.direction = obj.DIRECTION;
				result.history = obj.history(1:max(obj.pHistoryIndex-1,1));
			else
				result = [];
			end
		end
		
%___________________________________________________________________________________________________ appendHistory
% Appends a "number of cells treadmilled" value to the history array.
% * value		The number of cells to treadmill for the current iteration				int
		function appendHistory(obj,value)
			if (obj.pHistoryIndex > length(obj.history))
				obj.history = [obj.history; zeros(500,1)];
			end
			obj.history(obj.pHistoryIndex) = value;
			obj.pHistoryIndex = obj.pHistoryIndex + 1;
		end
			
	end%PUBLIC
    
%===================================================================================================	
	methods (Access = private) %												P R I V A T E    [M]
		
%___________________________________________________________________________________________________ TreadmillManager
% Creates a new TreadmillManager instance.
		function obj = TreadmillManager() 
            obj.pActive		  = false;
            obj.last		  = -1;
            obj.DIRECTION	  = 1;
			obj.pHistoryIndex = 1;
		end
		
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
		
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the TreadmillManager class, or creates one if none have
% been initialized yet.
		function singleObj = getInstance()
			persistent instance;
			if isempty(instance) || ~isvalid(instance) 
				instance = TreadmillManager();
			end
			singleObj = instance;
		end
		
	end%STATIC
	
    
end