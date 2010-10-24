classdef Paths < handle
% The storage class for all path related values needed to load/save data for an imogen run. 

%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
        DEFAULT      = 'def';		% ENUMERATION: "Defaulted" warning
        OVERRIDE     = 'over';		% ENUMERATION: "Override" warning
        RADIX_BUFFER = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'; % Base 36 Radix
    end%CONSTANT

%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %			P U B L I C  [P]
		imogen;				% Top-level imogen run path (where imogen.m is located).		str
		results;			% Top-level results path.										str
		indexPadding;		% Number of digits required for padded numbers.                 int
		hostName;			% Name of the host running Imogen.                              str
        saveFolder;         % Name of folder in which data is saved.                        str
        containerFolder;    % Name of container folder in which data will be saved.         str
        runCode;            % Run code for the type of simulation being executed.           str
        alias;              % Unique identifier for the run.                                str
	end%PUBLIC
	
%===================================================================================================
    properties (Dependent = true, SetAccess = public) %                        D E P E N D E N T [P]
        container;			% path to the container directory saving data.					str
        image;				% image subdirectory path for the save.							str
        save;				% path to the save directory for the run.						str
    end %DEPENDENT    
	
%===================================================================================================
	properties (SetAccess = private, GetAccess = private, Transient = true) %	P R I V A T E    [P]
    end
	
	
%===================================================================================================
	methods %																	G E T / S E T	 [M]	

%___________________________________________________________________________________________________ Paths
		function obj = Paths() 
            obj.alias = '';
        end
        
%___________________________________________________________________________________________________ GS: save
        function result = get.save(obj)
        % Access to the folder where the run data is to be stored.
            result = strcat( obj.container, filesep, obj.saveFolder);
        end
        
%___________________________________________________________________________________________________ GS: container
        function result = get.container(obj)
        % Access to the monthly container folder in the top level results directory.
            result = strcat( obj.results, filesep, obj.containerFolder);
        end
        
%___________________________________________________________________________________________________ GS: image
        function result = get.image(obj)
        % Access to the image subdirectory in the save directory.
            result = strcat(obj.save, filesep, 'images');
        end
        
    end%GET/SET
    
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]

%___________________________________________________________________________________________________ initialize
        function initialize(obj)
            timeManager             = TimeManager.getInstance();
            startTime               = timeManager.startTime;
            
            obj.containerFolder     = datestr(startTime,'mmmyy');
            obj.makePathUnique(startTime);
            obj.printHostVariables(startTime);
        end
        
%___________________________________________________________________________________________________ iterationToString
		function result = iterationToString(obj, iteration)
        % Pads an iteration value with zeros for consistent save length of intermediate slice data.
			result = Paths.paddedNumber(iteration, obj.indexPadding);
		end
		
	end%PUBLIC
	
%===================================================================================================
    methods (Access = private) %												P R I V A T E    [M]

%___________________________________________________________________________________________________ printHostVariables
		function printHostVariables(obj, startTime)
            fprintf('\n\nRun started at: %s (%s)\n', datestr(startTime), obj.saveFolder);
			fprintf('\tRunning on %s\n',            obj.hostName);
			fprintf('\tImogen directory: %s\n',     obj.imogen);
			fprintf('\tOutput directory: %s\n\n\n',	obj.results);
		end
        
%___________________________________________________________________________________________________ makePathUnique
% Creates a unique folder name on the path in cases where naming conflicts occur. This happens when
% multiple runs are started at the same time (within the same minute).
%
%>> folderName		fully resolved path to the folder to check for name conflicts			str
%<< updatedName		unique name for the folder differing from input if conflict existed		str
        function makePathUnique(obj, startTime)
            
            %--- Load the UID if it exists ---%
            try
                runID = evalin('base','alias;');
            catch ME
                runID = '';
            end
            
            if isempty(runID)
                runID = obj.alias;
            end
            
            timeCode = strcat(obj.encodeBase36(startTime(2)), obj.encodeBase36(startTime(3)), ...
                              obj.encodeBase36(startTime(4)));

            suffix = '';
            if ~isempty(runID)
                suffix = strcat('_', runID);
            end
            
            for i=0:10000
                id = obj.encodeBase36(i-1);

                folderName = strcat(obj.runCode, '_', timeCode, id, suffix);
                if ~exist(strcat( obj.container, filesep, folderName), 'dir')
                    obj.saveFolder = folderName;
                    return
                end
            end
        end 
    
%___________________________________________________________________________________________________ encodeBase36
% Encodes the number into a base 36 string.
        function result = encodeBase36(obj, number)
            %--- Initialize ---%
            if number == 0
                result = '0';
            else
                result = '';
            end

            %--- Create Base 36 Number ---%
            while number > 0
                index  = mod(number, length(obj.RADIX_BUFFER));
                result = strcat(obj.RADIX_BUFFER(index + 1), result);
                number = floor(number/36);
            end
        end
        
    end%PRIVATE
    
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]

%___________________________________________________________________________________________________ paddedNumber
% Converts a number to a string and pads it with leading zeros for consistent length for use in i/o
% naming conventions.
		function result = paddedNumber(number, padLength)
			result = num2str(number);
			zeroPad = padLength - length(result);
			if (zeroPad > 0)
				zeroStr = '00000000000000';
				result = strcat(zeroStr(1:zeroPad),result);
			end
		end
		
	end%PROTECTED
	
end
		