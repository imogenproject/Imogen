classdef SaveManager < handle
% The manager class responsible for handling saving/updating data actions. This is a singleton 
% class to be accessed using the getInstance() method and not instantiated directly.


%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
        SLICEFIELDS = {'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'cust'};   % Fields for slices.
		SLICELABELS = {'X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ', 'SP-XYZ'}; % Labels for slices.
	end
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %			P U B L I C  [P]
		%--- Results related logical states ---%
		updateUI;			% Specifies if the UI should be updated.					logical
		saveData;			% Specifies whether or not to save.							logical
		save1DData;			% Specifies whether or not to save 1D slices.				logical
		save2DData;			% Specifies whether or not to save 2D slices.				logical
		save3DData;			% Specifies whether or not to save 3D slices.				logical
		saveCustomData;		% Specifies whether or not to save custom slices.			logical
		done;				% Specifies if the code has reached the end for final save. logical
		previousUpdateTimes;% Saves the times of previous updates for each slice type.	double(5)
							%	1:	last time 1D data was saved.
							%	2:	last time 2D data was saved.
							%	3:  last time 3D data was saved.
							%	4:	last time Custom data was saved.
							%	5:	last time UI was updated.
		previousUpdateWallTimes; % Same as previousUpdateTimes except for wall time     double(5)
                                 % instead of simulation time.
        
        customSaveStr;		% An expression to be executed for custom saving			Str
        specialSaves1D;		% Array of specific iterations on which to save 1D data		int(N,1)
        specialSaves2D;		% Array of specific iterations on which to save 2D data		int(N,1)
        specialSaves3D;		% Array of specific iterations on which to save 3D data		int(N,1)

		FSAVE;				% Specifies if data should be saved to disk for the run.	logical
		PERSLICE;			% The percentages per save for each slice type:				double(4)
							%	1: % between 1D data saves.
							%	2: % between 2D data saves.
							%	3: % between 3D data saves.
							%	4: % between Custom data saves.
		
        ACTIVE;				% Specifies which slices should be saved.					logical(8,1)
        SLICE;				% Slice strings for each slice (e.g. ":,:,:")				Cell(8,1)

		SLICEINDEX;			% The indices at which to slice the grid for slice saves.  int(3,1)
        parent;				% parent manager                                           ImogenManager
        firstSave;          % Specifies if this is the first save action for the run.  logical
	end%PUBLIC

%===================================================================================================
	properties (SetAccess = private, GetAccess = private, Transient = true) %   P R I V A T E    [P]
        
	end%PRIVATE

	
	
	
	
	
	
%===================================================================================================
	methods %																	G E T / S E T	 [M]
		
%___________________________________________________________________________________________________ GS: SLICEINDEX
        function set.SLICEINDEX(obj,value)
			if (ndims(value) < 3),	value = [value 1]; end
            obj.SLICEINDEX = value;
		end

		function value = get.SLICEINDEX(obj); value = obj.SLICEINDEX; end
        
%___________________________________________________________________________________________________ GS: FSAVE
		function value = get.FSAVE(obj)
			value = obj.FSAVE;
		end
		
		function set.FSAVE(obj, value)
			obj.FSAVE = value;
			if (~obj.FSAVE), disp('Saving is disabled. Results in workspace only.'); end
		end
    
	end	
	
    
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]	

%___________________________________________________________________________________________________ preliminary
% Handles preliminary initialization of the SaveManager after all of the initialization settings 
% have been set. This function is meant to be called by the ImogenManager only.
        function preliminary(obj)
        
            % Skip if saving is inactive.
            if ~obj.FSAVE 
                return;
            end
            
            %--- Analyze grid directions for auto-slices ---%
            [~, indexMax] = max(obj.parent.gridSize);
            [~, indexMin] = min(obj.parent.gridSize);
           
            
            %--- Determine auto-slices ---%
            %           If no slices were specified by the user, the manager chooses the ones that
            %           appear most appropriate for the each dimension based on grid sizes.
            noSavesSpecified = ~any(obj.ACTIVE);
            
            if ~any(obj.ACTIVE(1:3)) && noSavesSpecified
                obj.ACTIVE(indexMax) = true;
            end

            if ~any(obj.ACTIVE(4:6))&& noSavesSpecified
                switch (indexMax + indexMin)
                    case 3; obj.ACTIVE(4) = true;
                    case 4; obj.ACTIVE(5) = true;
                    case 5; obj.ACTIVE(6) = true;
                end
            end
            
            if ~obj.ACTIVE(7) && noSavesSpecified
                obj.ACTIVE(7) = true;
            end

            obj.firstSave = true;
        end

%___________________________________________________________________________________________________ postliminary
        function postliminary(obj) %#ok<MANU>
        end
        
%___________________________________________________________________________________________________ logPrint
% Prints information to the standard output as well as a log file.
        function logPrint(obj, printLine, varargin)
            if labindex == 1
                fid       = fopen([obj.parent.paths.save, filesep, 'logfile.out'],'a');
                fprintf(fid, printLine, varargin{:});
                fprintf(printLine, varargin{:});
                fclose(fid);
            end
            labBarrier();
        end
        
%___________________________________________________________________________________________________ saveIniSettings
% Saves the ini structure to disk so that it can be reused later to restart a run if necessary.
        function saveIniSettings(obj, ini)
           if labindex == 1 && obj.FSAVE && ~isempty(ini)
               %--- Save ini Structure ---%
               %            Saves the entire ini structure to a mat file for later reloading and 
               %            reuse.
               save([obj.parent.paths.save filesep 'ini_settings.mat'],'ini');
               
                %--- Create run.log Entry ---%
                %               Writes basic run information to the run.log file for later 
                %               reference. This includes fixed basic information as well as any
                %               fields listed in the ini.logProperties.
                paths    = obj.parent.paths;
                run      = obj.parent;
                
                data     = paths.saveFolder;
                data     = [data, '\n\tSaved at: ', strrep(paths.save, '\', '/')];
                data     = [data, '\n\tVersion: ', run.detailedVersion];
                data     = [data, '\n\tStarted: ', datestr(run.time.startTime, ...
                                                           'HH:MM on mmm dd, yyyy')];
                data     = [data, '\n\tInfo: ', strrep(run.about, '---+ ', '')];

                for i=1:length(ini.logProperties)
                    prop = ini.logProperties{i};
                    data = [data, '\n\t', prop, ': ',  ImogenRecord.valueToString(ini.(prop))]; %#ok<AGROW>
                end
                
                data     = strcat(data, '\n');
                
                fid      = fopen(strcat(paths.results, '/run.log'), 'a');
                fprintf(fid, strcat(repmat('-', [1,80]), '\n', data));
                fclose(fid);
                
                fid      = fopen(strcat(paths.save, '/run.log'), 'w');
                fprintf(fid, data);
                fclose(fid);
           end
           labBarrier();
        end

%___________________________________________________________________________________________________ updateDataSaves
% Updates the data saves depending on whether the iterations or the time is expected to complete
% first.
        function updateDataSaves(obj)
            time = obj.parent.time;
            time.updateWallTime();
            
            [~, index] = max([time.iterPercent, time.timePercent, time.wallPercent]);
            switch index
                case 1
                    obj.updateIterationDataSaves(time);
                case 2
                    obj.updateTimeDataSaves(time);
                case 3
                    obj.updateWallDataSaves(time);
            end
            
			obj.updateOverrides(time.iteration);
			obj.scheduleSave();
        end
        		
%___________________________________________________________________________________________________ updateOverrides
		function updateOverrides(obj,iteration)
            special = (obj.done || obj.firstSave);
            i       = iteration;
            
			if (~isempty(obj.specialSaves1D) && any(~any(obj.specialSaves1D - i, 1)) || special)	
				obj.save1DData = true; 
			end
			
			if (~isempty(obj.specialSaves2D) && any(~any(obj.specialSaves2D - i, 1)) || special) 
				obj.save2DData = true; 
			end
			
			if (~isempty(obj.specialSaves3D) && any(~any(obj.specialSaves3D - i, 1)) || special)
				obj.save3DData = true; 
			end
		end
				
%___________________________________________________________________________________________________ getSaveSlice
% This function returns the sliced data from the input array. If the array is codistributed, the 
% slice is gathered to labindex 1 and returned only to that lab. All other labs see result = [].
%>> array		Data array to be sliced.                                        double(?)
%>> type		Type of slice to be made on the array (1-8).                    int
%<< result		Resulting slice (exists only on labindex 1).                    double(?)
		function result = getSaveSlice(obj,array,type)
        
			N = size(array);         
            if length(N) < 3;   N = [N 1]; end
            
			switch type
				case 1;			i = {1:N(1),obj.SLICEINDEX(2),obj.SLICEINDEX(3)};
				case 2;			i = {obj.SLICEINDEX(1),1:N(2),obj.SLICEINDEX(3)};
				case 3;			i = {obj.SLICEINDEX(1),obj.SLICEINDEX(2),1:N(3)};
					
				case 4;			i = {1:N(1),1:N(2),obj.SLICEINDEX(3)};
				case 5;			i = {1:N(1),obj.SLICEINDEX(2),1:N(3)};
				case 6;			i = {obj.SLICEINDEX(1),1:N(2),1:N(3)};
					
				case {7, 8};	i = { 1:N(1), 1:N(2), 1:N(3) };
			end

			%--- Gather slices to lab 1 ---%
			if iscodistributed(array)
				result = squeeze( gather(array(i{:}),1) );
			else
				result = squeeze( array(i{:}));
			end
		end
		
	end%PUBLIC
	
%===================================================================================================	
	methods (Access = private) %												P R I V A T E    [M]
	
%___________________________________________________________________________________________________ scheduleSave
% Tests save conditions and determines the state of the saveData property.
		function scheduleSave(obj)
			obj.saveData = ( obj.save1DData || obj.save2DData || obj.save3DData ...
                            || obj.saveCustomData || obj.done) && obj.FSAVE;
		end
		
%___________________________________________________________________________________________________ SaveManager
% Creates a new SaveManager instance.
        function obj = SaveManager() 
            obj.SLICE                   = cell(8,1);
            obj.SLICEINDEX              = ones(1,3);
            obj.ACTIVE                  = false(1,8);
			obj.previousUpdateTimes     = zeros(1,5);
            obj.previousUpdateWallTimes = zeros(1,5);
			obj.saveData                = true;
			obj.save1DData              = true;
			obj.save2DData              = true;
			obj.save3DData              = true;
			obj.saveCustomData          = false;
			obj.done                    = false;
        end

%___________________________________________________________________________________________________ updateWallDataSaves
        function updateWallDataSaves(obj, time)
			obj.done    = (time.wallTime >= time.WALLMAX);
			wallUpdates = (time.wallPercent - 100*obj.previousUpdateWallTimes/time.WALLMAX) ...
                            >= [obj.PERSLICE 10];
                        
			obj.save1DData      = wallUpdates(1);
            obj.save2DData      = wallUpdates(2);
            obj.save3DData      = wallUpdates(3);
            obj.saveCustomData	= wallUpdates(4);
            obj.updateUI        = wallUpdates(5);
            
            for i=1:length(wallUpdates)
                if wallUpdates(i)
                    obj.previousUpdateTimes(i)      = time.time;
                    obj.previousUpdateWallTimes(i)  = time.wallTime;
                end
            end
        end
        
%___________________________________________________________________________________________________ updateTimeDataSaves
		function updateTimeDataSaves(obj, time)
			obj.done    = (time.time >= time.TIMEMAX);
			timeUpdates = (time.timePercent - 100*obj.previousUpdateTimes/time.TIMEMAX) ...
                            >= [obj.PERSLICE 10];
                        
			obj.save1DData      = timeUpdates(1);
            obj.save2DData      = timeUpdates(2);
            obj.save3DData      = timeUpdates(3);
            obj.saveCustomData	= timeUpdates(4);
            obj.updateUI        = timeUpdates(5);
            
            for i=1:length(timeUpdates)
                if timeUpdates(i)
                    obj.previousUpdateTimes(i)      = time.time;
                    obj.previousUpdateWallTimes(i)  = time.wallTime;
                end
            end
		end
		
%___________________________________________________________________________________________________ updateIterationDataSaves
		function updateIterationDataSaves(obj, time)
			obj.done = (time.iteration >= time.ITERMAX);
			modVal = max(round(time.ITERMAX*([obj.PERSLICE/100 0.1])), 1);
			modVal = min(modVal, max(floor(time.ITERMAX/4),1)); % At least 4 saves if MAXITER >= 4.
			modVal = mod(time.iteration, modVal);               % Critical loop test.
			
			obj.save1DData		= ~logical(modVal(1));
            obj.save2DData		= ~logical(modVal(2));		
            obj.save3DData		= ~logical(modVal(3));
            obj.saveCustomData	= ~logical(modVal(4)); 
            obj.updateUI		= ~logical(modVal(5));
            
            for i=1:length(modVal)
                if ~logical(modVal(i))
                    obj.previousUpdateTimes(i)      = time.time;
                    obj.previousUpdateWallTimes(i)  = time.wallTime;
                end
            end
		end
        
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]
		
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the SaveManager class, or creates one if it doesn't exist.
		function singleObj = getInstance()
			persistent instance;
			if isempty(instance) || ~isvalid(instance) 
				instance = SaveManager();
			end
			singleObj = instance;
		end
	  
	end%STATIC
	

end%CLASS
