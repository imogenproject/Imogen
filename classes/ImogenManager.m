classdef ImogenManager < handle
% This is the primary management class for the Imogen code. It handles storage of all of the 
% un-typed variables and stores the other, specific managers. After being initialized, this object
% is passed throughout the code as the source for dynamically accessed data. This is a singleton 
% class to be accessed using the getInstance() method and not instantiated directly.
    
%===================================================================================================
    properties (Constant = true, Transient = true) %                            C O N S T A N T  [P]
        DEFAULT  = 'def';        % ENUMERATION: "Defaulted" warning
        OVERRIDE = 'over';        % ENUMERATION: "Override" warning
    end%CONSTANT
    
%===================================================================================================
    properties (SetAccess = public, GetAccess = public, Transient = true) %         P U B L I C  [P]
        about;          % Information about run (includes run code).                    str
        version;        % Imogen version for the run.                                   double
        detailedVersion;% Version of Imogen for the run including detailed information. str
        matlab;         % Information regarding matlab version being used.              struct
        info;           % Information generated during initialization and run.          cell(?,2)
        iniInfo;        % Information regarding data initial condition settings.        str
        DEBUG;          % Specifies if Imogen is run in debug mode.                     logical    
        DGRID;          % Grid spacings each spatial direction.                         cell(3)
        MINDGRID;       % Smallest grid spacing for each spatial direction.             double(3)
        GAMMA;          % Polytropic index for the run.                                 double
        PROFILE;        % Specifies state of profiler for the run.                      logical
        gridSize;       % Grid size for the arrays.                                     int(3)
        paths;          % Contains various paths needed for saving data.                Paths
        fades;          %

        useGPU;         % if true, Imogen is running on a GPU device - some changes in behavior necessary.
        pureHydro;      % if true, stores no magnetic information; 
        
        %--- Manager Classes ---%
        bc;             % Manages boundary conditions.                              BCManager
        image;          % Manages image generation and saving.                      ImageManager
        save;           % Manages updating and saving data.                         SaveManager
        time;           % Manages temporal actions.                                 TimeManager
        gravity;        % Manages potential solver.                                 GravityManager
        treadmill;      % Manages treadmill actions.                                TreadmillManager
        fluid;          % Manages fluid routines.                                   FluidManager
        magnet;         % Manages magnetic routines.                                MagnetManager
        parallel;       % Manages parallel processing.                              ParallelManager
    end%PUBLIC
    
%===================================================================================================
    properties (SetAccess = private, GetAccess = private, Transient = true) %    P R I V A T E   [P]
        infoIndex;      % Stores current index for the info cell array.             double
        pAbortTime;     % Time of last abort check.                                 serial date #
        pAbortFile;     % Name of the abort check file.                             str
        pWarnings;      % Warning information generated during ini and run.         str
        pNotes;         % Detailed user generated information about the run.        str       
    end%PRIVATE
    
%===================================================================================================
    properties (Dependent = true, SetAccess = public) %                        D E P E N D E N T [P]
        warnings;       % Accesses formatted warning statements logged during run.      str
        notes;          % String of information regarding all manner of run activity.   str
    end %DEPENDENT
    
    
    
    
    
    
%===================================================================================================
    methods %                                                                    G E T / S E T   [M]    
        
%___________________________________________________________________________________________________ GS: gridSize
% Stores the grid size parameters for the arrays.
        function set.gridSize(obj,value)
            if ( length(value) < 3), value = [value 1]; end
            obj.gridSize = value;
        end
        
%___________________________________________________________________________________________________ GS: notes
% The "notes" variable. On get, this variable is formatted for wiki syntax output.
        function result = get.notes(obj)
            if isempty(obj.pNotes)
                result = '\n   * ENTER INFORMATION ABOUT THIS RUN HERE.';
            else
                result = strrep(obj.pNotes,  '\n', '\n   * ' );
                result = strrep(result, sprintf('\n'), '\n   * ' );
                result = sprintf('\n   * %s',result);
            end
        end
        
        function set.notes(obj,value); obj.pNotes = value; end
        
%___________________________________________________________________________________________________ GS: warnings
% Returns the formatted warnings statements in wiki syntax
        function result = get.warnings(obj)
            if isempty(obj.pWarnings)
                result = '\n   * No warnings logged.';
            else
                result = strrep(obj.pWarnings,  '\n', '\n   * ' );
                result = strrep(result, sprintf('\n'), '\n   * ' );
            end     
        end
        
%___________________________________________________________________________________________________ G: infoIndex
        function value = get.infoIndex(obj); value = obj.infoIndex; end
        
    end%GET/SET    
    
%===================================================================================================
    methods (Access = public) %                                                     P U B L I C  [M]

%___________________________________________________________________________________________________ preliminary
% Function to be called after initialization is complete but before the run has begun.
        function preliminary(obj)
            
            %--- Write abort file ---%
            if (labindex == 1 && obj.save.FSAVE)
                fid = fopen([obj.paths.save filesep obj.pAbortFile],'w');
                fprintf(fid,'1');
                fclose(fid);
            end
            labBarrier();
            
            %--- Preliminary setup for children managers ---%
            obj.save.preliminary();
            obj.image.preliminary();
            obj.parallel.preliminary();
            obj.fluid.preliminary();
            
            %--- Start code profiling if requested by ini ---%
            if obj.PROFILE
                if (obj.parallel.ACTIVE), mpiprofile('on'); else profile('on'); end
            end
        end
        
%___________________________________________________________________________________________________ initialize
% Run pre-simulation initialization actions that require already initialized initial conditions for
% the primary array objects.
        function initialize(obj, mass, mom, ener, mag, grav)
            obj.gravity.initialize(obj, mass, grav);
            obj.gravity.solvePotential(obj, mass, grav);
            
            obj.fluid.radiation.initialize(obj, mass, mom, ener, mag);
        end
        
%___________________________________________________________________________________________________ postliminary
% Function to be called at the end of a run to deactivate and cleanup remaining resources before
% ending the run entirely.
        function postliminary(obj)
            
            %--- Copy the log file to the results directory ---%
            try
                logFile = evalin('base','logFile;');
                if ~isempty(logFile)
                    copyfile(logFile,[obj.paths.save '/logfile.out']);
                end
            catch ME
            end
            
            %--- Stop and save code profiling if active ---%
            if obj.PROFILE
                if obj.parallel.ACTIVE
                    mpiprofile('off');
                    proInfo = mpiprofile('info');
                    
                    if labindex == 1
                        pInfo = cell(1,numlabs);
                        pInfo{1} = proInfo;
                        for i=2:numlabs
                            [labInfo, index] = labReceive('any',100);
                            pInfo{index} = labInfo;
                        end
                        save(strcat(obj.paths.save, filesep, 'profile'),'pInfo');
                    else
                        labSend(proInfo,1,100);
                    end
                    labBarrier();
                    obj.save.logPrint('Save barrier passed for %g\n',labindex);
                else
                    profile('off');
                    proInfo = profile('info');
                    save(strcat(obj.paths.save, filesep, 'profile'),'proInfo');
                end
                
            end
            
            obj.save.logPrint('Run complete.\n');
            obj.save.postliminary();
        end
    
%___________________________________________________________________________________________________ setGridSpacing
% Creates the grid cell spacing arrays necessary for uniform OR non-uniform grids. In uniform cases,
% the grid spacing parameter is a single value for each dimension. For non-uniform cases the dgrid
% for each dimension is a 3D array with appropriate spacing values for each grid cell.
% * dgrid        grid spacing structure (or 1x3 double array) from run file input            *
% * gridSize    number of grid cells in each spatial dimension                                int(3)
        function setGridSpacing(obj, dgrid, gridSize)
            
            %--- Uniform Spacing ---%
            if ~isstruct(dgrid)
                obj.MINDGRID = dgrid;
                for i=1:3    
                    obj.DGRID{i} = dgrid(i); 
                end
                return
            end
            
            %--- Non-uniform Spacing ---%
            obj.MINDGRID    = zeros(1,3);
            fields          = {'x','y','z'};    
            for i=1:3
                if isfield(dgrid,fields{i})
                    if isstruct(dgrid.(fields{i}))
                        dg = dgrid.(fields{i});
                        
                        if isfield(dg, 'value')
                            value = dg.value;
                        elseif isfield(dgrid, 'value')
                            value = dgrid.value;
                        else
                            value = 1;
                        end
                        
                        if isfield(dg, 'points')
                            interps = dg.points;
                        else
                            interps = [0,1;100,1];
                        end
                        
                        %--- Rescale Interpolation Points ---%
                        %       Interpolation points are entered during initialization as 
                        %       percentages so that they are not explicitly expressed in terms of a
                        %       grid size. However, before interpolation they must be converted to
                        %       grid points.
                        interps(:,1)    = floor((gridSize(i)-1)*interps(:,1)/100) + 1;
                        
                        %--- Initialize 1D Interploted Array ---%
                        lineSize        = ones(1,3);
                        lineSize(i)     = gridSize(i);
                        deltaGridLine   = zeros(lineSize);
                        
                        %--- Interpolate and Replicate ---%
                        %       Uses Piecewise Hermite interpolation to spline interpolate smooth
                        %       data grid points over the specified interpolation values and then
                        %       replicate the 1D vector result as the full 3D DGRID points array.
                        lineSize        = num2cell(lineSize);
                        lineSize{i}     = 1:gridSize(i);
                        deltaGridLine(lineSize{:}) ...
                                        = interp1(interps(:,1), value*interps(:,2), ...
                                          1:gridSize(i), 'pchip');
                                      
                        gridReps        = gridSize; 
                        gridReps(i)     = 1;
                        obj.DGRID{i}    = repmat(deltaGridLine, gridReps);
                    else
                        obj.DGRID{i}    = dgrid.(fields{i});
                    end
                else
                    obj.DGRID{i}        = dgrid.value;
                end
                obj.MINDGRID(i)         = minFinderND(obj.DGRID{i});
            end
        end
        
%___________________________________________________________________________________________________ addFades
        function addFades(obj,iniFades)
            if isempty(iniFades); return; end
            
            ids = {ENUM.MASS, ENUM.MOM, ENUM.ENER, ENUM.MAG, ENUM.GRAV};
            obj.fades = cell(1,length(iniFades));
            for i=1:length(iniFades)
                switch iniFades.type
                    case ENUM.POINT_FADE;      
                        obj.fades{i} = PointFade(obj.gridSize, iniFades.location, iniFades.size);
                end
                
                obj.fades{i}.fluxes     = iniFades.fluxes;
                obj.fades{i}.activeList = iniFades.active;
                for n=1:length(iniFades.active)
                    if ~any(strcmp(iniFades.active{n},ids))
                        warning('ImogenManager:Fade', 'Unable to resolve fade for %s.', ...
                                    iniFades.active{n});
                    end
                end
            end
            
            
            
        end
        
        
%___________________________________________________________________________________________________ appendInfo
% Appends an info string and value to the info cell array.
% * info    the information string                                                        srt
% * value    the value corresponding to the information string                            *
        function appendInfo(obj, info, value)
            
            %--- Resize info cell array on overflow ---%
            if (obj.infoIndex > size(obj.info,1))
                obj.info = [obj.info; cell(30,2)];
            end
            
            %--- Append information to info cell array ---%
            obj.info{obj.infoIndex,1} = info;
            obj.info{obj.infoIndex,2} = ImogenRecord.valueToString(value);
            obj.infoIndex = obj.infoIndex + 1;
        end
        
%___________________________________________________________________________________________________ initializeMode
% Parses a mode structure to set values for the code.
% * modeStruct        the structure containing boolean mode fields                        struct
        function initializeMode(obj, modeStruct)
            if ( isempty(modeStruct) || ~isstruct(modeStruct) ); return; end
            
            modes = {'fluid','magnet','gravity'};
            for i=1:length(modes)
                if isfield(modeStruct, modes{i})
                    obj.(modes{i}).ACTIVE = modeStruct.(modes{i});
                end
            end
            
        end
        
%___________________________________________________________________________________________________ appendWarning
% Appends a warning string (with value if necessary) to the warning string.
% * warning        the warning string                                                        str
% * (type)        the kind of warning ( >DEFAULT< or OVERRIDE)                            ENUM
% * (value)        numeric value related warning                                            *
        function appendWarning(obj, warning, type, value)
            if (nargin < 3 || isempty(type) ),      type = obj.DEFAULT; end
            if (nargin < 4 || isempty(value) ),     value = '';            end
            switch (type)
                case 'over',    typeStr = 'OVERRIDE';
                case 'def',     typeStr = 'DEFAULT';
            end
            newWarning = sprintf('\n%s: %s',typeStr,warning);
            obj.pWarnings = [obj.pWarnings, newWarning];
            if ~isempty(value), obj.pWarnings = [obj.pWarnings, sprintf('\t\t%s', ...
                                                ImogenRecord.valueToString(value))]; end
        end
        

%___________________________________________________________________________________________________ abortCheck
% Determines if the state.itf file abort bit has been set, and if so adjusts the time manager so
% that the active run will complete on the next iteration.
        function abortCheck(obj)
            %--- Initialization ---%
            if obj.time.iteration >= obj.time.ITERMAX;    return; end %Skip if already at run's end
            testTime = rem(now,1);
            
            %--- Run abort check ---%
            if (labindex == 1)
                if abs(testTime - obj.pAbortTime) > 0.01 %Check every ~15 minutes
                    fid = fopen([obj.paths.save filesep obj.pAbortFile],'r');
                    res = fscanf(fid, '%s', 1);
                    fclose(fid);
                    res = labBroadcast(1, res);
                else
                    res = labBroadcast(1, []);
                end
            else
                res = labBroadcast(1);
            end
            if ~isempty(res) &&    strcmp(res,'0'); obj.time.ITERMAX = obj.time.iteration + 1; end
            obj.pAbortTime = testTime;
        end
        
    end%PUBLIC
    
%===================================================================================================
    methods (Access = private) %                                                  P R I V A T E  [M]
        
%___________________________________________________________________________________________________ ImogenManager
% Creates a new ImogenManager object and initializes default values.
        function obj = ImogenManager()
            obj.bc          = BCManager.getInstance();          obj.bc.parent           = obj;
            obj.time        = TimeManager.getInstance();        obj.time.parent         = obj;
            obj.save        = SaveManager.getInstance();        obj.save.parent         = obj;
            obj.image       = ImageManager.getInstance();       obj.image.parent        = obj;
            obj.gravity     = GravityManager.getInstance();     obj.gravity.parent      = obj;
            obj.treadmill   = TreadmillManager.getInstance();   obj.treadmill.parent    = obj;
            obj.fluid       = FluidManager.getInstance();       obj.fluid.parent        = obj;
            obj.magnet      = MagnetManager.getInstance();      obj.magnet.parent       = obj;
            obj.parallel    = ParallelManager.getInstance();    obj.parallel.parent     = obj;
            obj.paths       = Paths();
            obj.info        = cell(30,2);
            obj.infoIndex   = 1;
            obj.DEBUG       = false;
            obj.PROFILE     = false;
            obj.pAbortTime  = rem(now,1);
            obj.pAbortFile  = '/state.itf';
            obj.matlab      = ver('matlab');
            obj.DGRID       = num2cell(ones(1,3));
            obj.MINDGRID    = ones(1,3);
            

            obj.useGPU      = false;
            obj.pureHydro   = false;
        end
        
    end%PRIVATE
    
%===================================================================================================    
    methods (Static = true) %                                                      S T A T I C   [M]
        
%___________________________________________________________________________________________________ getInstance
% Accesses the singleton instance of the ImogenManager class, or creates one if none have
% been initialized yet.
        function singleObj = getInstance()
            persistent instance;
            if isempty(instance) || ~isvalid(instance) 
                instance = ImogenManager(); 
            end
            singleObj = instance;
        end
      
    end%STATIC
    
    
end %CLASS
