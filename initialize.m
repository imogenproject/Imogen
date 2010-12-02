function run = initialize(ini)
% This routine runs through all of the initialization variables testing to see if they've been
% set by user input. If not, then they are defaulted with a message to warn that a default
% value has been set. The values are parsed from the ini structure and then passed back as
% separate variables for explicit use in the code.
%
%>> ini          structure containing run related variable information              struct
%<< run          initialized results                                                ImogenManager  H
    
%% ============================= BEGIN RUNVAL PROPERTY INITIALIZATION =========================== %
   
    %--- Clear all manager singletons from possible previous runs ---%
    clear('ImogenManager','SaveManager','TimeManager','ImageManager','GravityManager', ...
            'FluidManager', 'MagnetManager', 'ParallelManager', 'BCManager', 'TreadmillManager');
    
    fclose all; % Prevent any lingering saves from disrupting run.
    
    run             = ImogenManager.getInstance();
    run.gridSize    = ini.grid;
    [run.version, run.detailedVersion]                        = versionInfo();
    [run.paths.hostName, run.paths.imogen, run.paths.results] = determineHostVariables();
        
%% .bcMode                      Edge condition modes

try
    if isa(ini.bcMode, 'struct')
        run.bc.modes = ini.bcMode;
        run.appendInfo('BC Mode', run.bc.modes);
    elseif isa(ini.bcMode, 'char')
        modes.x = ini.bcMode; modes.y = ini.bcMode; modes.z = ini.bcMode;
        run.bc.modes = modes;
    else
        error(['BoundaryConditionError: Boundary condition field of type %s is not recognized.' ... 
            ' bcMode recognizes string or structure input. Run aborted.'],class(ini.bcMode));
    end
    run.appendInfo('Boundary Conditions', run.bc.modes);
catch MERR, loc_initializationError('bcMode',MERR);
end

%% .bcInfinity                  # of cells for infinity of transparent edge condition

try
    run.bc.infinity = ini.bcInfinity;    
    run.appendInfo('Infinity Distance', run.bc.infinity);
catch MERR, loc_initializationError('infinity',MERR);
end

%% .fluxLimiter                 # type of flux limiter to use globally or in each dimension.

try
    run.fluid.setFluxLimiters(ini.fluxLimiter);
    run.magnet.limiter = run.fluid.limiter;
    run.appendInfo('Flux Limiters', run.fluid.limiter);
catch MERR, loc_initializationError('infinity',MERR);
end

%% .cfl                         CFL prefactor 

try
    run.time.CFL = ini.cfl;
    run.appendInfo('CFL Prefactor', run.time.CFL);
catch MERR, loc_initializationError('cfl',MERR);
end

%% .timeUpdateMode              How frequently the timestep should be updated.

try
    run.time.updateMode = ini.timeUpdateMode;
    run.appendInfo('Time Update Mode', run.time.updateMode);
catch MERR, loc_initializationError('timeUpdateMode',MERR);
end

%% .thresholdMass               Threshold value below which gravity will not act

try
    run.fluid.MASS_THRESHOLD = ini.thresholdMass;
    run.appendInfo('Mass Threshold', run.fluid.MASS_THRESHOLD);
catch MERR, loc_initializationError('thresholdmass',MERR);
end

%% .minMass                     Minimum allowed mass density value

try 
    run.fluid.MINMASS = ini.minMass;
    run.appendInfo('Minimum mass density', run.fluid.MINMASS);
catch MERR, loc_initializationError('minmass',MERR);
end

%% .profile                     Enable the profiler to record execution information

try 
    run.PROFILE = ini.profile;
    if (run.PROFILE); run.appendWarning('MATLAB profiler will be active for this run.'); end
    run.appendInfo('Profiler', run.PROFILE);
catch MERR, loc_initializationError('profile',MERR);
end

%% .iterMax                     Maximum number of iterations

try
    run.time.ITERMAX = ini.iterMax; 
    run.appendInfo('Maximum iteration',run.time.ITERMAX);
catch MERR, loc_initializationError('iterMax',MERR);
end

%% .timeMax                     Maximum simulation time

try
    run.time.TIMEMAX = ini.timeMax;
    run.appendInfo('Maximum simulation time.',run.time.TIMEMAX);
catch MERR, loc_initializationError('timeMax',MERR);
end

%% .wallMax                     Maximum wall time used for the run

try
    run.time.WALLMAX = ini.wallMax;
    run.appendInfo('Maximum allowed wall time (hours).',run.time.WALLMAX);
catch MERR, loc_initializationError('wallTimeMax',MERR);
end

%% .runCode                     Run code for the simulation

try
    run.paths.runCode = ini.runCode;
catch MERR, loc_initializationError('runCode',MERR);
end

%% .alias                       Alias for the simulation

try
    run.paths.alias = ini.alias;
catch MERR, loc_initializationError('alias',MERR);
end

%% .info                        Run information String

try
    run.about = ini.info;
catch MERR, loc_initializationError('info',MERR);
end

%% .notes                       Add notes (user generated)

try
    run.notes = ini.notes;    
catch MERR, loc_initializationError('notes',MERR);
end 

%% .iniInfo                     Add initialization information (procedurally generated)

try
    run.iniInfo = ini.iniInfo;     
catch MERR, loc_initializationError('iniinfo',MERR);
end 

%% .mode                        Activate code mode types (e.g. fluid, magnet, etc...)

try
    run.initializeMode(ini.mode);
catch MERR, loc_initializationError('mode',MERR);
end    

%% .dGrid                       Grid cell dimensions

try
    if isempty(ini.dGrid)
        run.setGridSpacing([1 1 1]);
        run.appendWarning('No dGrid parameter specified. Defaulting to isotropic unity values.');
    elseif isa(ini.dGrid,'double')
        len = length(ini.dGrid);
        switch (len)
            case 1;     run.setGridSpacing(ini.dGrid * ones(1,3), run.gridSize);
            case 2;     run.setGridSpacing([ini.dGrid, max(ini.dgrid)], run.gridSize);
            case 3;     run.setGridSpacing(ini.dGrid, run.gridSize);
            otherwise; error('Unable to parse ini.dGrid value.');
        end
    else
        run.setGridSpacing(ini.dGrid, run.gridSize);
    end
    run.appendInfo('Grid spacing', ImogenRecord.valueToString(run.DGRID));
catch MERR, loc_initializationError('dgrid',MERR);
end

%% .gamma                       Polytropic index for equation of state

try
    run.GAMMA = ini.gamma;
    run.appendInfo('Gamma', run.GAMMA);
catch MERR, loc_initializationError('gamma',MERR);
end

%% .debug                       Run the code in debug mode 

try
    run.DEBUG = ini.debug;
    if (run.DEBUG), run.appendWarning('Running in debug mode.'); end
catch MERR, loc_initializationError('debug',MERR);
end

%% .save                        Save data to files

try
    run.save.FSAVE = logical(ini.save);
catch MERR, loc_initializationError('save',MERR);
end

%% .ppSave                      Percentage executed between saves

try
    run.save.PERSLICE(1) = ini.ppSave.dim1;
    run.save.PERSLICE(2) = ini.ppSave.dim2;
    run.save.PERSLICE(3) = ini.ppSave.dim3;
    run.save.PERSLICE(4) = ini.ppSave.cust;
catch MERR, loc_initializationError('ppSave',MERR);
end

%% .slice                       Index Locations for slice and image save files

try       
    run.save.SLICEINDEX = ini.slice;
    run.appendInfo('Slices will be saved at',run.save.SLICEINDEX);
catch MERR, loc_initializationError('slice',MERR);
end

%% .activeSlices                Which data slices to save

try
    slLabels = {'x','y','z','xy','xz','yz','xyz','cust'};
    for i=1:8
        if ~isfield(ini.activeSlices,slLabels{i}); run.save.ACTIVE(i) = false;
        else run.save.ACTIVE(i) = logical(ini.activeSlices.(slLabels{i})); 
        end
        if run.save.ACTIVE(i)
            run.appendInfo('Saving slice', upper(slLabels{i}));
        end
    end
catch MERR, loc_initializationError('activeSlices',MERR);
end

%% .customSave                  Custom saving properties

try 
    saveStr = '''slTime'',''slAbout'',''version'',''slGamma'',''sldGrid''';
        
    custom = ini.customSave;
    if isstruct(custom)
        if (isfield(custom,'mass')  && custom.mass),    saveStr = [saveStr ',''slMass''']; end
        if (isfield(custom,'mom')   && custom.mom),     saveStr = [saveStr ',''slMom'''];  end 
        if (isfield(custom,'ener')  && custom.ener),    saveStr = [saveStr ',''slEner''']; end
        if (isfield(custom,'mag')   && custom.mag),     saveStr = [saveStr ',''slMag'''];  end
        run.save.customSaveStr = saveStr;
        customStr = 'Active';
    else
        run.save.customSaveStr = '';
        customStr = 'Inactive';
    end
    run.appendInfo('Custom save', customStr);
catch MERR, loc_initializationError('customSave',MERR);
end

%% .specSaves                   Specific save iterations

try
    if ~isempty(ini.specSaves) && isa(ini.specSaves,'double')
        run.save.specialSaves3D = ini.specSaves;
        run.appendInfo('Special save points 3D', run.save.specialSaves3D);
    else
        if isfield(ini.specSaves,'dim1')
            run.save.specialSaves1D = ini.specSaves.dim1; 
            run.appendInfo('Special save points 1D', run.save.specialSaves1D);
        end

        if isfield(ini.specSaves,'dim2')
            run.save.specialSaves2D = ini.specSaves.dim2;
            run.appendInfo('Special save points 2D', run.save.specialSaves2D);
        end

        if isfield(ini.specSaves,'dim3')
            run.save.specialSaves3D = ini.specSaves.dim3;
            run.appendInfo('Special save points 3D', run.save.specialSaves3D);
        end
    end
catch MERR, loc_initializationError('specSaves',MERR);
end

%% .image                       Image saving properties

try

    fields = ImageManager.FIELDS;
    for i=1:length(fields)
        if isfield(ini.image,fields{i})
            run.image.(fields{i}) = ini.image.(fields{i});
        end
        if isfield(ini.image,'logarithmic') && isfield(ini.image.logarithmic,fields{i})
            run.image.logarithmic.(fields{i}) = ini.image.logarithmic.(fields{i});
        end
    end
    run.image.activate();
    
    if run.image.ACTIVE
    
        if isfield(ini.image,'interval')
            run.image.INTERVAL = max(1,ini.image.interval); 
        else
            run.image.INTERVAL = 1;
            run.appendWarning('Image saving interval set to every step.');
        end

        if isfield(ini.image,'colordepth');    colordepth = ini.image.colordepth;
        else                                    colordepth = 256;
        end

        if isfield(ini.image,'colormap'); run.image.createColormap(ini.image.colormap, colordepth);
        else                                  run.image.createColormap('jet',colordepth); 
        end
        
        imageSaveState = 'Active';
    else imageSaveState = 'Inactive';
    end
    run.appendInfo('Image saving is', imageSaveState);
catch MERR, loc_initializationError('image',MERR);
end

%% .treadmill                   Grid treadmill action

try
    if ini.treadmill > 0
        run.treadmill.ACTIVE = true;
        run.treadmill.DIRECTION = ini.treadmill;
        run.appendInfo('Treadmill', sprintf('Active along %g',run.treadmill.DIRECTION));
    end
catch MERR, loc_initializationError('treadmill',MERR);
end
    
%% .fades                       Fade objects

try
    run.addFades(ini.fades);
catch MERR, loc_initializationError('fades',MERR);
end

%% .gravity.fixedPotential      Carries over any fixed
try
    if isfield(ini.gravity,'fixedPotential')
        run.gravity.fixedPotential = ini.gravity.fixedPotential;
    else
        run.gravity.fixedPotential = 0;
    end

    gpsize = size(run.gravity.fixedPotential);
    if gpsize(1) > 1; run.appendInfo('Including fixed gravitational potential',1); end;

catch MERR, loc_initializationError('gravity.vars',MERR);
end

%% .gravity.CONSTANT
try
    if isfield(ini.gravity,'constant')
        run.gravity.constant = ini.gravity.constant;
    else
        run.gravity.constant = 1;
    end

    if isfield(ini.gravity,'iterMax')
        run.gravity.iterMax = ini.gravity.iterMax;
    else
        run.gravity.iterMax = 100;
    end;

    if isfield(ini.gravity,'tolerance')
        run.gravity.tolerance = ini.gravity.tolerance;
    else
        run.gravity.tolerance = 1e-6;
    end;

        if isfield(ini.gravity,'bconditionSource')
        run.gravity.bconditionSource = ini.gravity.bconditionSource;
    else
        run.gravity.bconditionSource = GRAV_BCSOURCE_FULL;
    end

        run.gravity.mirrorZ = ini.gravity.mirrorZ;

catch MERR, loc_initializationError('gravity.constant',MERR);
end

%% .gravity.solver                Gravity solver
try   
    run.gravity.setSolver(ini.gravity.solver);
    run.appendInfo('Gravity solver', run.gravity.TYPE);
catch MERR, loc_initializationError('gravity.solver',MERR);
end

%% .viscosity                   Artificial viscosity settings
try
    run.fluid.viscosity.type                      = ini.viscosity.type;
    run.fluid.viscosity.linearViscousStrength     = ini.viscosity.linear;
    run.fluid.viscosity.quadraticViscousStrength  = ini.viscosity.quadratic;
catch MERR, loc_initializationError('viscosity', MERR);
end

%% .radiation                   Radiation settings
try
   run.fluid.radiation.type                      = ini.radiation.type;
   run.fluid.radiation.exponent                  = ini.radiation.exponent;
   run.fluid.radiation.initialMaximum            = ini.radiation.initialMaximum;
    
catch MERR, loc_initializationError('radiation', MERR);
end


end

function loc_initializationError(property, caughtError)
% Handles errors thrown by the try statements in the initialization routine.
%
%>> property               the property that threw the error                        str
%>> caughtError            The error captured by a try block                        error


    fprintf('\n\n--- Unable to parse property %s. Run aborted. ---\n', property);
    rethrow(caughtError);
end
