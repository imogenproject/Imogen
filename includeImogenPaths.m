function includeImogenPaths()
% Load the paths relative to the Imogen root path into the matlab path variable, giving access to
% all of the routines without having to change the current directory procedurally.

    %--- Add default directories ---%
    %       Iterates over the list of directories provided below and adds them and all their
    %       subdirectories to the Matlab environmental path to access standard Imogen functionality.
    directories = {'classes','equilibrium','experiment','fluid','enums','utils','gravity', ...
                   'io','magnet','parallel','pic','pressure','run','save','shifting', ...
                   'testing','time','treadmill','users','external','visualization', ...
                   'flux','gpuclass'};

    pathToHere = fileparts(mfilename('fullpath'));
    addpath(pathToHere); % Add base Imogen directory without subdirectories.
    
    for i=1:length(directories); 
        addpath(genpath([pathToHere filesep directories{i}])); 
    end
    fprintf('Imogen directories have been added to the path listing.\n');

    
    
    %--- Add compatibility directories ---%
    %       To deal with incompatibilities created by multiple Matlab versions, as well as missing
    %       toolboxes or related functionality, imogen includes a compat directory. Each 
    %       subdirectory in this path stores contains the compatibility functions necessary to 
    %       restore compatibility for Imogen in non-standard envrionmental situations. Add these
    %       paths only as needed.
    compat      = [filesep 'compat' filesep];
    compWarn    = 'Imogen:Compatibility';
    
    %--- Matlab version ---%
    version     = ver('matlab');
    verYear     = str2double( version.Release(3:6) );
	verLabel    = version.Release(7);
    if (verYear == 2009 && strcmp(verLabel,'a'));
        addpath(genpath([pathToHere compat 'r2009a']))
        warning(compWarn,['Parallel functionality in your version is limited. ', ...
                          'Compatibility functions have been added to facilitate.']);
    elseif (verYear < 2010 && ~strcmp(verLabel,'a'))
        addpath(genpath([pathToHere compat 'r2009a']));
        warning(compWarn,['Imogen uses the functionality of the latest version of '...
                          'Matlab. Your version may not work properly as result']);
    end
    
    
    %--- Parallel computing toolbox ---%
    version = ver('distcomp');
    if isempty(version)
        addpath(genpath([pathToHere compat 'distcomp']));
        warning(compWarn,['Unable to detect a parallel computing toolbox installtion. Imogen ', ...
                          'must be run in serial mode. Compatibility functions have been added ',...
                          'to compensate.']);
    end
    
    
end
