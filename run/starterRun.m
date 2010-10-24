function starterRun()
% This routine does the initialization and directory search to locate imogen root for starting a new
% run.

%---------------------------------------------------------------------------------------------------
% Locate Imogen root directory and source subfolder paths
%--------------------------------------------------------

	clear; % Clear existing variables from workspace to prevent conflict.

    % Attempt locate path from location of this m-file
    try
        rootPath = [strrep(mfilename('fullpath'),'starterRun','') '../'];
        cd(rootPath);
    catch MERR
    end
        
    found = false;
    for i=1:3
        
        files = dir(cd);
        for n=1:length(files)
            if strcmp(files(n).name,'imogen.m'), found = true; break; end
        end
        
        if found, break;
        else cd('..');
        end
    end

    try
        includeImogenPaths();    %Initialize directory structure
    catch MERR
        error('Imogen:starterRun:ImogenRootFailure','Unable to find Imogen root directory. Run aborted.');
    end
end
