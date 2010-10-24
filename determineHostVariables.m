function [host, imogenRootPath, resultPath] = determineHostVariables()
% This function identifies the computer system and user that are running imogen and returns the
% appropriate directories for the root imogen code and the results folder. If you are not a
% registered user the folders are defaulted and information regarding your system is returned. This
% can help you by supplying the values needed for the username and hostname as you create your new
% user file.
%
%<< host              name of the computer hosting imogen                             str
%<< imogenRootPath    path to the imogen.m and includeImogenPaths.m files             str
%<< resultPath        path to the root results storage folder                         str

    %--- Determine Host ---%
    try [ret1, host] = system('hostname');
    catch MErr, fprintf('Warning: hostname command unavailable.\n');
    end
    if (ret1 ~= 0)
        if ispc(); host = getenv('COMPUTERNAME');
        else       host = getenv('HOSTNAME');      
        end
    end
    host = strtrim(lower(host));
    
    %--- Determine User ---%
    try [ret2, user] = system('whoami');
    catch MErr, fprintf('Warning: User account unknown.\n');
    end
    user = strtrim(user);
    
    %--- Find user match ---%
    registeredUsers = addFileUsers();
    userFound       = false;
    for i=1:size(registeredUsers,1)
        if ~isempty(regexpi(host, registeredUsers{i,2})) ...
            && (strcmpi(user, registeredUsers{i,1}) ||  isempty(registeredUsers{i,1}) )
        
            imogenRootPath  = registeredUsers{i,3};
            resultPath      = registeredUsers{i,4};
            userFound       = true; 
            break;
        end
    end
    
    %-----------------------------------------------------------------------------------------------
    % If no user found: Set Defaults & Verbose 
    %-----------------------------------------
    if ~userFound
        if (isunix() || ismac())
            imogenRootPath = '~/imogen';
            resultPath = '~/Results';
        else
            imogenRootPath = 'C:\imogen';
            resultPath = 'C:\Results';
        end
        fprintf('Unable to ascertain host and user id. Environmental paths may be incorrect.\n');
        fprintf('Imogen registered the following values for your system:\nHost: %s\nUser: %s', ...
                    host, user);
    end    
end
