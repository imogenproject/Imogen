function userRegistry = addFileUsers()
% This function iterates over all of the files listed in Imogen's users directory and adds them to
% the registerd users list.
%
%<< userRegistry     Cell array of registered users.                             cell(N,4)

    pathToThisFile = fileparts(mfilename('fullpath'));  % Get poth of this file.
    files          = dir(pathToThisFile);               % Get listing of all files in directory.
    numFiles       = length(files) - 1;
    users          = cell(numFiles,1);
    numEntries     = 0;
    index          = 1;
    
    %--- Populate users cell array from files ---%
    for i=1:length(files)
        if ~isempty(regexpi(files(i).name,'user')) && ~strcmpi(files(i).name,'addFileUsers.m')
            users{index,1} = eval(strrep(files(i).name,'.m','()'));
            numEntries     = numEntries + size(users{index,1},1);
            index          = index + 1;
        end
    end

    %--- Create userRegistry output from users ---%
    userRegistry = cell(numEntries,4);
    index        = 1;
    
    for i=1:size(users,1)
        entry = users{i,1};
        for j=1:size(entry,1)
            for m = 1:4,    userRegistry{index,m} = entry{j,m}; end
            index = index + 1;
        end
    end
    
end