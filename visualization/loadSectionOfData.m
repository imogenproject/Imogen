function loadSectionOfData(filterExpression, startVec, endVec)
% This function loads subsections of an Imogen data set based on starting and ending times. If no
% start and end date-time vectors are specified, all data will be loaded. The date vectors use the
% standard Matlab format syntax which is:
%
%   Format:  [year month day hour minute second]
%   Example: [2009 11    17  23   15     0]
%
% In this function you can omit the year at the beginning and the function will assume that the
% data you want to load was created during the current year.
%
% The function will also accept date vectors with missing elements. The remaining elements are
% assumed to be zero. For example, if you wanted to load all the files on a specific day:
% 
%   loadSectionofData('3D*.mat', [11 17], [11 18]);
%
% which would load all files created between 00:00:00 on 11/17/2009 and 00:00:00 on 11/18/2009.
%
%
%>> filterExpression    Regular expression used to filter the data.                     string
%>> startVec            Date vector specifying earliest allowed time to load.           datevec
%>> endVec              Date vector specifying latest allowed time to load.             datevec

    if nargin < 1 || isempty(filterExpression)
        filterExpression = '3D*.mat';
        fprinft('\nWARNING: No filter expression supplied. Assuming filter is: 3D*.mat.\n');
    end
    
    if nargin < 2 || isempty(startVec)
        startVec = zeros(1,6);
        fprintf('\nWARNING: No starting date vector supplied. Assuming start is the beginning of time.\n');
    end
    if (startVec(1) < 2000), startVec = [year(now) startVec]; end
    if length(startVec) < 6
        startVec = [startVec zeros(1,6-length(startVec))];
    end
    
    if nargin < 3 || isempty(endVec)
        endVec = datevec(now);
        fprintf('\nWARNING: No ending date vector supplied. Assuming date as recent as now.\n');
    end
    if (endVec(1) < 2000), endVec = [year(now) endVec]; end
    if length(endVec) < 6
        endVec = [endVec zeros(1,6-length(endVec))];
    end
    
    files    = dir(filterExpression);

    for i=1:length(files)
        if (files(i).datenum >= datenum(startVec) && files(i).datenum <= datenum(endVec))
            fprintf('\nLoading: %s . . . . . . .',files(i).name);
            S = load(files(i).name);
            fields = fieldnames(S);
            for j=1:length(fields)
                assignin('base',fields{j},S.(fields{j}));
            end
            fprintf('complete');
        else
            fprintf('\nIgnored: %s',files(i).name);
        end
    end

    fprintf('\nOperation complete.\n');
end
