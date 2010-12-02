function [version, detailedVersion] = versionInfo()
% This routine is used for showing the command line UI display at the beginning of each run. It
% grabs the version information from the version.txt file and incorporates it into the header
% information.
%
%<< version         value of the version as read from the version.txt file      double  #

    %-----------------------------------------------------------------------------------------------
    % Read the version file
    %----------------------
    fid     = fopen('version.txt');
    version = deblank(fgetl(fid));
    modDate = fgetl(fid);
    fclose(fid);
   
    try
        [errorCode, status] = system('git status');
        if (errorCode ~= 0)
            detailedVersion = [version '-???'];
        else
            % Test for master branch
            if isempty(strfind(status, 'branch master'))
                branch = '?';
            else
                branch = 'MB';
            end
    
            pieces   = regexpi(status, sprintf('\n'), 'split');
            modified = false;
            
            for i=1:length(pieces)
                item = pieces{i};
                if length(item) > 1
                    modified = (modified || ~isempty(strfind(pieces{i}, 'modified')));
                end
            end
        
            if modified
                modified = 'MOD';
            else
                modified = 'BASE';
            end
            
            detailedVersion = [version, '.', branch, '-', modified];
        end
    catch e
        detailedVersion = [version '.???'];
    end

    header     = strcat('IMOGEN v', detailedVersion);
    headerLine = '===========================================================';
    index      = floor(0.5 * (length(headerLine) - (length(header) + 2)));
    headerLine(index:(index + length(header)+1)) = [' ' header ' '];
    
    fprintf('\n\n%s\n',headerLine);
    fprintf('   Updated: %s       Created: March, 2007               \n', modDate);
    fprintf('===========================================================\n');
end
