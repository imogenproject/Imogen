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
    version = cell2mat( textscan(fid, '%n', 1) );
    cellStr = textscan(fid, '%s %s', 1);
    month   = cell2mat(cellStr{1,1});
    year    = cell2mat(cellStr{1,2});
    fclose(fid);

    try 
        [status, revision] = system('svnversion');
        if ~status
            revision = deblank(regexpi(revision,'[0-9]*:','split'));
            for i=1:length(revision), if ~isempty(revision{i}), revision = revision{i}; break; end; end
        else
            revision = 'unknown';
        end
    catch MErr, revision = 'unknown';
    end
    
    detailedVersion = sprintf('%g.%s',version,revision);
    
    header     = strcat('IMOGEN v', detailedVersion);
    headerLine = '===========================================================';
    index      = floor(0.5 * (length(headerLine) - (length(header) + 2)));
    headerLine(index:(index + length(header)+1)) = [' ' header ' '];
    
    fprintf('\n\n%s\n',headerLine);
    fprintf('   Updated: %s %s       Created: March, 2007               \n', month, year);
    fprintf('===========================================================\n');
end
