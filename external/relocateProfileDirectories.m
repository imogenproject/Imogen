function object = relocateProfileDirectories(object, oldPath, newPath)
% When profiles are created on deployment servers, the directory listings for all of the m files
% will be incorrect for viewing on a data analysis workstation, causing errors in viewing the 
% saved profile statistics. This function recursively searches and replaces the server root imogen 
% path and replaces it with the one on your workstation as specified in the arguments.
%
% USAGE: object = relocateProfileDirectories(proIfno, '/home/user/imogen', 'C:Research/imogen');
%
%>< object      Profile object to be corrected.                                         struct|cell
%>> oldPath     Root imogen path on the server to be replaced in all entries.           str
%>> newPath     Root imogen path on your local workstation to be inserted.              str

    if isstruct(object)
        for n=1:length(object);
            fields = fieldnames(object(n));
            for i=1:length(fields)
                if isstruct(object(n).(fields{i}))
                    object(n).(fields{i}) = relocateProfileDirectories(object(n).(fields{i}), oldPath, newPath);
                elseif ischar(object(n).(fields{i}))
                    object(n).(fields{i}) = strrep(object(n).(fields{i}), oldPath, newPath);
                end
            end
        end
    elseif iscell(object)
        for n=1:length(object)
            if isstruct(object(n))
                object(n) = relocateProfileDirectories(object(n), oldPath, newPath);
            elseif ischar(object(n))
                object(n) = strrep(object(n), oldPath, newPath);
            end
        end
    elseif ischar(object)
        object = strrep(object, oldPath, newPath);
    end
end
