function notes = warnDefault(notes,str,type)
% This routine prints out warnings to prompt and adds them to the notes variable for later
% printing to a file for record keeping.
%
%>> notes     existing notes to concanenate with                                      str     *
%>> str       string describing the warning                                           str     *
%>> type      type of warning 'over' or 'default'                                     str     *
%<< notes     updated notes string                                                    str     *


    if (nargin < 3); type = ''; end
    switch (type)
        case 'over'
            str = sprintf('\tOVERRIDE: %s',str);
        otherwise
            str = sprintf('\tDEFAULT: %s',str);
    end
    disp(str);
    notes = strcat(notes,str);
    
end