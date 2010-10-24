function printDEBUG(str, mode)
% This routine is designed to display str information on the command line or logfile for help when
% debugging the code.
%
%>> str       primary information string to print								  str     *
%>> mode      boolean debug mode                                                  bool    t/f

    
    %--- Display only if debug mode is TRUE ---%
    if (mode); 
        
        %--- Get stack info to locate current line in code ---%
        [stack,I] = dbstack; 
        if (length(stack) > 1)
            strLen = 30 - length(str); if (strLen <= 0); strLen = 4; end
            str = strcat(str, sprintf('%s[', blanks(strLen)));
            for i=2:max( (length(stack)-1), 2)
                str = strcat(str, sprintf( ' (L%s:%s) ',num2str(stack(i).line),deblank(stack(i).name) )); 
            end
            str = strcat(str, ' ]');
        end
        
        %--- Output results ---%
        disp(str); 
    end
    
end