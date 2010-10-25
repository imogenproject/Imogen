function imogenLoad(runFile, logFile, alias)
% This script is the command and control manager for command line run Imogen scripts. It is
% designed to simplify the run script syntax and allow for greater extensibility.
% 
%>> runFile		run file function name to execute								str
%>> logFile		log file name for writing output information					str
%>> uid         unique identifier string for the run                            str

    %-- Initialize Imogen directory ---%
    starterRun();
    runFile = strrep(runFile,'.m','');
    assignin('base','logFile',logFile);
    assignin('base','alias',alias);
    try
        eval(runFile);
    catch ME
       rethrow(ME);
    end
    enderRun();
end