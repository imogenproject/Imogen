function rootPath = IMOGENROOT()
% Returns the root path to the Imogen code for accessing functionality with relative paths.

    pathToHere = fileparts(mfilename('fullpath'));
    rootPath   = fileparts(pathToHere);
    rootPath   = [rootPath filesep]; 

end