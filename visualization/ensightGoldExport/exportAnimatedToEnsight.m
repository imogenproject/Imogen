function exportAnimatedToEnsight(outBasename, inBasename, padlength, range, timeNormalization)
%>> outBasename:       Base filename for output Ensight files
%>> inBasename:        Input filename for Imogen .mat savefiles
%>> padlength:         Number of zeros in Imogen filenames
%>> range:             Set of .mats to export
%>> timeNormalization: Allows Imogen timestep-time to be converted into characteristic time units
	
%--- Interactively fill in missing arguments ---%
if nargin < 4
    fprintf('Not enough input arguments to run automatically.\n');
    outBasename = input('Base filename for exported files (e.g. "torus1"): ', 's');
    inBasename  = input('Base filename for source files, (e.g. "3D_XYZ", no trailing _):','s');
    padlength   = input('Length of frame #s in source files (e.g. 3D_XYZ_xxxx -> 4): ');
    range       = input('Range of frames to export; _START = 0 (e.g. 0:50:1000 to do every 50th frame from start to 1000): ');
    timeNormalization = input('Characteristic time to normalize by (e.g. alfven crossing time or characteristic rotation period. If in doubt hit enter): ');
    if timeNormalization == 0; timeNormalization = 1; end;
end

pertonly = input('Export perturbed quantities (1) or full (0)? ');

%--- Initialization ---%
fprintf('Beginning export of %i files\n', numel(range));
exportedFrameNumber = 0;

if max(round(range) - range) ~= 0; error('ERROR: Frame range is not integer-valued.\n'); end
if min(range) < 0; error('ERROR: Frame range must be nonnegative.\n'); end

range = removeNonexistantEntries(inBasename, padlength, range);
maxFrameno = max(range);

if nargin == 4; timeNormalization = 1; end;

equilframe = [];

%--- Loop over given frame range ---%
for ITER = 1:numel(range)
    % Take first guess; Always replace _START
    fname = sprintf('%s_%0*i.mat', inBasename, padlength, range(ITER));
    if range(ITER) == 0; fname = sprintf('%s_START.mat', inBasename); end

    % Check existance; if fails, try _FINAL then give up
    if exist(fname, 'file') == 0
        fname = sprintf('%s_FINAL.mat', inBasename);
        if exist(fname, 'file') == 0
            % Weird shit is going on. Run away.
            error('UNRECOVERABLE: File existed when checked but is not openable.\n');
        end
    end
    
    fprintf('Exporting %s as frame %i... ', fname, exportedFrameNumber);
    load(fname);

    structName = who('sx_*');
    structName = structName{1};

    if ITER == 1
        eval(sprintf('equilframe = %s;', structName));
    end

    eval(sprintf('dataframe = %s;', structName));
    clear -regexp 'sx_';
    
    if pertonly == 1
        dataframe = subtractEquil(dataframe, equilframe);
    end

    writeEnsightDatafiles(outBasename, exportedFrameNumber, dataframe);
    if range(ITER) == maxFrameno
        writeEnsightMasterFiles(outBasename, range, dataframe, timeNormalization);
    end

    exportedFrameNumber = exportedFrameNumber + 1;
    fprintf('done.\n');
end

end

function out = subtractEquil(in, eq)
out = in;

out.mass = in.mass - eq.mass;
out.ener = in.ener - eq.ener;

out.momX = in.momX - eq.momX;
out.momY = in.momY - eq.momY;
out.momZ = in.momZ - eq.momZ;

out.magX = in.magX - eq.magX;
out.magY = in.magY - eq.magY;
out.magZ = in.magZ - eq.magZ;

end


function newrange = removeNonexistantEntries(inBasename, padlength, range)

existrange = [];

for ITER = 1:numel(range)
    % Take first guess; Always replace _START
    fname = sprintf('%s_%0*i.mat', inBasename, padlength, range(ITER));
    if range(ITER) == 0; fname = sprintf('%s_START.mat', inBasename); end

    % Check existance; if fails, try _FINAL then give up
    doesExist = exist(fname, 'file');
    if doesExist == 0;
        fname = sprintf('%s_FINAL.mat', inBasename);
        doesExist = exist(fname, 'file');
    end
    
    if doesExist ~= 0; existrange(end+1) = ITER; end
end

newrange = range(existrange);
if numel(newrange) ~= numel(range);
    fprintf('WARNING: Removed %i entries that could not be opened from list.\n', numel(range)-numel(newrange));
end

if numel(newrange) == 0;
   error('UNRECOVERABLE: No files indicated existed. Perhaps remove trailing _ from base name?\n'); 
end

end
