function appendAnimatedToEnsight(outBasename, inBasename, padlength, oldrange, addrange, timeNormalization)
%>> outBasename:       Base filename for output Ensight files
%>> inBasename:        Input filename for Imogen .mat savefiles
%>> padlength:         Number of zeros in Imogen filenames
%>> range:             Set of .mats to export
%>> timeNormalization: Allows Imogen timestep-time to be converted into characteristic time units
	
%--- Initialization ---%
fprintf('Beginning export of %i files\n', numel(addrange));
exportedFrameNumber = numel(oldrange);

if max(round(addrange) - addrange) ~= 0; error('ERROR: Frame range is not integer-valued.\n'); end
if min(addrange) < 0; error('ERROR: Frame range must be nonnegative.\n'); end

addrange = removeNonexistantEntries(inBasename, padlength, addrange);
maxFrameno = max(addrange);

if nargin == 4; timeNormalization = 1; end;

%--- Loop over given frame range ---%
for ITER = 1:numel(addrange)
    % Take first guess; Always replace _START
    fname = sprintf('%s_%0*i.mat', inBasename, padlength, addrange(ITER));
    if addrange(ITER) == 0; fname = sprintf('%s_START.mat', inBasename); end

    % Check existance; if fails, try _FINAL then give up
    if exist(fname, 'file') == 0
        fname = sprintf('%s_FINAL.mat', inBasename);
        if exist(fname, 'file') == 0
            % Weird shit is going on. Run away.
            error('UNRECOVERABLE: File previously checked out but now does not exist.\n');
        end
    end
    
    fprintf('Exporting %s as frame %i... ', fname, exportedFrameNumber);
    load(fname);

    structName = who('sx_*');
    structName = structName{1};

    eval(sprintf('dataframe = %s;', structName));
    clear -regexp 'sx_';
    
    writeEnsightDatafiles(outBasename, exportedFrameNumber, dataframe) 
    if addrange(ITER) == maxFrameno
        writeEnsightMasterFiles(outBasename, [oldrange addrange], dataframe, timeNormalization);
    end

    exportedFrameNumber = exportedFrameNumber + 1;
    fprintf('done.\n');
end

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
