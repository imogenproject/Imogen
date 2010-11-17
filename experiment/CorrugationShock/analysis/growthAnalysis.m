function [grorate phaserate groR phaseR timeVals fftVals] = growthAnalysis(inBasename, padlength, range, timeNormalization)
%>> inBasename:        Input filename for Imogen .mat savefiles
%>> padlength:         Number of zeros in Imogen filenames
%>> range:             Set of .mats to export
%>> timeNormalization: Allows Imogen timestep-time to be converted into characteristic time units
	
%--- Interactively fill in missing arguments ---%
if nargin < 4
    fprintf('Not enough input arguments to run automatically.\n');
    inBasename  = input('Base filename for source files, (e.g. "3D_XYZ", no trailing _):','s');
    padlength   = input('Length of frame #s in source files (e.g. 3D_XYZ_xxxx -> 4): ');
    range       = input('Range of frames to export; _START = 0 (e.g. 0:50:1000 to do every 50th frame from start to 1000): ');
    timeNormalization = input('Characteristic time to normalize by (e.g. alfven crossing time or characteristic rotation period. If in doubt hit enter): ');
    if timeNormalization == 0; timeNormalization = 1; end;
end

%--- Initialization ---%
fprintf('Beginning analysis of %i files...\n', numel(range));
analyzedFrameNumber = 1;

if max(round(range) - range) ~= 0; error('ERROR: Frame range is not integer-valued.\n'); end
if min(range) < 0; error('ERROR: Frame range must be nonnegative.\n'); end

range = removeNonexistantEntries(inBasename, padlength, range);
maxFrameno = max(range);

if nargin == 4; timeNormalization = 1; end;

fftVals = [];
timeVals = [];

seedTime = 0;

%--- Loop over given frame range ---%
for ITER = 1:numel(range)
    % Take first guess; Always replace _START
    fname = sprintf('%s_%0*i.mat', inBasename, padlength, range(ITER));
    if range(ITER) == 0; fname = sprintf('%s_START.mat', inBasename); end

    % Check existance; if fails, try _FINAL then give up
    if (exist(fname, 'file') == 0) && (ITER == numel(range))
        fname = sprintf('%s_FINAL.mat', inBasename);
        if exist(fname, 'file') == 0
            % Weird shit is going on. Run away.
            error('UNRECOVERABLE: File existed when checked but is not openable.\n');
        end
    end
    
%    fprintf('Analyzing frame %i... ', analyzedFrameNumber);
    load(fname);

    fprintf('*');

    structName = who('sx_*');
    structName = structName{1};

    eval(sprintf('dataframe = %s;', structName));
    clear -regexp 'sx_';

    fftVals(analyzedFrameNumber,:,:) = computeFrameFFT(dataframe);
    timeVals(analyzedFrameNumber) = sum(dataframe.time.history);
    
    analyzedFrameNumber = analyzedFrameNumber + 1;
%    fprintf('done.\n');

    if ITER ==1;
        velX = dataframe.momX(1,1,1) / dataframe.mass(1,1,1);
	seedTime = 10*numel(find(dataframe.mass(1:floor(end/2),1,1) ~= 1));
    end

    if ITER == numel(range)

	OUTF = fopen('bestmodes.txt','w');

        fprintf('Identifying onset of nonlinearity...\n');
        u = mean(abs(diff(dataframe.time.history(1000:min(2000,end))))); % get jumps in linear regime (i.e. small)
        OoNL = min(find(abs(diff(dataframe.time.history(1000:end))) > 3*u));
%figure(); plot(dataframe.time.history);
        if isempty(OoNL);
            fprintf(OUTF,'Run does not appear to enter nonlinear regime based on dt.\n');
            Tnonlinear = 1e18;
        else;
            OoNL = OoNL + 1000;
            Tnonlinear = sum(dataframe.time.history(1:OoNL));
            fprintf(OUTF, 'Run appears to go nonlinear at step %i, time %g\n', OoNL, Tnonlinear);
        end

	seedTime = sum(dataframe.time.history(1:seedTime));

        fprintf(OUTF, 'Linear values extracted from saveframes %i to %i\n', 1+numel(find(timeVals < seedTime)), numel(find(timeVals < Tnonlinear)));
	fclose(OUTF);

    end
end

fprintf('\nRunning FFT perturbation analyis.\n')

analyzedFrameNumber = analyzedFrameNumber - 2;
fftVals = fftVals(2:end,:,:); timeVals = timeVals(2:end);

bigtime = ones(size(fftVals));
for t = 1:analyzedFrameNumber; bigtime(t,:,:) = timeVals(t); end

N = numel(timeVals(timeVals < Tnonlinear));
M = numel(timeVals(timeVals < seedTime  ))+1;

logfftAmp = log(abs(fftVals(M:N,:,:)));
fftPhase = unwrap(angle(fftVals(M:N,:,:)), 1.5*pi, 1);

[grorate groR] = linearValueFit(bigtime(M:N,:,:), logfftAmp, N-M+1);
[phaserate phaseR] = linearValueFit(bigtime(M:N,:,:), fftPhase, N-M+1);

printoutBestModes(grorate, phaserate, groR, phaseR, velX); 

end

function newrange = removeNonexistantEntries(inBasename, padlength, range)

existrange = [];

for ITER = 1:numel(range)
    % Take first guess; Always replace _START
    fname = sprintf('%s_%0*i.mat', inBasename, padlength, range(ITER));
    if range(ITER) == 0; fname = sprintf('%s_START.mat', inBasename); end

    % Check existance; if fails, try _FINAL then give up
    doesExist = exist(fname, 'file');
    if (doesExist == 0) && (ITER == numel(range))
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

