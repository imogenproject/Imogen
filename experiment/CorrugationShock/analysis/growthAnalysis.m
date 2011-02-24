function [timeVals front selectRegion] = growthAnalysis(inBasename, padlength, range, timeNormalization)
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

xran = input('X range to analyze: ');
yran = input('Y range to analyze: ');
zran = input('Z range to analyse: ');

front.X = [];
lastframe = [];
selectRegion = [];

rhoside(1) = 1;
rhoside(2) = 1;
whichside = 0;

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


    % Load the next frame into workspace; Assign it to a standard variable name.    
    load(fname); fprintf('*');
    structName = who('sx_*');
    structName = structName{1};

    eval(sprintf('dataframe = %s;', structName));
    clear -regexp 'sx_';

    % Grab the mass density postshock.
    if ITER == 1;
        rhoside(2) = dataframe.mass(end,1,1);
        if max(xran) < size(dataframe.mass,1)/2; whichside = 1; end
        if min(xran) > size(dataframe.mass,1)/2; whichside = 2; end
    end

    %  Acquire mode and time data for the block requested to be examined
    if whichside > 0
        [selectRegion.damp(ITER,:,:) selectRegion.corr(ITER,:,:) selectRegion.mdr(ITER)] = analyzeDampRates(dataframe, xran, yran, zran, rhoside(whichside), dataframe.dGrid{1}(round(end/2),1,1) );
    end

%    fftVals(ITER,:,:) = computeFrameFFT(dataframe, xran, yran, zran);
    timeVals(ITER) = sum(dataframe.time.history);
    
    % This uses a linear extrapolation to track linear-regime sub-cell variations in the shock front's position
    % We define that position as being when density is exactly halfway between analytic equilibrium pre & post values
    front.X(:,:,ITER) = squeeze(trackFront(dataframe));

    % Hold onto final data frame - it's got all the timesteps stored for us to look at
    if ITER == numel(range); lastframe = dataframe; end
end

% Use the Grad Student Algorithm to find when the run stops its initial transients and when it goes nonlinear
figno = figure(); plot(diff(dataframe.time.history(50:end)));
tl = input('Input frame 2x that where dt stops wobbling at the start: ');
th = input('Input frame where dt plunges or last frame: ');
close(figno);

tl = 1+numel(find(timeVals < sum(dataframe.time.history(1:tl))));
th = numel(find(timeVals < sum(dataframe.time.history(1:th))));
linearFrames = tl:th;

fprintf('Run indicated as being in linear regime for frames %i to %i inclusive.\n', tl, th);

%fprintf('\nRunning FFT perturbation analyis of selected region.\n')
%analyzedFrameNumber = ITER - 2;
%fftVals = fftVals(2:end,:,:); timeVals = timeVals(2:end);

%bigtime = ones(size(fftVals));
%for t = 1:analyzedFrameNumber; bigtime(t,:,:) = timeVals(t); end

%logfftAmp = log(abs(fftVals(linearFrames,:,:)));
%fftPhase = unwrap(angle(fftVals(linearFrames,:,:)), 1.5*pi, 1);

%[grorate groR] = linearValueFit(bigtime(linearFrames,:,:), logfftAmp, numel(linearFrames));
%[phaserate phaseR] = linearValueFit(bigtime(linearFrames,:,:), fftPhase, numel(linearFrames));

%printoutBestModes(grorate, phaserate, groR, phaseR, velX); 

fprintf('\nDoing FFT analysis and linear fit of shock front\n');

front.FFT = zeros(size(front.X));
for ITER = 1:size(front.X,3)
    front.FFT(:,:,ITER) = fft2(front.X(:,:,ITER));
    front.rms(ITER) = sum(sum(sqrt( (front.X(:,:,ITER) - mean(mean(front.X(:,:,ITER)))).^2  ))) / numel(front.X(:,:,ITER));
end

[front.growthRate front.residualNorm] = analyzeFront(front.FFT, timeVals, linearFrames);


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

