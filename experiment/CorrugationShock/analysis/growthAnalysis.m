function [timeVals front kxInfo lastframe] = growthAnalysis(inBasename, padlength, range, timeNormalization)
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

%xran = input('X range to analyze: ');
yran = input('# of Y modes to analyze: ');
zran = input('# of Z modes to analyze: ');

front.X = [];
lastframe = [];
kxInfo = [];

%--- Loop over given frame range ---%
for ITER = 1:numel(range)
    
    dataframe = frameNumberToData(inBasename, padlength, range(ITER) );
    fprintf('*');
    if mod(ITER, 100) == 0; fprintf('\n'); end

    timeVals(ITER) = sum(dataframe.time.history);
    
    % This uses a linear extrapolation to track linear-regime sub-cell variations in the shock front's position
    % We define that position as being when density is exactly halfway between analytic equilibrium pre & post values
    % This is used to calculate growth rates & omega.
    % It can presumably remain meaningful into the nonlinear regime as long as the shock's position is still functional in Y and Z.
    front.X(:,:,ITER) = squeeze(trackFront(dataframe));
end



% Use the Grad Student Algorithm to find when the run stops its initial transients and when it goes nonlinear
figno = figure(); plot(diff(dataframe.time.history(50:end)));
tl = input('Input frame 2x that where dt stops wobbling at the start: ');
th = input('Input frame where dt plunges or last frame: ');
close(figno);

tl = 1+numel(find(timeVals < sum(dataframe.time.history(1:tl))));
th = numel(find(timeVals < sum(dataframe.time.history(1:th))));
linearFrames = tl:th;

fprintf('Run indicated as being in linear regime for saveframes %i to %i inclusive.\n', tl, th);

lastframe = frameNumberToData(inBasename, padlength, range(th) );

fprintf('\nFourier analyzing shock front (eta):\n');

front.FFT = zeros(size(front.X));
for ITER = 1:size(front.X,3)
    front.FFT(:,:,ITER) = fft2(front.X(:,:,ITER));
    front.rms(ITER) = sum(sum(sqrt( (front.X(:,:,ITER) - mean(mean(front.X(:,:,ITER)))).^2  ))) / numel(front.X(:,:,ITER));
end

[front.growthRate front.residualNorm] = analyzeFront(front.FFT, timeVals, linearFrames);

fprintf('\nUsing indicated last linear frame to calculate pre & postshock damping rates\n');



xd = size(lastframe.mass,1);
xpre = round(xd/2 - xd/6):round(xd/2 - 6);
xpost = round(xd/2 + 6):round(xd/2 + xd/6);
middleDx = lastframe.dGrid{1}(round(end/2),1,1);

[kxInfo.dampPre kxInfo.corrPre kxInfo.KY, kxInfo.KZ]   = analyzeDampRates(lastframe, xpre,  yran, zran, middleDx);
[kxInfo.dampPost kxInfo.corrPost] = analyzeDampRates(lastframe, xpost, yran, zran, middleDx);

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

function dataframe = frameNumberToData(inBasename, padlength, frameno)
    % Take first guess; Always replace _START
    fname = sprintf('%s_%0*i.mat', inBasename, padlength, frameno);
    if frameno == 0; fname = sprintf('%s_START.mat', inBasename); end

    % Check existance; if fails, try _FINAL then give up
    if exist(fname, 'file') == 0
        fname = sprintf('%s_FINAL.mat', inBasename);
        if exist(fname, 'file') == 0
            % Weird shit is going on. Run away.
            error('UNRECOVERABLE: File existed when checked but is not openable.\n');
        end
    end

    % Load the next frame into workspace; Assign it to a standard variable name.    
    load(fname);
    structName = who('sx_*');
    structName = structName{1};

    eval(sprintf('dataframe = %s;', structName));
end
    
