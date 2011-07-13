function ANALYSIS  = growthAnalysis2(inBasename, padlength, range, timeNormalization)
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

ANALYSIS = [];

%--- Initialization ---%
fprintf('Beginning analysis of %i files...\n', numel(range));
analyzedFrameNumber = 1;

if max(round(range) - range) ~= 0; error('ERROR: Frame range is not integer-valued.\n'); end
if min(range) < 0; error('ERROR: Frame range must be nonnegative.\n'); end

range = removeNonexistantEntries(inBasename, padlength, range);
maxFrameno = max(range);

if nargin == 4; timeNormalization = 1; end;

% Store this to enable the analysis routine to extend a given analysis
ANALYSIS.about.inBasename = inBasename;
ANALYSIS.about.padlength = padlength;
ANALYSIS.about.range = range;
ANALYSIS.about.timeNormalization = timeNormalization;

is2d = [];

%xran = input('X range to analyze: ');
yran = input('# of Y modes to analyze: ');
zran = input('# of Z modes to analyze: ');

ANALYSIS.nFrames = numel(range);
ANALYSIS.nModes = [yran zran];

%--- Loop over given frame range ---%
for ITER = 1:numel(range)
    
    dataframe = frameNumberToData(inBasename, padlength, range(ITER) );
    fprintf('*');
    if mod(ITER, 100) == 0; fprintf('\n'); end

    ANALYSIS.frameTimes(ITER) = sum(dataframe.time.history);
    
    if ITER == 1
        ANALYSIS.equil.rho = dataframe.mass(:,1,1);
        ANALYSIS.equil.ener= dataframe.ener(:,1,1);

        ANALYSIS.equil.mom(1,:) = dataframe.momX(:,1,1)';
        ANALYSIS.equil.mom(2,:) = dataframe.momY(:,1,1)';
        for tp = 1:size(dataframe.mass,1); ANALYSIS.equil.vel(:,tp) = ANALYSIS.equil.mom(:,tp) / ANALYSIS.equil.rho(tp); end
    
        ANALYSIS.equil.B(1,:) = dataframe.magX(:,1,1)';
        ANALYSIS.equil.B(2,:) = dataframe.magY(:,1,1)';

        xd = size(dataframe.mass,1);
        xpre = round(xd/2 - xd/6):round(xd/2 - 6);
        xpost = round(xd/2 + 6):round(xd/2 + xd/6);

        ANALYSIS.gridXvals = cumsum(dataframe.dGrid{1}(:,1,1));

        if size(dataframe.mass,3) == 1; is2d = true; else; is2d = false; end

        ANALYSIS.kyValues   = (0:(ANALYSIS.nModes(1)-1))' * (2*pi/(size(dataframe.mass,2)*dataframe.dGrid{2}));
        ANALYSIS.kyWavenums =  0:(ANALYSIS.nModes(1)-1)';
        ANALYSIS.kzValues   = (0:(ANALYSIS.nModes(2)-1))' * (2*pi/(size(dataframe.mass,3)*dataframe.dGrid{3}));
        ANALYSIS.kzWavenums =  0:(ANALYSIS.nModes(2)-1)';
    end

    xd = size(dataframe.mass,1);
    xpre = round(xd/2 - xd/6):round(xd/2 - 10);
    xpost = round(xd/2 + 10):round(xd/2 + xd/6);

    % This uses a linear extrapolation to track the shock front's position
    % We define that position as being when density is exactly halfway between analytic equilibrium pre & post values
    % This is used to calculate growth rates & omega.
    % It can remain meaningful into the nonlinear regime as long as the shock's position is still functional in Y and Z.
    ANALYSIS.front.X(:,:,ITER) = squeeze(trackFront2(dataframe.mass, ANALYSIS.gridXvals));
    if ITER == 1
        ANALYSIS.gridXvals = ANALYSIS.gridXvals - ANALYSIS.front.X(1,1,1); % place initial shock X at 0.
        ANALYSIS.front.X(1,1,1) = 0; % Reset to be consistent
    end

    ANALYSIS.pre.X = ANALYSIS.gridXvals(xpre);
    ANALYSIS.post.X = ANALYSIS.gridXvals(xpost);

    selectY = 1:ANALYSIS.nModes(1);
    selectZ = 1:ANALYSIS.nModes(2);

    for xi = 1:numel(xpre)
        dq = fft2(squeeze(shiftdim(dataframe.mass(xpre(xi),:,:),1)) - ANALYSIS.equil.rho(xpre(xi)) );
        ANALYSIS.pre.drho(:,:,xi,ITER)= dq(selectY, selectZ);

        dq = fft2(squeeze(shiftdim(dataframe.momX(xpre(xi),:,:)./dataframe.mass(xpre(xi),:,:),1)) - ANALYSIS.equil.mom(1,xpre(xi))/ANALYSIS.equil.rho(xpre(xi)) );
        ANALYSIS.pre.dvx(:,:,xi,ITER) = dq(selectY, selectZ);

        dq = fft2(squeeze(shiftdim(dataframe.momY(xpre(xi),:,:)./dataframe.mass(xpre(xi),:,:),1)) - ANALYSIS.equil.mom(2,xpre(xi))/ANALYSIS.equil.rho(xpre(xi)) );
        ANALYSIS.pre.dvy(:,:,xi,ITER) = dq(selectY, selectZ);

        if size(dataframe.mass,3) > 1 
        dq = fft2(squeeze(shiftdim(dataframe.momZ(xpre(xi),:,:)./dataframe.mass(xpre(xi),:,:),1)) );
        ANALYSIS.pre.dvz(:,:,xi,ITER) = dq(selectY, selectZ);
        end

        dq = fft2(squeeze(shiftdim(dataframe.magX(xpre(xi),:,:),1)) - ANALYSIS.equil.B(1,xpre(xi)) );
        ANALYSIS.pre.dbx(:,:,xi,ITER) = dq(selectY, selectZ);

        dq = fft2(squeeze(shiftdim(dataframe.magY(xpre(xi),:,:),1)) - ANALYSIS.equil.B(2,xpre(xi)) );
        ANALYSIS.pre.dby(:,:,xi,ITER) = dq(selectY, selectZ);

        if size(dataframe.mass,3) > 1
        dq = fft2(squeeze(shiftdim(dataframe.magZ(xpre(xi),:,:),1)));
        ANALYSIS.pre.dbz(:,:,xi,ITER) = dq(selectY, selectZ);
        end
    end
    for xi = 1:numel(xpost)
        dq = fft2(squeeze(shiftdim(dataframe.mass(xpost(xi),:,:),1)) - ANALYSIS.equil.rho(xpost(xi)) );
        ANALYSIS.post.drho(:,:,xi,ITER) = dq(selectY, selectZ);

        dq = fft2(squeeze(shiftdim(dataframe.momX(xpost(xi),:,:)./dataframe.mass(xpost(xi),:,:),1)) - ANALYSIS.equil.mom(1,xpost(xi))/ANALYSIS.equil.rho(xpost(xi)) );
        ANALYSIS.post.dvx(:,:,xi,ITER) = dq(selectY, selectZ);

        dq = fft2(squeeze(shiftdim(dataframe.momY(xpost(xi),:,:)./dataframe.mass(xpost(xi),:,:),1)) - ANALYSIS.equil.mom(2,xpost(xi))/ANALYSIS.equil.rho(xpost(xi)) );
        ANALYSIS.post.dvy(:,:,xi,ITER) = dq(selectY, selectZ);

        if size(dataframe.mass,3) > 1
        dq = fft2(squeeze(shiftdim(dataframe.momZ(xpost(xi),:,:)./dataframe.mass(xpost(xi),:,:),1)) );
        ANALYSIS.post.dvz(:,:,xi,ITER) = dq(selectY, selectZ);
        end

        dq = fft2(squeeze(shiftdim(dataframe.magX(xpost(xi),:,:),1)) - ANALYSIS.equil.B(1,xpost(xi)) );
        ANALYSIS.post.dbx(:,:,xi,ITER) = dq(selectY, selectZ);

        dq = fft2(squeeze(shiftdim(dataframe.magY(xpost(xi),:,:),1)) - ANALYSIS.equil.B(2,xpost(xi)) );
        ANALYSIS.post.dby(:,:,xi,ITER) = dq(selectY, selectZ);

        if size(dataframe.mass,3) > 1
        dq = fft2(squeeze(shiftdim(dataframe.magZ(xpost(xi),:,:),1) ));
        ANALYSIS.post.dbz(:,:,xi,ITER) = dq(selectY, selectZ);
        end
    end

end

for ITER = 1:size(ANALYSIS.front.X,3)
    ANALYSIS.front.FFT(:,:,ITER) = fft2(ANALYSIS.front.X(:,:,ITER));
    ANALYSIS.front.rms(ITER) = sum(sum(sqrt( (ANALYSIS.front.X(:,:,ITER) - mean(mean(ANALYSIS.front.X(:,:,ITER)))).^2  ))) / numel(ANALYSIS.front.X(:,:,ITER));
end


% Use the Grad Student Algorithm to find when the run stops its initial transients and when it goes nonlinear
figno = figure(); plot(log(ANALYSIS.front.rms));
diff(log(ANALYSIS.front.rms(3:end)'),1,1)'
linearFrames = input('Set of frames where line is straight or numbers vaguely constant: ');
close(figno);

fprintf('Run indicated as being in linear regime for saveframes %i to %i inclusive.\n', min(linearFrames), max(linearFrames));

ANALYSIS.lastLinearFrame = frameNumberToData(inBasename, padlength, linearFrames(end) );
ANALYSIS.linearFrames = linearFrames;

fprintf('\nAnalyzing shock front (eta)...\n');

[growthrates growresidual phaserates phaseresidual] = analyzeFront(ANALYSIS.front.FFT, ANALYSIS.frameTimes, linearFrames);
ANALYSIS.linear.omega = phaserates + 1i*growthrates;
ANALYSIS.linear.omegaResidual = phaseresidual + 1i*growresidual;

fprintf('kx/w from post drho: ');
[ANALYSIS.post.drhoKx ANALYSIS.omega_fromdrho2] = analyzePerturbedQ(ANALYSIS.post.drho, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);
fprintf('kx/w from post dv: ');
[ANALYSIS.post.dvxKx ANALYSIS.omega_fromdvx2]   = analyzePerturbedQ(ANALYSIS.post.dvx, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);
[ANALYSIS.post.dvyKx ANALYSIS.omega_fromdvy2]   = analyzePerturbedQ(ANALYSIS.post.dvy, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);
fprintf('kx/w from post db: ');
[ANALYSIS.post.dbxKx ANALYSIS.omega_fromdbx2]   = analyzePerturbedQ(ANALYSIS.post.dbx, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);
[ANALYSIS.post.dbyKx ANALYSIS.omega_fromdby2]   = analyzePerturbedQ(ANALYSIS.post.dby, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);

if is2d == 0
    fprintf('kx/w from dvz/dbz: ');
    [ANALYSIS.post.dvzKx ANALYSIS.omega_fromdvz2] = analyzePerturbedQ(ANALYSIS.post.dvz, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);
    [ANALYSIS.post.dbzKx ANALYSIS.omega_fromdbz2] = analyzePerturbedQ(ANALYSIS.post.dbz, ANALYSIS.post.X, ANALYSIS.frameTimes, 1e-5);
end

fprintf('kx/w from perturbed pre: ');
[ANALYSIS.pre.drhoKx ANALYSIS.omega_fromdrho1] = analyzePerturbedQ(ANALYSIS.pre.drho, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);
[ANALYSIS.pre.dvxKx ANALYSIS.omega_fromdvx1]   = analyzePerturbedQ(ANALYSIS.pre.dvx, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);
[ANALYSIS.pre.dvyKx ANALYSIS.omega_fromdvy1]   = analyzePerturbedQ(ANALYSIS.pre.dvy, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);
[ANALYSIS.pre.dbxKx ANALYSIS.omega_fromdbx1]   = analyzePerturbedQ(ANALYSIS.pre.dbx, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);
[ANALYSIS.pre.dbyKx ANALYSIS.omega_fromdby1]   = analyzePerturbedQ(ANALYSIS.pre.dby, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);

if is2d == 0
    [ANALYSIS.pre.dvzKx ANALYSIS.omega_fromdvz2] = analyzePerturbedQ(ANALYSIS.pre.dvz, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);
    [ANALYSIS.pre.dbzKx ANALYSIS.omega_fromdbz2] = analyzePerturbedQ(ANALYSIS.pre.dbz, ANALYSIS.pre.X, ANALYSIS.frameTimes, 1e-5);
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
    
