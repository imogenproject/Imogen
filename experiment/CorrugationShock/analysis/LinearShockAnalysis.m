classdef LinearShockAnalysis < handle

    properties (Constant = true)
        version = 1.0;
    end

    properties(SetAccess = public, GetAccess = public)

    end

    properties(SetAccess = protected, GetAccess = public)
        nModes; % number of transverse modes in [Y Z] direction
        nFrames;
        is2d;   % True if run was of z extent 1, otherwise false.

        originalSrcDirectory; % Set at analyze time, the directory the analyzed data resided in.

        frameTimes;
        gridXvals;
        kyValues; kyWavenums; % Transverse mode values and integer wavenumbers
        kzValues; kzWavenums;

        linearFrames;
        lastLinearFrame;

        equil; % Structure containing equilibrium data
        front; % Structure containing information about the shock front
        pre;   % Structures containing pre and post shock spectral data
        post;
        omega;
    end

    properties(SetAccess = protected, GetAccess = protected)
        inputBasename;
        inputPadlength;
        inputFrameRange;
        maxFrameno;

        manfit_state;
    end

methods % SET

end

methods (Access = public)

    function obj = LinearShockAnalysis()

        obj.selectFileset();

        an = input('Run fourier analysis now? ','s');
        if strcmp(an,'y') || strcmp(an,'yes'); obj.performFourierAnalysis(); end

        an = input('Perform automated mode analysis? ','s');
        if strcmp(an,'y') || strcmp(an,'yes'); obj.curveFit_automatic(); end

    end

    function help(obj)

    fprintf('This is the Imogen corrugation shock analyzer.\n\nWhen you created me, I took you through gathering a set of files; Access that again with selectFileset().\nTo gather the fourier-space data, run performFourierAnalysis().\nTo have me try to automagically fit all the wavevectors and growth rates, run curveFit_automatic().\nI''m bad at determining what data''s actually valid; run curveFit_manual() on every mode you intend to actually believe, if only to check.\n\nYou can save me using "save(''filename.mat'',''myname'')" at any time and load me back the same way.\n');

    end

    function selectFileset(obj)
        obj.inputBasename  = input('Base filename for source files, (e.g. "3D_XYZ", no trailing _):','s');
        obj.inputPadlength   = input('Length of frame #s in source files (e.g. 3D_XYZ_xxxx -> 4): ');
        obj.inputFrameRange       = input('Range of frames to export; _START = 0 (e.g. 0:50:1000 to do every 50th frame from start to 1000): ');
%    timeNormalization = input('Characteristic time to normalize by (e.g. alfven crossing time or characteristic rotation period. If in doubt hit enter): ');
%    if timeNormalization == 0; timeNormalization = 1; end;

        if max(round(obj.inputFrameRange) - obj.inputFrameRange) ~= 0; error('ERROR: Frame obj.inputFrameRange is not integer-valued.\n'); end
        if min(obj.inputFrameRange) < 0; error('ERROR: Frame obj.inputFrameRange must be nonnegative.\n'); end

        obj.inputFrameRange = obj.removeNonexistantEntries(obj.inputBasename, obj.inputPadlength, obj.inputFrameRange);
        obj.maxFrameno = max(obj.inputFrameRange);

        obj.originalSrcDirectory = pwd();
    end

    function performFourierAnalysis(obj)

        obj.nModes(1) = input('# of Y modes to analyze: ');
        obj.nModes(2) = input('# of Z modes to analyze: ');

        yran = obj.nModes(1); zran = obj.nModes(2); % For shortform convenience

        obj.nFrames = numel(obj.inputFrameRange);


        %--- Loop over given frame obj.inputFrameRange ---%
        for ITER = 1:numel(obj.inputFrameRange)

            dataframe = obj.frameNumberToData(obj.inputBasename, obj.inputPadlength, obj.inputFrameRange(ITER) );
            fprintf('*');
            if mod(ITER, 50) == 0; fprintf('\n'); end

            obj.frameTimes(ITER) = sum(dataframe.time.history);

            if ITER == 1 % Extract equilibrium data from the first frame
                obj.equil.rho = dataframe.mass(:,1,1);
                obj.equil.ener= dataframe.ener(:,1,1);

                obj.equil.mom(1,:) = dataframe.momX(:,1,1)';
                obj.equil.mom(2,:) = dataframe.momY(:,1,1)';
                for tp = 1:size(dataframe.mass,1); obj.equil.vel(:,tp) = obj.equil.mom(:,tp) / obj.equil.rho(tp); end

                obj.equil.B(1,:) = dataframe.magX(:,1,1)';
                obj.equil.B(2,:) = dataframe.magY(:,1,1)';

                xd = size(dataframe.mass,1);
                xpre = round(xd/2 - xd/6):round(xd/2 - 6);
                xpost = round(xd/2 + 6):round(xd/2 + xd/6);

                obj.gridXvals = cumsum(dataframe.dGrid{1}(:,1,1));

                if size(dataframe.mass,3) == 1; obj.is2d = true; else; obj.is2d = false; end

                obj.kyValues   = (0:(obj.nModes(1)-1))' * (2*pi/(size(dataframe.mass,2)*dataframe.dGrid{2}));
                obj.kyWavenums =  0:(obj.nModes(1)-1)';
                obj.kzValues   = (0:(obj.nModes(2)-1))' * (2*pi/(size(dataframe.mass,3)*dataframe.dGrid{3}));
                obj.kzWavenums =  0:(obj.nModes(2)-1)';
            end

            xd = size(dataframe.mass,1);
            xpre = round(xd/2 - xd/6):round(xd/2 - 4);
            xpost = round(xd/2 + 4):round(xd/2 + xd/6);

            % This uses a linear extrapolation to track the shock front's position
            % We define that position as being when density is exactly halfway between analytic equilibrium pre & post values
            % This is used to calculate growth rates & omega.
            % It can remain meaningful into the nonlinear regime as long as the shock's position is still functional in Y and Z.
            obj.front.X(:,:,ITER) = squeeze(trackFront2(dataframe.mass, obj.gridXvals));
            if ITER == 1
                obj.gridXvals = obj.gridXvals - obj.front.X(1,1,1); % place initial shock X at 0.
                obj.front.X(1,1,1) = 0; % Reset to be consistent
            end

            obj.pre.X = obj.gridXvals(xpre);
            obj.post.X = obj.gridXvals(xpost);

            selectY = 1:obj.nModes(1);
            selectZ = 1:obj.nModes(2);

            for xi = 1:numel(xpre)
                dq = fft2(squeeze(shiftdim(dataframe.mass(xpre(xi),:,:),1)) - obj.equil.rho(xpre(xi)) );
                obj.pre.drho(:,:,xi,ITER)= dq(selectY, selectZ);

                dq = fft2(squeeze(shiftdim(dataframe.momX(xpre(xi),:,:)./dataframe.mass(xpre(xi),:,:),1)) - obj.equil.mom(1,xpre(xi))/obj.equil.rho(xpre(xi)) );
                obj.pre.dvx(:,:,xi,ITER) = dq(selectY, selectZ);

                dq = fft2(squeeze(shiftdim(dataframe.momY(xpre(xi),:,:)./dataframe.mass(xpre(xi),:,:),1)) - obj.equil.mom(2,xpre(xi))/obj.equil.rho(xpre(xi)) );
                obj.pre.dvy(:,:,xi,ITER) = dq(selectY, selectZ);

                if size(dataframe.mass,3) > 1
                dq = fft2(squeeze(shiftdim(dataframe.momZ(xpre(xi),:,:)./dataframe.mass(xpre(xi),:,:),1)) );
                obj.pre.dvz(:,:,xi,ITER) = dq(selectY, selectZ);
                end

                dq = fft2(squeeze(shiftdim(dataframe.magX(xpre(xi),:,:),1)) - obj.equil.B(1,xpre(xi)) );
                obj.pre.dbx(:,:,xi,ITER) = dq(selectY, selectZ);

                dq = fft2(squeeze(shiftdim(dataframe.magY(xpre(xi),:,:),1)) - obj.equil.B(2,xpre(xi)) );
                obj.pre.dby(:,:,xi,ITER) = dq(selectY, selectZ);

                if size(dataframe.mass,3) > 1
                dq = fft2(squeeze(shiftdim(dataframe.magZ(xpre(xi),:,:),1)));
                obj.pre.dbz(:,:,xi,ITER) = dq(selectY, selectZ);
                end
            end
            for xi = 1:numel(xpost)
                dq = fft2(squeeze(shiftdim(dataframe.mass(xpost(xi),:,:),1)) - obj.equil.rho(xpost(xi)) );
                obj.post.drho(:,:,xi,ITER) = dq(selectY, selectZ);

                dq = fft2(squeeze(shiftdim(dataframe.momX(xpost(xi),:,:)./dataframe.mass(xpost(xi),:,:),1)) - obj.equil.mom(1,xpost(xi))/obj.equil.rho(xpost(xi)) );
                obj.post.dvx(:,:,xi,ITER) = dq(selectY, selectZ);

                dq = fft2(squeeze(shiftdim(dataframe.momY(xpost(xi),:,:)./dataframe.mass(xpost(xi),:,:),1)) - obj.equil.mom(2,xpost(xi))/obj.equil.rho(xpost(xi)) );
                obj.post.dvy(:,:,xi,ITER) = dq(selectY, selectZ);

                if size(dataframe.mass,3) > 1
                dq = fft2(squeeze(shiftdim(dataframe.momZ(xpost(xi),:,:)./dataframe.mass(xpost(xi),:,:),1)) );
                obj.post.dvz(:,:,xi,ITER) = dq(selectY, selectZ);
                end

                dq = fft2(squeeze(shiftdim(dataframe.magX(xpost(xi),:,:),1)) - obj.equil.B(1,xpost(xi)) );
                obj.post.dbx(:,:,xi,ITER) = dq(selectY, selectZ);

                dq = fft2(squeeze(shiftdim(dataframe.magY(xpost(xi),:,:),1)) - obj.equil.B(2,xpost(xi)) );
                obj.post.dby(:,:,xi,ITER) = dq(selectY, selectZ);

                if size(dataframe.mass,3) > 1
                dq = fft2(squeeze(shiftdim(dataframe.magZ(xpost(xi),:,:),1) ));
                obj.post.dbz(:,:,xi,ITER) = dq(selectY, selectZ);
                end
            end
        end

        for ITER = 1:size(obj.front.X,3)
            obj.front.FFT(:,:,ITER) = fft2(obj.front.X(:,:,ITER));
            obj.front.rms(ITER) = sum(sum(sqrt( (obj.front.X(:,:,ITER) - mean(mean(obj.front.X(:,:,ITER)))).^2  ))) / numel(obj.front.X(:,:,ITER));
        end
    

    end

    function curveFit_automatic(obj)
        % Use the Grad Student Algorithm to find when the run stops its initial transients and when it goes nonlinear
        figno = figure(); plot(log(obj.front.rms));
        diff(log(obj.front.rms(3:end)'),1,1)'
        linearFrames = input('Set of frames where line is straight or numbers vaguely constant: ');
        close(figno);

        fprintf('Run indicated as being in linear regime for saveframes %i to %i inclusive.\n', min(linearFrames), max(linearFrames));
        obj.lastLinearFrame = obj.frameNumberToData(obj.inputBasename, obj.inputPadlength, linearFrames(end) );
        obj.linearFrames = linearFrames;

        fprintf('\nAnalyzing shock front (eta)...\n');

        [growthrates growresidual phaserates phaseresidual] = analyzeFront(obj.front.FFT, obj.frameTimes, linearFrames);
        obj.omega.front = phaserates + 1i*growthrates;
        obj.omega.frontResidual = phaseresidual + 1i*growresidual;

        fprintf('kx/w from post drho: ');
        [obj.post.drhoKx obj.omega.fromdrho2] = analyzePerturbedQ(obj.post.drho, obj.post.X, obj.frameTimes, 1e-5);
        fprintf('kx/w from post dv: ');
        [obj.post.dvxKx obj.omega.fromdvx2]   = analyzePerturbedQ(obj.post.dvx, obj.post.X, obj.frameTimes, 1e-5);
        [obj.post.dvyKx obj.omega.fromdvy2]   = analyzePerturbedQ(obj.post.dvy, obj.post.X, obj.frameTimes, 1e-5);
        fprintf('kx/w from post db: ');
        [obj.post.dbxKx obj.omega.fromdbx2]   = analyzePerturbedQ(obj.post.dbx, obj.post.X, obj.frameTimes, 1e-5);
        [obj.post.dbyKx obj.omega.fromdby2]   = analyzePerturbedQ(obj.post.dby, obj.post.X, obj.frameTimes, 1e-5);

        if obj.is2d == 0
            fprintf('kx/w from dvz/dbz: ');
            [obj.post.dvzKx obj.omega.fromdvz2] = analyzePerturbedQ(obj.post.dvz, obj.post.X, obj.frameTimes, 1e-5);
            [obj.post.dbzKx obj.omega.fromdbz2] = analyzePerturbedQ(obj.post.dbz, obj.post.X, obj.frameTimes, 1e-5);
        end

        fprintf('kx/w from perturbed pre: ');
        [obj.pre.drhoKx obj.omega.fromdrho1] = analyzePerturbedQ(obj.pre.drho, obj.pre.X, obj.frameTimes, 1e-5);
        [obj.pre.dvxKx obj.omega.fromdvx1]   = analyzePerturbedQ(obj.pre.dvx, obj.pre.X, obj.frameTimes, 1e-5);
        [obj.pre.dvyKx obj.omega.fromdvy1]   = analyzePerturbedQ(obj.pre.dvy, obj.pre.X, obj.frameTimes, 1e-5);
        [obj.pre.dbxKx obj.omega.fromdbx1]   = analyzePerturbedQ(obj.pre.dbx, obj.pre.X, obj.frameTimes, 1e-5);
        [obj.pre.dbyKx obj.omega.fromdby1]   = analyzePerturbedQ(obj.pre.dby, obj.pre.X, obj.frameTimes, 1e-5);

        if obj.is2d == 0
            [obj.pre.dvzKx obj.omega.fromdvz2] = analyzePerturbedQ(obj.pre.dvz, obj.pre.X, obj.frameTimes, 1e-5);
            [obj.pre.dbzKx obj.omega.fromdbz2] = analyzePerturbedQ(obj.pre.dbz, obj.pre.X, obj.frameTimes, 1e-5);
        end

    end

    function curveFit_manual(obj)

        % Determine which mode and structure we're trying to fit

        ymode = input('Y mode # to fit: ');
        zmode = input('Z mode # to fit: ');
        % Throw up a contoured surf with registered callback

% memory:
% .kx = the stored kx value guess
% .w  = the stored omega value guess
% df = the change if the arrows are pressed
% whofit = which curve we're trying to fit (toggle with space
        obj.manfit_state.y = ymode;
        obj.manfit_state.z = zmode;

        obj.manfit_state.kx = obj.post.drhoKx(ymode, zmode);
        obj.manfit_state.w  = obj.omega.fromdrho2(ymode, zmode);
        obj.manfit_state.df = .001;
        obj.manfit_state.whofit = 1;
        obj.manfit_state.a0 = 0;

        qty = input('Quantity: (1) drhopost (2) dvxpost (3) dvypost (4) dbxpost (5) dbypost (6) drhopre (7) dvxpre (8) dvypre (9) dbxpre (10) dbypre: ');
        obj.manfit_state.qty = qty;
        dq = [];
        switch(qty)
            case 1;  dq = squeeze(obj.post.drho(ymode, zmode,:,:));
            case 2;  dq = squeeze( obj.post.dvx(ymode, zmode,:,:));
            case 3;  dq = squeeze( obj.post.dvy(ymode, zmode,:,:));
            case 4;  dq = squeeze( obj.post.dbx(ymode, zmode,:,:));
            case 5;  dq = squeeze( obj.post.dby(ymode, zmode,:,:));
            case 6;  dq = squeeze( obj.pre.drho(ymode, zmode,:,:));
            case 7;  dq = squeeze(  obj.pre.dvx(ymode, zmode,:,:));
            case 8;  dq = squeeze(  obj.pre.dvy(ymode, zmode,:,:));
            case 9;  dq = squeeze(  obj.pre.dbx(ymode, zmode,:,:));
            case 10; dq = squeeze(  obj.pre.dby(ymode, zmode,:,:));
        end
        
        figure('KeyPressFcn',{@manualfitter_callback, @obj.manfit_setKW, @obj.manfit_memory, obj, obj.post.X, obj.frameTimes, dq});
        %manualfitter_callback(src, eventdata, setterfcn, memfcn, memory, xaxis, yaxis, datain
    end

    function manfit_memory(obj, s)
        % Store manual fit's data in our memory
        obj.manfit_state = s;
    end

    function manfit_setKW(obj, y, z, w, kx, qty)
        switch(qty);
            case 1 ; obj.post.drhoKx(y,z) = kx; obj.omega.fromdrho2(y,z) = w;
            case 2 ; obj.post.dvxKx(y,z) = kx; obj.omega.fromdvx2(y,z) = w;
            case 3 ; obj.post.dvyKx(y,z) = kx; obj.omega.fromdvy2(y,z) = w;
            case 4 ; obj.post.dbxKx(y,z) = kx; obj.omega.fromdbx2(y,z) = w;
            case 5 ; obj.post.dbyKx(y,z) = kx; obj.omega.fromdby2(y,z) = w;
            case 6 ; obj.pre.drhoKx(y,z) = kx; obj.omega.fromdrho1(y,z) = w;
            case 7 ; obj.pre.dvxKx(y,z) = kx; obj.omega.fromdvx2(y,z) = w;
            case 8 ; obj.pre.dvyKx(y,z) = kx; obj.omega.fromdvy2(y,z) = w;
            case 9 ; obj.pre.dbxKx(y,z) = kx; obj.omega.fromdbx2(y,z) = w;
            case 10; obj.pre.dbyKx(y,z) = kx; obj.omega.fromdby2(y,z) = w;
        end
    end
end % Public methods

methods % SET

end

methods (Access = protected)

    function newframeranges = removeNonexistantEntries(obj, namebase, padsize, frameranges)

    existframeranges = [];

        for ITER = 1:numel(frameranges)
            % Take first guess; Always replace _START
            fname = sprintf('%s_%0*i.mat', namebase, padsize, frameranges(ITER));
            if frameranges(ITER) == 0; fname = sprintf('%s_START.mat', namebase); end

            % Check existance; if fails, try _FINAL then give up
            doesExist = exist(fname, 'file');
            if (doesExist == 0) && (ITER == numel(frameranges))
                fname = sprintf('%s_FINAL.mat', namebase);
                doesExist = exist(fname, 'file');
            end
        
            if doesExist ~= 0; existframeranges(end+1) = ITER; end
        end

        newframeranges = frameranges(existframeranges);
        if numel(newframeranges) ~= numel(frameranges);
            fprintf('WARNING: Removed %i entries that could not be opened from list.\n', numel(frameranges)-numel(newframeranges));
        end

        if numel(newframeranges) == 0;
            error('UNRECOVERABLE: No files indicated existed. Perhaps remove trailing _ from base name?\n'); 
        end

    end

    function dataframe = frameNumberToData(obj, namebase, padsize, frameno)
        % Take first guess; Always replace _START
        fname = sprintf('%s_%0*i.mat', namebase, padsize, frameno);
        if frameno == 0; fname = sprintf('%s_START.mat', namebase); end

        % Check existance; if fails, try _FINAL then give up
        if exist(fname, 'file') == 0
            fname = sprintf('%s_FINAL.mat', namebase);
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

end % protected methods;

end % class    
