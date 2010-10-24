classdef DiskAnalyzer < Analyzer
% Analyzer subclass for Kojima disk analysis.
	
%===================================================================================================
	properties (Constant = true, Transient = true) %							C O N S T A N T	 [P]
    end%CONSTANT
	
%===================================================================================================
    properties (SetAccess = public, GetAccess = public) %							P U B L I C  [P]
    end %PUBLIC

%===================================================================================================
    properties (SetAccess = protected, GetAccess = protected) %				   P R O T E C T E D [P]
    end %PROTECTED	
	
%===================================================================================================
    methods %																	  G E T / S E T  [M]
	end%GET/SET
	
%===================================================================================================
    methods (Access = public) %														P U B L I C  [M]

%___________________________________________________________________________________________________
        function obj = DiskAnalyzer()
            obj.path = 'workspace'; % use 'workspace' or '/path/to/matfiles'
            obj.plotType = 'none';
        end

%___________________________________________________________________________________________________ modeAnalyze
        function result = modeAnalyze(obj,radialIndex, optString)
            % This function sequentially loads disk MAT files from a list, unwraps them, and performs a mode
            % spectrum analysis.
            %>> radialIndex: The radial cell # to examine modes at; set to zero for radial averaging
            %<< modeinfo   : The 12-by-N returned information on the first 12 modes

            if ischar(radialIndex)
                if strcmpi(radialIndex,'help')
                    fprintf('Mode analysis help: Argument 1 can be either\n	[x] to examine radial index x,\n [x y] to examine a range\n	a cell array of [x] and [x y].\n\nIf the second argument is "eigenfunction", modeanalyze will examine regularly spaced bins and argument 1 should be either\n	[x y] to look at (y-x) bins of size 1\n	[x s y] to look at bins of size s stepping from x to y.\n');
                end
                return;
            end

            switch nargin;
                case 0; fprintf('Averaging over all R; Normalizing by m=0 @ t=0\n'); radialindex = 0; %#ok<NASGU>
                case 1;
                case 2;
                    if strcmpi(optString,'eigenfunction')

                        switch numel(radialIndex)
                            case 1; fprintf('Error: require arg #1 = either [in step out] or [in out] to create eigenfunction\n'); result = 0; return;
                            case 2; radBins = radialIndex(1):radialIndex(2); dstep = 1;
                            case 3; radBins = radialIndex(1):radialIndex(2):radialIndex(3); dstep = radialIndex(2);
                        end

                        radialIndex = cell(numel(radBins),1);
                        result.zaxis = radBins + .5*dstep;

                        for y = 1:numel(radBins);
                            radialIndex{y} = [radBins(y) radBins(y)+dstep-1];
                        end
                    end
            end

            X = length(obj.matNames); %<----- number of files to be analyzed

            % Get storage for all the modes we will be analyzing
            result.yaxis = 1:12;
            result.xaxis = zeros(1,X);

            if ~iscell(radialIndex)
                q = radialIndex;
                radialIndex = cell(1);
                radialIndex{1} = q;
            end

            result.modeinfo = zeros(12,X,numel(radialIndex));
            result.growthrate = zeros(12,2,numel(radialIndex));

            for x = 1:X 

                obj.index = x;   %<----- set data index
                sxcur = obj.data; %<----- retrieve data

                result.xaxis(x) = sum(sxcur.time.history) / (2*pi); % Normalize to chirps

                % Perform polar transform, rescale by radius while radially summing
                mpolar = diskUnwrap(sxcur.mass);

                % Calculate radius vector; Rescale z axis appropriately if it exists
                if x == 1;
                    origrad = size(sxcur.mass, 1) * sxcur.dGrid{1} * .5;
                    newdr = origrad / size(mpolar, 1);
                    newr = cumsum(ones(1, size(mpolar,1)) * newdr);
                    if isfield(result, 'zaxis')
                        result.zaxis = result.zaxis * sxcur.dGrid{1} * .5;
                    end
                end

                % Perform the Fourier analysis on all requested bins
                for y = 1:numel(radialIndex)
                    [result.modeinfo(:,x,y) nc] = obj.fftSlice(mpolar, newr, radialIndex{y});
                    if x == 1; result.nConst(y) = nc; end
                end

            end

            for y = 1:numel(radialIndex);
                result.modeinfo(:,:,y) = result.modeinfo(:,:,y) ./ result.nConst(y);

                if numel(find(result.xaxis >= 1)) > 2
                    for m = 1:12;
                        result.growthrate(m,:,y) = obj.findGrowthRate(result.xaxis, result.modeinfo(m,:,y), 1);
                    end
                else
                    fprintf('Insufficient time elapsed to guess at growth rate\n');
                end
            end

        end

    end%PUBLIC
	
%===================================================================================================	
	methods (Access = protected) %											P R O T E C T E D    [M]
	end%PROTECTED
		
%===================================================================================================	
	methods (Static = true) %													  S T A T I C    [M]

%___________________________________________________________________________________________________ fftSlice    
% Takes unwrapped disk, radius vector and slice selection; Returns FFT
        function [result ncst] = fftSlice(fpolar, r, slice)        
            switch numel(slice)
                case 1;
                    if slice ~= 0; r(1:(slice-1)) = 0; r((slice+1):end) = 0; end
                case 2; r(1:(slice(1)-1)) = 0; r((slice(2)+1):end) = 0;
                otherwise;
            end     
            result = fft(r*fpolar);       
            ncst = abs(result(1));
            result = abs(result(2:13)); % return only m=0 and first 12 harmonics        
        end

%___________________________________________________________________________________________________ findGrowthRate        
% Given time and amplitude, calculate growth rate characteristic time
        function result = findGrowthRate(t, a, tmin)        
            x = numel(find(t < tmin))+1;        
            t = t(x:end);
            a = a(x:end);        
            fcn = polyfit(t, log(a), 1);       
            result = fcn;
        end       
        
%___________________________________________________________________________________________________ analyze_KojimaDiskEquilibrium     
        function [ imbalanceMetric toverw drhodt] = analyze_KojimaDiskEquilibrium(mats)
            % analyze_KojimaDiskEquilibrium: A function designed to highlight the stability or instability
            % of actually any flow situation - pass it an Imogen 2D save file and it will
            % return the magnitude of the time derivative of the momentum conservation equation.

            rho = mats.mass;
            momX = mats.momX;
            momY = mats.momY;
            phi = mats.grav;

            griddim = size(rho); griddim = griddim(1);

            gamma = mats.gamma;
            delta = mats.dGrid;
            %deevee = delta{1}*delta{2}*delta{3};

            %%% Compute magnitude of net force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            pressure = rho .^ (gamma); % polytrope EOS

            % fivePtDerivative doesn't seem to work on plain arrays - ???
            % Get pressure and potential gradients
            %pressx = fivePtDerivative(pressure, 1, delta{1});
            %pressy = fivePtDerivative(pressure, 2, delta{2});
            [pressy pressx] = gradient(pressure, delta{1}, delta{2});

            %gravx = fivePtDerivative(pressure, 1, delta{1});
            %gravy = fivePtDerivative(pressure, 2, delta{2});
            [gravy gravx] = gradient(phi, delta{1}, delta{2});

            % Get radial distances and angles
            oneddist = (1:griddim) - floor(griddim/2) - .5; % one-by-n
            oneddist = ( ones(griddim, 1) * oneddist ) .* delta{1}; % n-by-n, distance in Y dir.
            arrdist = sqrt( oneddist.^2 + (oneddist').^2);

            % Compute angles WRT center. cosine = x/r, sine = y/r
            costheta = oneddist' ./ arrdist;
            sintheta= oneddist  ./ arrdist;

            magp = momX.^2 + momY.^2;

            % Get force due to acceleration
            fcent = magp ./ ( rho .* arrdist  );
            fcentX = fcent .* costheta;
            fcentY = fcent .* sintheta;

            % Calculate Fnet = mv^2/r - rho grad phi - grad P
            ForceX = fcentX - rho .* gravx - pressx;
            ForceY = fcentY - rho .* gravy - pressy;

            imbalanceMetric = sqrt(ForceX.^2 + ForceY.^2);

            %%% Compute T/W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kinetic energy = integral(p^2 / 2rho)
            T = .5 * sum(sum(sum(magp ./ rho)));
            % Potential energy = integral(-rho / r)
            W = -sum(sum(sum(rho ./ arrdist)));

            toverw = abs(T / W);

            %%% Compute drho/dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [drhody drhodx] = gradient(rho, delta{1}, delta{2});
            [dvdy dvdx] = gradient(momY ./ rho, delta{1}, delta{2});

            % Compute rho dot div(v) + v dot grad(rho)
            drhodt = rho .* (dvdx + dvdy) + (momY ./ rho) .* drhody + (momX ./ rho) .* drhodx;

        end

%___________________________________________________________________________________________________ diskUnwrap             
        function polmass = diskUnwrap(mass)
            % Function applies a polar to cartesian transform on a disk simulation,
            % Effectively turning x and y into r and theta
            % >> mass : The input square array (not necessarily mass) to be converted
            % << polmass: The transformed result

            grid = size(mass);

            % R and Phi components for transform
            rho = 0:.5:floor(grid(1)/2)-.5;
            phi = (0:.5*pi/(floor(grid(2)/2)):2*pi) - .5*pi/floor(grid(2)/2);

            % Map to cartesian components
            mu = exp(1i*phi);
            polarrayx = rho' * real(mu) + grid(1)/2 + .5;
            polarrayy = rho' * imag(mu) + grid(2)/2 + .5;

            %polarrayx(polarrayx < 1)'
            %polarrayx(polarrayx > 2048)'

            % Generate interpolation
            polmass = interp2(mass,polarrayx,polarrayy);

        end
        
%___________________________________________________________________________________________________ diskRewrap        
        function mass = diskRewrap(polmass)
            % Function applies a polar to cartesian transform on a disk simulation,
            % Effectively turning x and y into r and theta
            % >> mass : The input square array (not necessarily mass) to be converted
            % << polmass: The transformed result


            grid = size(polmass,1)*[1 1];

            % R and Phi components for transform
            [X Y] = ndgrid(1:grid(1),1:grid(2));

            X = X - grid(1)/2 - .5;
            Y = Y - grid(1)/2 - .5;

            rho = sqrt(X.^2 + Y.^2)*2 + .5;
            phi = (-angle(rot90(X + 1i*Y,-1)) + pi)*2*grid(1)/(2*pi) + 1.1 ;

            %figure(3); imagesc(phi)
            %figure(4); imagesc(rho)
            %figure(5); imagesc(polmass)
            % Generate interpolation
            mass = interp2(polmass, phi, rho);

            mass(isnan(mass)) = 1e-10;

        end

%___________________________________________________________________________________________________ plotTimeSpectrum
        function ts = plotTimeSpectrum(stepby, theend, numzeros, DIR) %#ok<INUSD>
            % Started in a run directory,
            % 1. Loads next 2d slice
            % 2. Transforms mass to polar coordinates and sums
            % 3. Generates FFT
            % 4. Clear slice's sx_... and goes to 1

            stepnum = 0;

            ts = zeros([round(theend/stepby + 1) 16]);

            while stepnum <= theend
                % Get frame
                if stepnum == 0; frmstring='START'; else
                    %if stepnum == theend; frmstring='END'; else
                        frmstring = sprintf('%05i',stepnum);
                    %end
                end

                eval(['load 2D_XY_' frmstring '.mat']);
                eval(['dataunit = sx_XY_' frmstring ';']);

                mpol = cartesianToPolar(dataunit.mass);
                s = size(mpol);
                r = ones([1 s(1)]) * dataunit.dGrid{1};
                r = cumsum(r);

                mpol = r * mpol; % Radially integrate each slice to find mass(theta)
                mpol = mpol - mean(mpol);

                fourier = real(fft(mpol));
                % Compute azimuthal mass perturbations
                ts((stepnum / stepby) + 1,:) = fourier(1:16);
                stepnum = stepnum + stepby;

                fprintf('*');
            end
            fprintf('\n');

        end

%___________________________________________________________________________________________________ trackArmPhase        
        function [R Theta] = trackArmPhase(pdata, r0)
            % >> pdata: unwrapped disk information
            % << R: Radius of arm
            % << Theta: phase

            % Find how many radial elements there are to examine
            nRelements = size(pdata, 1) - r0(1) + 1;

            % Create the result holders
            R = [];
            Theta = [];

            % Set the initial 'phase' and the width to search for maxima
            presphase = r0(2); 
            pitch = 40;

            for rhat = 1:nRelements
                % Create a selection that wraps around the y direction
                thetasel = mod([ (presphase-pitch):(presphase+pitch)]-1, size(pdata,2))+1; %#ok<NBRAK>

                % Grab it and compute backwards finite difference
                checkarea = pdata(r0(1)+rhat- 1, thetasel);
                diff = abs(checkarea - circ_shift(checkarea,2,-1));

                diff(1) = 0; diff(end) = 0; % Kill the boundary cells since they don't represent
                % actual jumps

                % Set R
                R(rhat) = r0(1)+rhat-1; %#ok<AGROW>

                % Identify the location of the maximum jump
                md = 0; mdi = 1;
                for q = 1:2*pitch;
                    if diff(q) > md; md = diff(q); mdi = q; end
                end

            % debug stuff
            %	if rhat < 15
            %		figure(3); plot(diff); hold on; scatter(mdi, diff(mdi)); xxx = input('continue: ');
            %		fprintf('%3.3g ', mdi - pitch);
            %	end

                % Set the phase value
                Theta(rhat) = presphase - pitch + mdi ; %#ok<AGROW>
                presphase = Theta(rhat);
            end

            % Wrap around
            Theta = mod(Theta,4096);

        end
        
    end%STATIC
	
end%CLASS