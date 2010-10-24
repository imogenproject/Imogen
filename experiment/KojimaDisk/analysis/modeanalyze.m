function result = modeanalyze(radialIndex, optString)
% This function sequentially loads disk MAT files from a list, unwraps them, and performs a mode
% spectrum analysis.
%>> radialIndex: The radial cell # to examine modes at; set to zero for radial averaging
%<< modeinfo   : The 12-by-N returned information on the first 12 modes

if ischar(radialIndex)
    if lower(radialIndex) == 'help'
        fprintf('Mode analysis help: Argument 1 can be either\n	[x] to examine radial index x,\n	[x y] to examine a range\n	a cell array of [x] and [x y].\n\nIf the second argument is "eigenfunction", modeanalyze will examine regularly spaced bins and argument 1 should be either\n	[x y] to look at (y-x) bins of size 1\n	[x s y] to look at bins of size s stepping from x to y.\n');
    end
    return;
end

% Get a range of mat files to analyze
inr = input('.mat range to load: ', 's');
inr = sprintf('irange = %s;', inr);
eval(inr);

filenamepad = input('How many numbers in filename pad? ');

X = numel(irange);

switch nargin;
    case 0; fprintf('Averaging over all R; Normalizing by m=0 @ t=0\n'); radialindex = 0;
    case 1;
    case 2;
        if optString == 'eigenfunction'

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

% Get storage for all the modes we will be analyzing
result.yaxis = [1:12];
result.xaxis = zeros(1,X);

if ~iscell(radialIndex)
    q = radialIndex;
    radialIndex = cell(1);
    radialIndex{1} = q;
end

result.modeinfo = zeros(12,X,numel(radialIndex));
result.growthrate = zeros(12,2,numel(radialIndex));

for x = 1:X
    s = sprintf('%i ', x);
    fprintf('%s ', s);
    % Load a mat file and move it to our known structure name
    evstr = ['2D_XY_%' num2str(filenamepad) '.' num2str(filenamepad) 'i.mat'];
    if irange(x) == 0
        load 2D_XY_START.mat;
        sxcur = sx_XY_START;
        clear sx_XY_START;
    else
        evstr = sprintf(evstr, irange(x));
        eval(['load ' evstr]);   

        evstr = ['sx_XY_%' num2str(filenamepad) '.' num2str(filenamepad) 'i'];
        evstr = sprintf(evstr, irange(x));
    
        eval(['sxcur = ' evstr ';']);
        eval(['clear ' evstr]);
    end

    result.xaxis(x) = sum(sxcur.time.history) / (2*pi); % Normalize to chirps
    
    % Perform polar transform, rescale by radius while radially summing
    mpolar = diskUnwrap(sxcur.mass);

    %--- Calculate radius vector; Rescale z axis appropriately if it exists ---%
    if x == 1;
        origrad = size(sxcur.mass, 1) * sxcur.dGrid{1} * .5;
        newdr = origrad / size(mpolar, 1);
        newr = cumsum(ones(1, size(mpolar,1)) * newdr);
        if isfield(result, 'zaxis')
            result.zaxis = result.zaxis * sxcur.dGrid{1} * .5;
        end
    end

    %--- Perform the Fourier analysis on all requested bins ---%
    for y = 1:numel(radialIndex)
        [result.modeinfo(:,x,y) nc] = fftSlice(mpolar, newr, radialIndex{y});
        if x == 1; result.nConst(y) = nc; end
    end

end

for y = 1:numel(radialIndex);
    result.modeinfo(:,:,y) = result.modeinfo(:,:,y) ./ result.nConst(y);

    if numel(find(result.xaxis >= 1)) > 2
        for m = 1:12;
            result.growthrate(m,:,y) = findGrowthRate(result.xaxis, result.modeinfo(m,:,y), 1);
        end
    else
        fprintf('Insufficient time elapsed to guess at growth rate\n');
    end
end

end

%--- Takes unwrapped disk, radius vector and slice selection; Returns FFT ---%
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

%--- Given time and amplitude, calculate growth rate characteristic time ---%
function result = findGrowthRate(t, a, tmin)

x = numel(find(t < tmin))+1;

t = t(x:end);
a = a(x:end);

fcn = polyfit(t, log(a), 1);

result = fcn;
end
