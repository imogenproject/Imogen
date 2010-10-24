function ts = plotTimeSpectrum(stepby, theend, numzeros, DIR)
% Started in a run directory,
% 1. Loads next 2d slice
% 2. Transforms mass to polar coordinates and sums
% 3. Generates FFT
% 4. Clear slice's sx_... and goes to 1

stepnum = 0;

ts = zeros([round(theend/stepby + 1) 16]);

while stepnum <= theend
    % Get frame
    frmstring = '';
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