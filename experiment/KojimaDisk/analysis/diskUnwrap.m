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
