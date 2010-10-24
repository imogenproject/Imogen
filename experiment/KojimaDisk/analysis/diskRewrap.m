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
