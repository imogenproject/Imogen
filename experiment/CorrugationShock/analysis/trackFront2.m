function frontX = trackFront(mass, x)
% Given a saveframe from Imogen, Identifies the shock's position.

% Take dRho/dx and find the indices of the max
d = diff(mass, 1, 1);
[dRho_dx ind] = max(d, [], 1);

halfval = (mass(end,1,1) + mass(1,1,1))/2;

% Get rho in the cell before the max derivative
% Linearly extrapolate to where it's between equilibrium values
f0 = mass(ind);
dx = (halfval - f0) ./ dRho_dx;

% Convert the linear indices to X indices and add the delta we found before
a = x(mod(ind, size(mass, 1)))';
b = x(mod(ind, size(mass, 1))+1)';

frontX = a + (b-a).*dx;

end
