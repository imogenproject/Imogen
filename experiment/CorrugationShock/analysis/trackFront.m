function frontX = trackFront(frame)
% Given a saveframe from Imogen, Identifies the shock's position.

% Take dRho/dx and find the indices of the max
d = diff(frame.mass, 1, 1);
[dRho_dx ind] = max(d, [], 1);

halfval = (frame.mass(end,1,1) + frame.mass(1,1,1))/2;

% Get rho in the cell before the max derivative
% Linearly extrapolate to where it's between equilibrium values
f0 = frame.mass(ind);
dx = (halfval - f0) ./ dRho_dx;

% Convert the linear indices to X indices and add the delta we foudn before
frontX = mod(ind, size(frame.mass, 1)) + dx;

end
