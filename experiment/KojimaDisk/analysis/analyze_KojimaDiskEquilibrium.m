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
deevee = delta{1}*delta{2}*delta{3};

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

indisk = (rho > rho(1,1));

T = .5 * sum(sum(sum(magp(indisk) ./ rho(indisk))));
% Potential energy = integral(-rho / r)
W = -sum(sum(sum(rho(indisk) ./ arrdist(indisk))));

toverw = abs(T / W);

%%% Compute drho/dt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[drhody drhodx] = gradient(rho, delta{1}, delta{2});
[dvdy zip] = gradient(momY ./ rho, delta{1}, delta{2});
[zip dvdx] = gradient(momX ./ rho, delta{1}, delta{2});

% Compute rho dot div(v) + v dot grad(rho)
drhodt = rho .* (dvdx + dvdy) + (momY ./ rho) .* drhody + (momX ./ rho) .* drhodx;

end

