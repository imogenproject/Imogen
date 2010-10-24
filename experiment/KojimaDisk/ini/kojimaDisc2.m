function [mass, momX, momY, dR, info] = kojimaDisc2(q, GAMMA, radiusRatio, grid, GRAVITY, ...
                                          starMass, edgePadding, pointRadius, minDensityCoeff, momregions, useZMirror)
% Generates an equilibrium differentially-rotating disc in 2D or 3D.
% My convention regarding the cell-centered vs corner-centered question while writing this is such:
% Consider the cell indexed by grid(1,1,1) as filling the space (0,0,0) x (dx,dy,dz);
% then it has center (dx/2, dy/2, dz/2)
% Then with matlab's indexing starting at 1,
% an arbitrary cell's center is ([grid index] - [.5, .5, .5]) .* [dx, dy, dz]
% With index-from-0 like C the center is (idx + .5)*dR instead.
%
% 
%
% Original test conditions:  q = 2, GAMMA = 1.5, radiusRatio = .75,
%                            grid = [128, 128, 32], GRAVITY = 1, starMass = 1
%
%>> q               angular velocity exponent; omega = omega_0 * r^-q.                  double
%>> GAMMA           polytropic index used for equation of state.                        double
%>> radiusRatio     inner radius / radius of density max.                               double
%>> grid            [x y z] size of the grid.                                           int(3)
%>> GRAVITY         gravitational constant in units desired.                            double
%>> starMass        specified mass of central star for normalization.                   double
%>> edgepadding     The distance between the outer edge of the disk and the grid        double
%                   edge in polytrope units, typically perhaps .2
%>> pointRadius     Distance from the star at which the point potential is softened     double
%>> minDensityCoeff What fraction of the density max to assign to zero-mass regions     double
%>> momregions      4 element array that assigns momentum distributions

%Return values:
%<< mass            equilibrium mass density.                                           double(GRID)
%<< momX            equilibrium x momentum density.                                     double(GRID)
%<< momY            equilibrium y momentum density.                                     double(GRID)
%<< dR              grid cell dimensions for proper scaling of the model.               double
%<< info            information string containing specifics about the run.              str

    %___________________________________________________________________________
    %--- Initialization ---%
    mass = zeros(grid);
    momX = zeros(grid);
    momY = zeros(grid);

    info        = '';
    discMass    = 0;
    xq          = 2 - 2*q; %Useful constant

    if useZMirror ~= 1; useZMirror = 0; end

    if momregions == 0
        % Default to zero momentum except in the disk region
    % This becomes clear around line 125
    fprintf('Defaulting on momentum distribution.\n');
        momregions = [1 1 2 1];
    end

    pointMG           = GRAVITY * starMass;
    rdinf.ptrad       = pointRadius;

    %--- Find outer disk radius ---%
    %       Computed by solving for density equal to zero numerically
    rdinf.rin         = radiusRatio;
    rdinf.rout        = kojimaFindOuterRadius(q, rdinf.rin);
    rdinf.rmax        = rdinf.rout + edgePadding;

    % Let the number of radial cells be half the grid size, rounded down
    rdinf.nRadialZones = floor(grid(1)/2);

    % Calculate the cell indices of the radial zones
    rdinf.ptradIndex  = round(rdinf.nRadialZones * rdinf.ptrad/rdinf.rmax);
    rdinf.rinIndex    = round(rdinf.nRadialZones * rdinf.rin / rdinf.rmax);
    rdinf.routIndex   = round(rdinf.nRadialZones * rdinf.rout/ rdinf.rmax);

    strcat(info, sprintf('Softening, inner, outer, edge radii: %g | %g | %g | %g\n', rdinf.ptrad, rdinf.rin, rdinf.rout, rdinf.rmax));

    %--- Determine the grid step size and let it be uniform ---%
    rdinf.dr = kojimaDiskDr(q, rdinf.rin, grid, edgePadding);
    strcat(info, sprintf('Determined grid cell spacing to be: %.8f\n', rdinf.dr));

    %--- Calculate the cell-centered distances ---%
    %        axialRadius is the polar coordinate radius looking down on the top of the disk
    %        centerRadius is the spherical coordinate radius from the point mass
    %        If useZMirror is 1, we're simulating only the top half of a disk
    [X Z] = ndgrid(1:2*rdinf.nRadialZones, 1:grid(3));

    X = X - .5;

    if useZMirror == 1
        Z = Z - .5;
    else
        Z = Z - floor(grid(3)/2) - .5;
    end

    rdinf.axialRadius  = X * rdinf.dr;
    rdinf.centerRadius = sqrt(X.^2 + Z.^2) * rdinf.dr;
    clear X Z;

    %___________________________________________________________________________
    %--- Calculation for an azimuthal slice of a disk (2D) ---%

    %--- Constants of integration ---%
    %       hsq is useless here because of our unit scaling.
    %       c1 is the integration constant from solving the Bernoulli equation.

%%% EXPERIMENTAL CODE
sharpness = 30;
v = @(r) .5*(exp(sharpness*(r-rdinf.rin)).*(r.^(1-q)) + exp(sharpness*(rdinf.rin-r)).*(r.^(-.5)) ) ./ cosh(sharpness*(r-rdinf.rin));


    hsq                = xq * (1/rdinf.rout - 1/rdinf.rin)/(rdinf.rin^xq - rdinf.rout^xq);
    c1                 = -1/rdinf.rin - (hsq/xq)*(rdinf.rin^xq);

    peakCenterDensity  = findMaxDensity(q, rdinf.rin, GAMMA);
    lomass             = minDensityCoeff * peakCenterDensity;

    %--- Calculates the pressure integral (Bernoulli equation), clamps #s < 0, solves for rho ---%
    bernoulli          = c1 + (1 ./ rdinf.centerRadius) + (hsq / xq) * rdinf.axialRadius.^xq;
    isPartOfDisk       = (bernoulli > 0);
    bernoulli(~isPartOfDisk) = 0;

    rho                = max((bernoulli * (GAMMA-1)/GAMMA).^(1/(GAMMA-1)), lomass);

%%% EXPERIMENTAL CODE

for q = 1:numel(rho); rho(q) = .4*quad(@(r) (-r.^(-2) + (v(r).^2)./r), 1,rdinf.centerRadius(q)) + rhomax^(2/3); end;
%isPartOfDisk = (real(rho) > 0);

isPartOfDisk = (rdinf.centerRadius < 1.9) & (rdinf.centerRadius > .6);
rho = max(real(rho.^(3/2)), lomass);

innerEdge = (rdinf.centerRadius < rdinf.rin) & (rdinf.centerRadius > .6);
rho(innerEdge) = max(rho(innerEdge), 7*lomass);



    borderRegion       = simpleBlur(isPartOfDisk * 1.0, .5, 4,1);
    borderRegion       = (borderRegion < 1) & (borderRegion > 0);
  
 %   rhoblur            = simpleBlur(rho, .4, 4, 1);
%    rho(borderRegion)  = rhoblur(borderRegion);

    %--- Compute momentum distribution ---%
    %        Apply Kojima momentum to the vast bulk of the disk, apply a blurred function to the
    %        inner and outer perimeters.
%    momentumKojima = rdinf.axialRadius.^(1-q) .* rho;
    momentumKojima = v(rdinf.axialRadius) .* rho;
    momentumKepler = rdinf.axialRadius.^(-.5) .* rho;

    sliceMomentum(isPartOfDisk)  = momentumKojima(isPartOfDisk);
    sliceMomentum(~isPartOfDisk) = momentumKepler(~isPartOfDisk);
    sliceMomentum                = simpleBlur(sliceMomentum, .2, 3, 1);

    mom = zeros(size(rho));
    mom(isPartOfDisk) = momentumKojima(isPartOfDisk);
%   mom(borderRegion) = sliceMomentum(borderRegion);

    %--- Lathe radial info onto cubical grid --- %
    %        I experimented with doing this in one interpolation with not much luck.
    for zct = 1:grid(3)
        [mass(:,:,zct) momX(:,:,zct) momY(:,:,zct)] = cyl2rect(rdinf.axialRadius(:,zct), rho(:,zct), mom(:,zct), grid(1)/2, rdinf.dr);
    end

    fprintf('Ratio of Disk/Star mass from radial integration found to be: %g\n', ...
                              (1+useZMirror)*2*pi*sum(sum(rdinf.axialRadius.*rho))*rdinf.dr*rdinf.dr);
    fprintf('Radio of Disk/Star mass from grid summation found to be: %g\n', ...
                              (1+useZMirror)*sum(sum(sum(mass)))*rdinf.dr^3);

dR = rdinf.dr;

end

% This is supposed to vertically integrate a disk's density
% So it can be pancaked to run in 2d but that doesn't work
function Psigma = verticallyIntegrateRho(r, q, rin, gamma)
Psigma = zeros(size(r));

xq = 2-2*q;
eta = (gamma-1)/gamma;
nu = 1/(gamma-1);
c1 = 1/rin + (rin^xq)/xq;

for g = 1:numel(r)
        r0 = r(g);
        anonrho = @(z)  (( (r0.^xq)/xq + 1./sqrt(r0.*r0+z.*z) - c1)*eta).^nu;
        h0 = real(sqrt( ((r0^xq)/xq - c1)^-2 - r0^2));

    if (r0 > rin) && (h0 > 0)
        Psigma(g) = 2*quadgk(anonrho, 0, h0);
    else
        Psigma(g) = 0;
    end
end

Psigma = real(Psigma);
% Outside the valid disk radii mass density and by extension its integral become complex.
% We try to avoid this as much as possible but near edges we typically end up with a small (1e-9 like) imaginary part which is nonphysical.

end

%--- Given the 3 parameters that characterize a disk, returns the extact max density ---%
function rhomax = findMaxDensity(q, rin, gamma)
xq = 2-2*q;
nu = 1/(gamma-1);

rhomax = ((1 - 1/gamma)*(1/xq + 1 - (rin^xq)/xq - 1/rin))^nu;
end

