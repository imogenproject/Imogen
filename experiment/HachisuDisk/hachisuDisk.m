function [rho, linMomDensity, dr, R, info, phi] = hachisuDisk(b, dtype, q, d, n, gridres, starMass, padding, tol, sgThreshold)

% Computes an axisymmetric self-gravitating structure utilizing the Hachisu '86 self-consistnet
% field method (HSCF) with a generalized rotational law.
% b           Inner radius; Outer rad fixed at 1 defines aspect ratio          double
% dtype       'spheroid' to place b on z axis, 'torus' to place b on r axis    string
% q, d        Define the rotational angular velocity w \prop (r^2 + d^2)^-q    double
% n           Polytropic index (gamma = n / (n-1))
% gridres     [Nx Ny Nz] cells. Disk rotational axis is z-hat.                 double [3]
% centerMG    Optionally places a point mass at the center of the grid         double
% padding     Amount of radial pad given to edge of torus                      double
% tol         Tolerance to reduce change per iteration by before returning     double
% sgThreshold Self gravitation threshold to be applied

% This function is implemented using the algorithm presented in Hachisu '86
% There are a few sign changes regarding h0^2. They do not affect the outcomes.

%--- Define basic variables: density, potential, enthalpy arrays ---%

rho = zeros([gridres(1) gridres(2)]);
phi = zeros(size(rho));
H   = zeros(size(rho));

a = 1; % b defines radius ratio. Set initial outer radius to 1.

dr = (a+padding) / gridres(1);
dz = dr; % for now
da = dr * dz;

[R Z] = ndgrid((0:(gridres(1)-1))*dr, ((-gridres(2)/2+.5):(gridres(2)/2 -.5))*dz );

indexA = round(.9/(1.1*dr));
indexB = round((b+.1)/(1.1*dr));

%--- Table of parameterized elliptic integrals for cylindrical potential solver ---%
ellipticTable = generateEllipticTable(.002);

%--- Generate initial guess density ---%
%        Time to convergence depends on the quality of this guess; These work OK ---%
if strcmp(dtype, 'spheroid')
    metric = sqrt(1 - (R/a).^2 - (Z/b).^2);
    metric(abs(imag(metric)) > 0) = -1;

    if starMass > 0
        error('HSCF ERROR: Spheroid cannot contain point mass!\n');
        return;
    end
end
if strcmp(dtype, 'torus')
    mu = .5*(a+b);
    ep = .5*(a-b);
    metric = sqrt(1 - ((R-mu)/ep).^2 - (Z/ep).^2);
end
rho(metric > 0) = metric(metric > 0);

phiExtern = zeros(size(rho));
if starMass > 0
    phiExtern = -starMass ./ sqrt(R.^2+Z.^2);
end

%--- Provide function to solve Omega integral ---%
%        We know the exact solution to all the radius-power-law curves considered in Hachisu '86:
%        omega^2 \prop (r^2 + d^2)^(-q)
switch q
    case 0; omegaIntegral = @(r) r.^2 / 2;
    case 1; omegaIntegral = @(r) .5 * log(r.^2 + d^2);
    otherwise; omegaIntegral = @(r) (1/(2-2*q))*(r.^2 + d^2).^(1-q);
end

potSolver = @getCylsymmPotential_h1;
diskMass = @(rho, R, dr) 2*pi*sum(sum(rho .* R)) * dr^2;

H = zeros(size(rho));
diffMetric(0,0,0,0,1); % Reset difference function
hasConverged = 0;
nIterations = 0;

% --- Run the HSCF loop on the dimensionless system (rhomax = r_outer = G = 1)
while hasConverged == 0;
    %--- 1: Solve gravitational potential ---%
    phi = potSolver(rho, dr, dz, ellipticTable) + phiExtern * diskMass(rho, R, dr);

    %--- 2: Compute h0 ---%
    % \Psi = \int r_{ax} \omega^2 dr_{ax}
    % omega^2 \prop (r^2 + d^2)^-q
    % h0^2 = [Phi(b) - Phi(a)] / [Psi(a) - Psi(b)]
    Phi_a = interp2(R', Z', phi', a, 0);
    if strcmp(dtype,'spheroid')
        Phi_b = interp2(R', Z', phi', 0, b);
        Psi_b = omegaIntegral(0);
    end;
    if strcmp(dtype, 'torus');
        Phi_b = interp2(R', Z', phi', b, 0);
        Psi_b = omegaIntegral(a*b);
    end
    h0sqr = (Phi_a - Phi_b) / (omegaIntegral(a) - Psi_b);

    %--- 3. Compute C ---%
    % C = Phi(a) + h0^2 Psi(a)
    C = Phi_a - h0sqr * omegaIntegral(a);

    %--- 4. Compute H ---%
    % H = C - Phi + h0^2 Psi
    H = C - phi + h0sqr * omegaIntegral(R);

    if strcmp(dtype,'torus') && (starMass > 0)
        H(R < .5*b) = -1;
    end
    
    %--- 5: Compute Rho from H and renormalize ---%
    % rho = (H / (K*(1+N)))^N = (H * gamma / (K*(gamma-1)))^(1/(gamma-1)
    rho = (H / (1+n));
    rho(rho < 0) = 0;
    rho = rho .^ n;
    rho = rho / max(rho(:));

    %--- 6: Test for self-consistency and break or goto 1 ---%
    diff = diffMetric(h0sqr, C, H, rho, 0);

    fprintf('iter %i: h0sqr = %g, C = %g, diff=%g\n', nIterations, h0sqr, C, diff);

    %--- Switches to slower 2nd-order solver near end ---%
%    if diff < 20*tol; potSolver = @getCylsymmPotential_h2; end
    
    if diff < tol; hasConverged = 1; end
    nIterations = nIterations + 1;
end

%--- Recompute potential one more time using final density distribution ---%
phi = potSolver(rho, dr, dz, ellipticTable) + phiExtern;

switch q
    case 0; linMomDensity = sqrt(h0sqr) * R .* rho;
    otherwise; linMomDensity = sqrt(h0sqr) * R .* rho .* (R.^2 + d^2).^(-q/2);
end
linMomDensity(isnan(linMomDensity)) = 0;
angMomDensity = linMomDensity .* R;

fprintf('Converged with tolerance %g in %i iterations. Dimensionless results:\n', diff, nIterations);
fprintf('Rb: %g		h0^2: %g\n', b, h0sqr);
fprintf('M : %g		V   : %g\n', 2*pi*da*sum(sum(rho.*R)), 2*pi*da*sum(sum(R .* (rho > 0) )) );
% compute J and T
T = linMomDensity.^2 ./ rho; T(isnan(T)) = 0;
fprintf('J : %g		T   : %g\n', 2*pi*da*sum(sum(R.*angMomDensity)), pi*da*sum(sum(R.*T)));
fprintf('-W: %g		\n', -2*pi*da*sum(sum(phi .* rho .* R)));

%--- Now we must convert to a physical system ---%
%        We choose to fix G and K at 1, which compels a relation between rho and Req:
%        Req = sqrt[ (1+n)rho^(1/n - 1) / max(H) ]
%        We normalize by setting disk mass to 1

epsilon = 2*pi*dr*dr*sum(sum(rho.*R));

tpower = 1/n - 1/3;
rhoScale = (epsilon^(-2/3)*max(H(:))/(1+n))^(1/tpower);
rScale = (1/(epsilon*rhoScale))^(1/3);

fprintf('dist scale: %g\tdensity scale: %g\n', rScale, rhoScale);
info = '';
R = R * rScale;
dr = dr * rScale;
dz = dz * rScale;
phi = phi * rScale^2 * rhoScale;
rho = rho .* rhoScale;

%--- Re-solve the rotation curve the easy but bad-at-the-edges way ---%
%        Because I'm too lazy to calculate how the rotation equation scales
P = rho.^(1+1/n);
linMomSqr = -(R.*rho.*(circ_shift(P,1,-1)-circ_shift(P,1,1)) + R.*(rho.^2).*(circ_shift(phi,1,-1)-circ_shift(phi,1,1)))/(2*dr);
linMomDensity = real(sqrt(linMomSqr));

%--- Choose 2 points and recover the analytical curve ---%
w1    = linMomDensity(indexA,gridres(2)/2) / (rho(indexA,gridres(2)/2)*R(indexA, gridres(2)/2));
w2    = linMomDensity(indexB,gridres(2)/2) / (rho(indexB,gridres(2)/2)*R(indexB, gridres(2)/2));

wrat  = (w2 / w1)^(2/q);
beta  = (R(indexB, gridres(2)/2)^2 * wrat - R(indexA, gridres(2)/2)^2) / (1-wrat);
alpha = w1 * (R(indexA, gridres(2)/2)^2 + beta^2)^(q/2);
linMomDensity = alpha * R .* rho .* (R.^2 + beta^2).^(-q/2);

%--- Applies a diffusion operator to soften speed gradient ---%
speed = linMomDensity ./ rho;
speed(isnan(speed)) = 0;
diffme = (speed ~= 0);

rhomaxidx = find(rho == max(rho(:)));
fprintf('MIRP = %g\n', 2*pi*R(rhomaxidx) / speed(rhomaxidx));

%for q = 1:2;
%    speed = hscf_diffuse(speed,.25);
%end
linMomDensity = smoothVelocity(rho, linMomDensity);

rho(speed > 0) = max(rho(speed > 0), 1.5*sgThreshold*rhoScale);
%linMomDensity = rho .* speed;

end

% Apply a diffusion operator
function result = hscf_diffuse(f, D)
lap = constant_shift(f,1,1) + constant_shift(f,1,-1) + constant_shift(f,2,1) + constant_shift(f,2,-1) - 4*f;
result = f + D*lap;
end

function result = diffMetric(h0, C, H, rho, doreset)

persistent h0_old;
persistent C_old;
persistent H_old;
persistent rho_old;

if isempty(h0_old) || doreset == 1
    h0_old = 1;
    C_old = 1;
    H_old = zeros(size(rho));
    rho_old = zeros(size(rho));
end

result = max([abs(h0/h0_old - 1) abs(C / C_old - 1) abs(max(max(H - H_old))) abs(max(max(rho - rho_old))) ]) ; % lol

h0_old = h0;
C_old = C;
H_old = H;
rho_old = rho;

end

