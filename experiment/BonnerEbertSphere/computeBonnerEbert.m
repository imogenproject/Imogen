function [R rho ener] = computeBonnerEbert(rho0, GAMMA, dr, rmax, k, rhobg)
% This function computes with extreme accuracy the mass distribution of any hydrostatic
% self-gravitating isothermal gas sphere obeying an equation of state P = k rho^gamma, which
% includes the ideal gas (gamma=1). It returns arrays of radii and distance.
%
%>> rho0     Density of the center of the globule                              double
%>> GAMMA    Pressure scaling exponent                                         double
%>> dr       Radial step size                                                  double
%>> rmax     Integrates the equation to this radius                            double
%>> k        Pressure scaling constant; Related to temperature                 double
%>> rhobg    The background density                                            double
%
%<< R        Distance array                                                    double [1xn]
%<< rho      Density array                                                     double [1xn]
%<< ener     Energy array                                                      double [1xn]

%--- These anonymous functions are used by the Adams-Bashforth integrator ---%
f = @(MASS, RHO, RADIUS) -MASS * RHO^(2-GAMMA) / (k * RADIUS^2 * GAMMA);
dm = @(RHO0, RHOPRIME, R0) 4*pi*( ((RHO0 - RHOPRIME*R0)*((R0+dr)^3 - R0^3))/3 + ...
                                              (RHOPRIME*((R0+dr)^4 - R0^4))/4 );

%--- Some startup checks ---%
%        The convergence parameter measures how near the "defined" steps come to the point
%        where the Picard iterations unravel
if nargin < 6
    rhobg = rho0 / 10000;
end

convergenceParameter = exp(.5*log(k/rho0) - 1)/(5*dr);

fprintf('Convergence parameter: %g\n', convergenceParameter);
if convergenceParameter < 1;
    if convergenceParameter < .5;
        fprintf('WARNING: Convergence parameter is < .5. Results will most likely be disastrous.\n');
    else
        fprintf('WARNING: Convergence parameter is < 1. Results are not trustworthy.\n');
    end
end

%--- Evaluate the first few steps using Picard series ---%
%        Expansion of density to 4th order is calculated and used to seed the Adams-Bashforth
%        integrator routine. Extreme accuracy here is _crucial_ to obtaining precision in the
%        final distribution. These remain converged until just past the 'knee'
beta = pi/(k*GAMMA);
theta = 2-GAMMA;
alpha = beta * rho0^theta;

for q = 1:6
x      = q - 1;
R(q)   = x*dr;

% Second order accurate for small r, gauranteed right
%rho(q) = rho0;
%M(q) = 4*pi*R(q)^3/3;

% Fourth order accurate for small r, gauranteed right
%rho(q) = rho0 * (1 - (2/3)*alpha*R(q)^2);
%M(q)   = 4*pi*rho0*(R(q)^3/3 - .2*alpha*R(q)^5);

% Sixth  order accurate for small r, something fishy but it works
rho(q) = rho0 * (1 - (2/3)*alpha*R(q)^2 + 4*(8/3)*(theta/12 - .05)*alpha^2*R(q)^4);
M(q)   = 4*pi*rho0*(R(q)^3/3 - .2*alpha*R(q)^5 + (32/21)*(theta/12 - .05)*alpha^2*R(q)^7);
end

%--- BEGIN NORMAL NUMERICAL INTEGRATION ---%
x = q+1;
hasReachedZeroDensity = 0;

while (R(end) < rmax) && (hasReachedZeroDensity == 0)
    % Predict next density using 5th order Adams-Bashforth
    % This is a high-order routine and _will_ destabilize if the convergence parameter is bad!
    R(x) = (x-1)*dr;
    rho(x) = rho(x-1) + dr*( (1901/720)*f(M(x-1), rho(x-1), R(x-1)) ...
                            -(1387/360)*f(M(x-2), rho(x-2), R(x-2)) ...
                            +(109/30)  *f(M(x-3), rho(x-3), R(x-3)) ...
                            -(637/360) *f(M(x-4), rho(x-4), R(x-4)) ...
                            +(251/720) *f(M(x-5), rho(x-5), R(x-5)) );
    
    % Integrate mass; straightforward 3nd order piecewise-linear integration
    % FIXME: Improve this
    if (real(rho(x)) < rhobg) || (abs(imag(rho(x))) > 0)
        hasReachedZeroDensity = 1;
        rho(x) = rhobg;

        rho(x+1) = rhobg;
        R(x+1) = rmax;
    end

    M(x) = M(x-1) + dm(rho(x-1), (rho(x)-rho(x-1))/dr, R(x-1));
    x = x + 1;

end

R(1) = 0;
fprintf('Isothermal self-gravitating sphere with gamma=%g, k=%g computed: R = %g, Total mass = %g\n', GAMMA, k, R(end), M(end));

fprintf('Rcrit = %g, Pcrit = %g, Pouter = %g\n', .411*M(end)/k, 1.4*k^4/M(end)^2, k*rhobg);

if hasReachedZeroDensity == 0
    fprintf('WARNING: Sphere did not reach background density before integration ended!\n');
end

ener = k * rho;

end

