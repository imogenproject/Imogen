function vanleerLimiter(flux, dLeft, dRight)
% This function uses the Van-Leer flux limiter to average out any non-monotonic artifacts in the
% input flux value and return an appropriate Total Variation Diminishing (TVD) result.
%
%>< flux     Array of current flux values.                                  FluxArray
%>> dLeft    Differences between left fluxVals.                             double(Nx,Ny,Nz)
%>> dRight   Differences between right fluxVals.                            double(Nx,Ny,Nz)

    signTest = max(dLeft .* dRight, 0) ./ (dLeft + dRight); % 1. Harmonic average.
    signTest(isnan(signTest)) = 0;                          % 2. Remove NaN.
    flux.array = flux.array + signTest;                   % 3. Impose monotonicity.

end
