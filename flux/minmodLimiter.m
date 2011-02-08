function minmodLimiter(flux, dLeft, dRight)
% This function uses the MinMod flux limiter to average out any non-monotonic artifacts in the
% input flux value and return an appropriate Total Variation Diminishing (TVD) result.
%
%>< flux     Array of current flux values.                                  FluxArray
%>> dLeft    Differences between left fluxVals.                             double(Nx,Ny,Nz)
%>> dRight   Differences between right fluxVals.                            double(Nx,Ny,Nz)

    correction = 0.125*(sign(dLeft) + sign(dRight)) .* min(abs(dLeft), abs(dRight));
    flux.array = flux.array + correction;
end
