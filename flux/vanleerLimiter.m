function vanleerLimiter(flux, dLeft, dRight)
% This function uses the Van-Leer flux limiter to average out any non-monotonic artifacts in the
% input flux value and return an appropriate Total Variation Diminishing (TVD) result.
%
%>< flux     Array of current flux values.                                  FluxArray
%>> dLeft    Differences between left fluxVals.                             double(Nx,Ny,Nz)
%>> dRight   Differences between right fluxVals.                            double(Nx,Ny,Nz)



    if isa(dLeft,'GPUdouble') == 1
%	signTest = ( max(double(dLeft .* dRight), 0) ./ double(dLeft + dRight));
%	signTest(isnan(signTest)) = 0;
%	signTest = GPUdouble(signTest);

%        signTest = dLeft.*dRight;
%	cudaArrayAtomic(signTest, 0, ENUM.CUATOMIC_SETMIN);
%	signTest = double(signTest ./(dLeft+dRight));
%	signTest(isnan(signTest)) = 0;
%        signTest = GPUdouble(signTest);
%	signTest  = GPUdouble(); setReal(signTest); setSize(signTest, size(dLeft)); GPUallocVector(signTest);
	%signTest = zeros(size(dLeft), GPUdouble);
	signTest = fluxLimiterKernels(dLeft, dRight, 1);

%	sb = max(double(dLeft .* dRight), 0) ./ double(dLeft + dRight);
%	sb(isnan(sb))=0;
	%double(signTest)-2*sb

	flux.array = flux.array + signTest;                   % 3. Impose monotonicity.
    else
        signTest = max(dLeft .* dRight, 0) ./ (dLeft + dRight); % 1. Harmonic average.
        signTest(isnan(signTest)) = 0;                          % 2. Remove NaN.
        flux.array = flux.array + signTest;                   % 3. Impose monotonicity.
    end

end
