function phi = getCylsymmPotential_h1(rho, dr, dz, ellipticTable)
% Computes the potential for an axisymmetric mass distribution
% rho(1,:) is the central axis, r=0
% Adopt convention such that first index is radial and second index is z
%

%[BoverA Eints] = generateEllipticTable(.002);

phi = zeros(size(rho));
subSel = (rho > 0);

%--- Get size of array and construct X/Z grids ---%
% X runs from 0 to (Nx-1) * h
% Z runs from -Nz*h/2 to +Nz*h/2
D = size(rho);
[R Z] = ndgrid( (0:(D(1)-1))*dr, ((-D(2)/2 ):(D(2)/2 - 1))*dz + dz/2);

da = dr * dz;

radrho = R.*rho;

for a = 1:D(1)
    for b = 1:D(2)/2
        aParameter = (Z - Z(a,b)).^2 + R.^2 + (a*dr - dr)^2; aParameter(aParameter == 0) = dr/2;
        bParameter = 2*R.*(a*dr - dr);

        aParameter = aParameter(subSel);
        bParameter = bParameter(subSel);
        % Table lookop of elliptical integrals
        ellipticIntegrals = interp1(ellipticTable.BoverA, ellipticTable.F, bParameter ./ aParameter) ./ sqrt(aParameter);
        % Integrate with g rho r dr dz
        phi(a,b) = sum(sum( radrho(subSel) .* ellipticIntegrals ));
    end
end
% mirror top-bottom symmetry
phi(:,(D(2)/2+1):end) = phi(:,(D(2)/2):-1:1);

phi = -da*phi;
   
end

