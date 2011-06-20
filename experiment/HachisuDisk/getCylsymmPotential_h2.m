function phi = getCylsymmPotential_h2(rho, dr, dz, ellipticTable)
% Computes the potential for an axisymmetric mass distribution
% rho(1,:) is the central axis, r=0
% Adopt convention such that first index is radial and second index is z
%

%[BoverA Eints] = generateEllipticTable(.002);
phi = zeros(size(rho));

subSel = (rho > 0);

drho_dr        = (circshift(rho, [-1 0]) - circshift(rho, [1 0])) / (2*dr);
drho_dr(1,:)   = (-3*rho(1,:) + 4*rho(2,:) - rho(3,:)) / (2*dr);
drho_dr(end,:) = (3*rho(end,:) - 4*rho(end-1,:) + rho(end-2,:)) / (2*dr);
drho_dr = drho_dr(subSel);

%--- Get size of array and construct X/Z grids ---%
% X runs from 0 to (Nx-1) * h
% Z runs from -Nz*h/2 to +Nz*h/2
D = size(rho);
[R Z] = ndgrid( (0:(D(1)-1))*dr, ((-D(2)/2 ):(D(2)/2 - 1))*dz + dz/2);

da = dr * dz;
%nleft = 20;

for a = 1:D(1)
    for b = 1:D(2)/2
        
        %if ((a >= nleft) && (a < D(1))) 

	%	if (b > 1) && (b < D(2)); continue; end
	%end

        r0 = R(a,b);
        aParameter = (Z - Z(a,b)).^2 + R.^2 + r0^2; aParameter(aParameter == 0) = dr;
        bParameter = 2*R.*r0;

        aParameter = aParameter(subSel);
        bParameter = bParameter(subSel);

        % Table lookop of elliptical integrals
        E =     interp1(ellipticTable.BoverA, ellipticTable.F, bParameter ./ aParameter) ./ sqrt(aParameter);
        de_dq = interp1(ellipticTable.BoverA, ellipticTable.df_dq, bParameter ./ aParameter);

        dq_dr = 2*r0./aParameter - 2*bParameter .* R(subSel) ./ (aParameter.^2);

        % Integrate with g rho r dr dz
        phi(a,b) = sum( rho(subSel) .* E .* R(subSel) ...
                      + de_dq .* dq_dr .* R(subSel) .* rho(subSel) ./ sqrt(aParameter) * dr ...
                      + E .* ( - R(subSel).^2 .* rho(subSel) ./ aParameter + rho(subSel) + R(subSel).*drho_dr) * dr );
    end
end

phi(:,(D(2)/2+1):end) = phi(:,(D(2)/2):-1:1);

phi = -da * phi;

%rhoSubset = rho(nleft:(end-1),2:(end-1));

%rhoSubset(1,:) = rhoSubset(1,:) - phi(nleft-1,2:(end-1)) / (dr^2) + phi(nleft-1,2:(end-1)) ./ (2*dr*R(nleft-1,2:(end-1)));
%rhoSubset(end,:) = rhoSubset(end,:) - phi(end,2:(end-1)) / (dr^2) - phi(end,2:(end-1)) ./ (2*dr*R(end,2:(end-1)));
%rhoSubset(:,1) = rhoSubset(:,1) - phi(nleft:(end-1),1) / (dr^2);
%rhoSubset(:,end) = rhoSubset(:,end) - phi(nleft:(end-1),end) / (dr^2);

%rhoSubset = reshape(rhoSubset,[numel(rhoSubset) 1]);

%phiSubset = bicgstab(@(x) doCylLaplacian(phi, x, nleft, dr, R), 4*pi*rhoSubset, 1e-8, 1000);

%phi(nleft:(end-1),2:(end-1)) = reshape(phiSubset, [(D(1)-nleft) D(2)-2]);

%phi = -da * phi;
    
end

function lap = doCylLaplacian(phibord, phi, nleft, h, R)
phiDims = size(phibord) - [nleft 2];
phi = reshape(phi, phiDims);
phibord(:)=0;
phibord((nleft):(end-1),2:(end-1)) = phi;

lap = zeros(size(phibord));

lap = lap - 4*phibord;
lap = lap + circshift(phibord,[-1 0]);
lap = lap + circshift(phibord,[1 0]);
lap = lap + circshift(phibord,[0 1]);
lap = lap + circshift(phibord,[0 -1]);

lap = lap / (h^2);

lap = lap + (circshift(phibord, [-1 0]) - circshift(phibord,[1 0])) ./ (2*R*h);

lap = reshape(lap(nleft:(end-1),2:(end-1)), [numel(phi) 1]);
end
