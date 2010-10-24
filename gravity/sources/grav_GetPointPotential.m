function phi = grav_GetPointPotential(grid, DGRID, pointpos, GM, softenRadius)
%===================================================================================================
% This function calculates the potential of a point mass
%===================================================================================================
% Example call: phi = grav_GetPointPotential(grid, dgrid, grid/2, 1, 8);
%===================================================================================================
% Variables:
% grid        double [3] - size of grid to solve over
% dgrid        double [3] or struct - grid step size
% pointpos    double [3] - cell index of point
%===================================================================================================

if any(numel(DGRID{:}) > 1)
    for i = 1:3
        if numel(DGRID{i}) == 1
        DGRID{i} = DGRID{i} * ones(grid);
        end
    end
end

%--- This can handle rectilinear coordinates (h(xi) = f(x)) or uniform step size ---%
if numel(DGRID{1}) > 1
    X = cumsum( DGRID{1}, 1 );
    Y = cumsum( DGRID{2}, 2 );
    Z = cumsum( DGRID{3}, 3 );

    gridr = floor(pointpos);
    delz = 1; if grid(3) == 1; delz = 0; end;

    % Linearly interpolate to find where point is
    posapprox = ([X(gridr(1)+1) Y(gridr(2)+1) Z(gridr(3)+delz)] - [X(gridr(1)) Y(gridr(2)) Z(gridr(3))]) .* (pointpos - gridr);
    posapprox = posapprox + [X(gridr(1)) Y(gridr(2)) Z(gridr(3))];

    X = X - posapprox(1);
    Y = Y - posapprox(2);
    Z = Z - posapprox(3);
else

    [X Y Z] = ndgrid(1:grid(1), 1:grid(2), 1:grid(3));
    X = (X - pointpos(1))*DGRID{1};
    Y = (Y - pointpos(2))*DGRID{2};
    Z = (Z - pointpos(3))*DGRID{3};
end

%--- Compute the potential and soften near-field ---%
%        First takes 1/R for all points. Then find the potential associated with the softening
%        radius and apply the piecewise-continuous softened parabolic potential
phi = -GM ./ sqrt(X.^2 + Y.^2 + Z.^2);

critPhi = - GM ./ softenRadius;

paramA = 1.5 * GM / (softenRadius);
paramB = GM^3 / (2*(softenRadius)^3);
tooclose = (phi < critPhi);
phi(tooclose) = -paramA + paramB * phi(tooclose).^-2;

end

