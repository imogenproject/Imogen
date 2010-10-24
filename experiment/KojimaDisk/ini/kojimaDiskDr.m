function dr = kojimaDiskDr(q, radiusRatio, grid, edgePadding)
% Given the intrinsic properties q and radioRatio and the numerical quantities grid and
% edgePadding, calculates the appropriate grid step size to use for a run
%
%>> q           Angular momentum exponent. Omega = r^-q                        double
%>> radiusRatio The ratio of inner edge radius to density max radius           double
%>> grid        Grid resolution                                                double [3,1]
%>> edgePadding How much further to go past the outer edge                     double
%
%<< dr          The step size                                                  double

%--- Initialization ---%
info        = '';
xq          = 2 - 2*q; %Useful constant

%--- Find outer disk radius ---%
rout = kojimaFindOuterRadius(q, radiusRatio);
numRadialCells = floor(grid(1)/2);

dr = (rout + edgePadding) / numRadialCells;

end
