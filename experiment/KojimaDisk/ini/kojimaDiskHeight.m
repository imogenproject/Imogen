function [ hmax nzmin ] = kojimaDiskHeight(q, rin, gamma, grid, padding)
% This function finds the minimum number of Z zones needed given uniform dx = dy = dz with no
%padding.
%
%>> q         Angular momentum exponent                                        double
%>> rin       Ratio of inner radius to density max radius                      double
%>> gamma     Pressure P = rho^(gamma)                                         double
%>> grid      Grid resolution                                                  double [3,1]
%>> padding   Amount of padding on the grid's edge                             double

%--- Compute useful factors ---%
xq = 2-2*q;
eta = (gamma-1)/gamma;
nu = 1/(gamma-1);
c1 = 1/rin + (rin^xq)/xq;

%--- Find height of disk exactly ---%
%        For details, see Kojima model papers
hmax = real(sqrt( (1/xq - c1)^-2 - 1));

%--- If a radius and amount of padding are given, return a number of cells too ---%
nzmin = [];
if nargin > 3
	delr = kojimaDiskDr(q, rin, grid, padding);
	nzmin = 2 * round(hmax / delr);
end

end
