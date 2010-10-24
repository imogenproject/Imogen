function [massGrid momGridX momGridY] = cyl2rect(r, rho, rhoV, griddim, deltaR)
% This function currently converts a line of dist, mass, and mom into a plane, which is the line 
% revolved around the z-axis. This function will convert a meridional slice of distance, mass 
% density and momentum density, into a 3D array which is the slice revolved about the z-axis. 
% Momentum is input in the phi direction, pos by right hand rule about the z-axis.
%
%>> r          Distance from the star's center                                 double  [1 nr]
%>> rho        Mass density at each r value                                    double  [1 nr]
%>> rhoV       Momentum density                                                double  [1 nr]
%>> griddim    Number of radial cells (= floor(x resoluion/2)                  double
%>> deltaR     Step size (uniform in x/y/z)
% 
%<< massGrid   Mass density array result (cell)                                double  [nx ny]
%<< momGridX   X momentum density array                                        double  [nx ny]
%<< momGridY   Y momentum density array                                        double  [nx ny]

%--- Prealocate grids and some constants. ---%
r0         = r(1);
massGrid   = zeros(2*griddim, 2*griddim);
momGrid    = zeros(2*griddim, 2*griddim, 2);

[i, j]     = ndgrid(1:2*griddim, 1:2*griddim); 

rr         = refinedGrid(i,j,r, griddim, deltaR);
u          = (rr == 0);

i          = deltaR * (i - griddim - .5)';
j          = deltaR * (j - griddim - .5)';

%--- Lathe mass and (scalar) momentum onto 2D grids ---%
massGrid   = pchip(r,rho,rr);
rrhoV      = pchip(r,rhoV,rr);

%--- Multiply by -cos and sin to get vector momentum ---%
rr(u)      = 1.0;
i(u)       = 0;
j(u)       = 0;

i          = i ./ rr;
j          = j ./ rr;

momGridX   = -i .* rrhoV;
momGridY   = j  .* rrhoV;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rr] = refinedGrid(i,j,r, nRadialZones, deltaR)
%    This function calculates the distance from the grid center to each
%    cell, then creates an array of indices associate with that cell
%    distance.

    r0    = r(1);

    rr = deltaR * sqrt( (i- nRadialZones - .5).^2+ (j- nRadialZones - .5).^2);

    test = (rr<=max(r));
    rr = test.*rr;
end
