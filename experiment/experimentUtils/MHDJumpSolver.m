function [result allvx] = MHDJumpSolver(ms, ma, theta, GAMMA)
% Solves the 3D MHD Jump Conditions for an equilibrium shock wave assuming that the pre-shock region
% conforms to zero z components for both velocity and magnetic field (this
% can always be acheived by a frame transform) as well as unity values for 
% mass density and pressure in the pre-shock region. Gamma has also been set to 5/3 in the equations to 
% simplify the numeric complexity, but that could be changed if necessary. Once the solution of the
% system has been found, the solver looks through the numerical results structure, determines which
% solution is the physical one and returns the data structure for that solution. It will also save
% a mat file with the solution for subsequent loading in the Corrugation shock initializer.
%
%><  run      Manager object for current solver            LinogenManager object
%><  linVars  Prototype data class for solver              LinogenVariables object

    %--- Define jump conditions ---%
    %   Knowns:
    %            pre: pr, rho, ms, ma, vz=0, bz=0
    %           post: BX = bx (divergance eqn), VZ=0, BZ=0
    %   Uknowns:
    %            pre: vx, vj, bx, bj
    %           post: RHO, VX, VJ, PR, BJ 

    theta = theta * pi / 180; 

    % Solve the preshock flow in terms of the 3 free quantities
    vxpre = ms * sqrt(GAMMA);
    vypre = vxpre * tan(theta);
    
    bx = vxpre / ma;
    bypre = bx * tan(theta);
    
    rhopre = 1;
    Ppre = 1;
    g = GAMMA;
  
    pxpre = rhopre*vxpre;
    txpre = rhopre*vxpre^2;

    % The following equations could be rewritten in terms of the 3 free quantities only (since they
    % completely determine the flow) but writing in terms of physical quantities gives some idea where
    % the terms come from.

    % These define a quartic equation sum(vxpost^n a_n) = 0 that determines vxpost
    a4 = (1+g)*pxpre^3;
    a3 = -pxpre^2*(2*bx^2*(g+1) + g*(bypre^2 + 2*Ppre + 2*txpre));
    a2 = pxpre*(bx^4*(g+1) + txpre*(2*bypre^2*(g-1) + 2*g*Ppre +(g-1)*txpre) + bx^2*(bypre^2*(g+1) + 4*g*(Ppre+txpre)));
    a1 = -bypre^2*(g-2)*txpre^2 - 2*bx^4*g*(Ppre+txpre) - 2*bx^2*txpre*(bypre^2*g + 2*g*Ppre +(g-1)*txpre);
    a0 = bx^2*vxpre*(bypre^2*(g-1)*txpre + bx^2*(2*g*Ppre +(g-1)*txpre));
    vpost = solveQuartic(a4, a3, a2, a1, a0);
    allvx = vpost; % For those who wish to examine

    vpost = real(vpost(imag(vpost) < 1e-11)); % There is a potential numerical instability in the equation solver that this avoids.

    vxpost = min(vpost); % The lesser is the one containing a discontinuity; The other corresponds to no jump and the complex options are unphysical.
    bypost = (-bx^2*bypre + bypre*rhopre*vxpre^2)/(rhopre*vxpost*vxpre - bx^2);
    vypost = (bx*bypost - bx*bypre + rhopre*vxpre*vypre) / (rhopre*vxpre);
    Ppost = .5*(-bypost^2 + bypre^2 + 2*Ppre - 2*rhopre*vxpost*vxpre + 2*rhopre*vxpre^2);
    rhopost = rhopre *vxpre / vxpost;

    result.mass       = [1; rhopost];
    result.pressure   = [1; Ppost];
    result.velocity   = [vxpre vxpost; vypre vypost; 0 0;];
    result.magnet     = [bx bx; bypre bypost; 0 0];
    result.theta      = theta;
    result.sonicMach  = ms;
    result.alfvenMach = ma;

    err = evalRankineHugoniotConditions(rhopre, vxpre, vypre, bx, bypre, Ppre, rhopost, vxpost, vypost, bypost, Ppost, GAMMA);

%    fprintf('Norm for obediance of Rankine-Hugoniot equations (lower=better): %g', norm(err));

end

function f = evalRankineHugoniotConditions(rho1, vx1, vy1, bx, by1, P1, rho2, vx2, vy2, by2, P2, g)

f(1) = vx2*by2 - vy2*bx - vx1*by1 + vy1*bx;
f(2) = rho2*vx2 - rho1*vx1;
f(3) = rho2*vx2*vx2 + P2 + by2*by2/2 - rho1*vx1*vx1 - P1 - by1*by1/2;
f(4) = rho2*vx2*vy2 - bx*by2 - rho1*vx1*vy1 + bx*by1;
f(5) = .5*rho2*(vx2^2+vy2^2)*vx2 + g*P2*vx2/(g-1) + by2*(vx2*by2-bx*vy2) -  .5*rho1*(vx1^2+vy1^2)*vx1 - g*P1*vx1/(g-1) - by1*(vx1*by1-bx*vy1);

f=f';

end

