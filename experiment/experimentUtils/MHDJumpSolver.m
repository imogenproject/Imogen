function result = MHDJumpSolver(ms, ma, theta, GAMMA)
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
    %            pre: pr, rho, ms, ma, vz=0, bz=0, tang=tan(theta) gamma=5/3
    %           post: BX = bx (divergance eqn), VZ=0, BZ=0
    %   Uknowns:
    %            pre: vx, vj, bx, bj
    %           post: RHO, VX, VJ, PR, BJ 
    
    % Translate these into names used in Mathematica formulation
    vxpre = ms * sqrt(5/3);
    vypre = vxpre * tan(theta);
    
    bx = vxpre / ma;
    bypre = bx * tan(theta);
    
    rhopre = 1;
    Ppre = 1;
    g = GAMMA;
    tt = tan(theta);
    
    % Solve for postshock vx: Use quartic solver
    a4 = (rhopre^3 *vxpre^3 + g *rhopre^3 *vxpre^3);
    a3 = (-2 *bx^2 *rhopre^2 *vxpre^2 - 2 *bx^2 *g *rhopre^2 *vxpre^2 - 2 *g *Ppre *rhopre^2 *vxpre^2 - 2 *bx^2 *g *rhopre^2 *tt^2 *vxpre^2 - 2 *g *rhopre^3 *vxpre^4);
    a2 = (bx^4 *g *rhopre* vxpre + 4 *bx^2 *g *Ppre *rhopre *vxpre + 3 *bx^4 *g *rhopre *tt^2 *vxpre + bx^4 *rhopre *(1 + tt^2)* vxpre + 4* bx^2* g *rhopre^2 *vxpre^3 + 2 *g *Ppre* rhopre^2 *vxpre^3 - 2 *bx^2 *rhopre^2 *tt^2 *vxpre^3 + 2 *bx^2 *g *rhopre^2 *tt^2* vxpre^3 - rhopre^3 *vxpre^5 + g* rhopre^3 *vxpre^5);
    a1 = (-2 *bx^4 *g *Ppre - 2 *bx^4* g *rhopre *vxpre^2 - 4 *bx^2* g *Ppre* rhopre* vxpre^2 - 4* bx^4 *g *rhopre *tt^2 *vxpre^2 + 2 *bx^2* rhopre^2* vxpre^4 - 2 *bx^2 *g *rhopre^2 *vxpre^4 + 2 *bx^2 *rhopre^2 *tt^2 *vxpre^4);
    a0 = (2 *bx^4 *g *Ppre *vxpre + bx^4 *g *rhopre *vxpre^3 + bx^4 *g *rhopre *tt^2* vxpre^3 - bx^4 *rhopre* (1 + tt^2) *vxpre^3);
    
    vpost = solveQuartic(a4, a3, a2, a1, a0);
    vpost = vpost(imag(vpost) == 0); % Two real solutions should exist
    vxpost = min(vpost); % The lesser is the one containing a discontinuity; The other corresponds to no jump.
    
    bypost = (-bx^2*bypre + bypre*rhopre*vxpre^2)/(rhopre*vxpost*vxpre - bx^2);
    vypost = (bx*bypost - bx*bypre + rhopre*vxpre*vypre) / (rhopre*vxpre);
    Ppost = -bypost^2 + bypre^2 + Ppre - rhopre *vxpost *vxpre + rhopre *vxpre^2;
    rhopost = rhopre *vxpre / vxpost;

    result.mass       = rhopost;
    result.pressure   = Ppost;
    result.velocity   = [vxpost vypost];
    result.magnet     = [bx bypost];
    result.theta      = theta;
    result.sonicMach  = ms;
    result.alfvenmach = ma;

end
