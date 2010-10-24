function result = multigridPotentialSolver(run, mass)
% This function computes gravitational potential using a multigrid coarsening algorithm
% Execution occurs in O(n^3 ln(n^3)) time on a 3d grid and O(n^2 ln(n^2)) time on a 2d grid.
%
%>< run         Data manager                                                    ImogenManager
%>< mass        Mass density                                                    FluidArray
%<< result      Gravitational potential                                         GravityManager

    %--- Compute coarsenings ---%
    %       Utilize pre-stored top-level arrays. This speeds by ~20% on 256^3 tests.
    [rhos poss] = massQuantization(mass.array, run.gridSize, run.DGRID, run.gravity.MG_TOPPOS);
    a           = run.DGRID;
    b           = run.gridSize;

    % This is a set of parameters required by the compiled (mg_bc) routine that define the volume it computes potential over
    % The following values will be correct for constant dx (though not necessarily equal in all directions)
    bvec(1) = a{1}*.5;
    bvec(2) = a{2}*.5;
    bvec(3) = a{3}*.5;

    bvec(4) = a{1}*b(1) - bvec(1);
    bvec(5) = a{2}*b(2) - bvec(2);
    bvec(6) = a{3}*b(3) - bvec(3);

    bvec(7) = b(1) - 1;
    bvec(8) = b(2) - 1;
    bvec(9) = b(3) - 1;

    %--- Approximate a cell's average distance from its center. ---%
    %       This matches the HOC4 stencil's point response and is correct even for non-uniform grid 
    %       spacing. The multigrid summer below IS NOT.
    self_pot_rad = sqrt(3) / sqrt(a{1}.^2 + a{2}.^2 + a{3}.^2);

    %--- Execute Multigrid Solver ---%
    %       Requires 4 arguments in this order:
    %           1. Cell array of hierarhcially coarsened density
    %           2. Cell array of hierarchially coarsened mass centers
    %           3. Realspace locations to sum at (9x1 vector)
    %           4. Coarsening algorithm's constant
    result = mg_bc_matlab(rhos, poss, bvec, sqrt(a{1}.^2 + a{2}.^2 + a{3}.^2) )
    
    if (run.gravity.constant ~= 1)
        result = run.gravity.constant*result;
    end
    
end
