function result = dumbPotentialSolver(run, mass)
% This routine calculates the potential at every point by directly summing rho/|r|
% DO NOT USE IT except for test comparison purposes, it executes in O(n^6) time for a 3d grid!
% I'm not even making this an option from imogen runfiles
% It would take a month to compute on a typical grid
%
%>< run         Data manager                                                    ImogenManager
%>< mass		Mass density                                                    FluidArray
%<< result      Gravitational potential                                         GravityArray

    grid     = mass.gridSize;
    glen     = prod(grid);
    allcells = 1:glen;
    [i j k]  = indexToijk(allcells, grid, 0);
    dgridMat = ones(glen, 1) * run.DGRID;

    % Might as well use loops on the outer parts, they're already down to seconds/iteration
    for u = 1:grid(1)
        fprintf('Dumb solver working on %i/%i\n', u, grid(1));
        for v = 1:grid(2) 
            for w = 1:grid(3)

                % For all cells in the grid, run a parllel calculation of the 1/r potential and sum
                pointMass = mass.array(u,v,w) * run.gravity.CONSTANT;

                % Create an x vector with the grid index minus the point position
                x = ( [i',j',k'] - ones(glen,1) * run.gravity.vars.test_pointGridPosition) .* dgridMat;

                % Get result at every point as mg/r and dump it back to a matrix
                gravPLine = gravPLine + pointMass ./ abs(x);

            end
        end
    end

    %--- Convert 1D Result to 3D Grid ---%
    result = convert1DLineto3D(gravPLine, grid, [], 0);

end
