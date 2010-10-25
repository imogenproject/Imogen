function run = fauxImogenManager(dGrid, gridSize)
% This function simulates the ImogenManager by creating default values needed for common ImogenArray
% operations.
%
%>> dGrid       Grid spacing cell array or 3 entry array for the data.          cell{3}/double(3)
%>> gridSize    Size of the grid in each spatial dimension.                     int(3)

    run.dGrid                = dGrid;
    run.gridSize             = gridSize;
    run.matlab               = ver('matlab');
    run.bc                   = BCManager.getInstance();
    run.bc.infinity          = 20;
    run.bc.modes.x           = 'circ';
    run.bc.modes.y           = 'circ';
    run.bc.modes.z           = 'circ';
    run.fades                = [];
    run.parallel.ACTIVE      = false;
    run.fluid.MASS_THRESHOLD = 1e-5;
    run.fluid.MINMASS        = 1e-5;

end
