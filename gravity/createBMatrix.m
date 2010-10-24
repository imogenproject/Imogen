function result = BOperator(grid)

Nx         = grid(1);
    Ny         = grid(2);
    Nz         = grid(3);
    Nxy        = Nx*Ny;
    numCells = prod(grid);

    [X Y Z]    = ndgrid(1:Nx, 1:Ny, 1:Nz);
    indices    = reshape(1:numCells, grid);

    %--- Specify finite difference coefficients ---%
%    CEN = -128; % center cell
    FC  = 1/6;   % touching face of center cell
%    EC  = 3/128;    % touching edge of center cell
%    CC  = 1/128;    % touching corner of center cell

    %--- Begin creating the matrix ---%
    %       Each row in the laplacian matrix will nominally have 27 nonzero entries (less around the
    %       edges). To add to a matrix M of size i by j, coefficientMatrixSubset takes a set of i
    %       indices, a set of j indices, a constant to insert, an amount to shift by, and a set of
    %       locations to not act on (due to shifting ending up outside the grid) in that order.
    % Fourth argument: x_shift + y_shift + z_shift

%   result =     coefficientMatrixSubset(indices, indices, CEN,  0+0*Nx+0*Nxy, (X <= Nx ), numCells);

result =          coefficientMatrixSubset(indices, indices, FC,  1+0*Nx+0*Nxy, (X < Nx), numCells);
result = result + coefficientMatrixSubset(indices, indices, FC, -1+0*Nx+0*Nxy, (X > 1 ), numCells);

result = result + coefficientMatrixSubset(indices, indices, FC,  0+1*Nx+0*Nxy, (Y < Ny), numCells);
result = result + coefficientMatrixSubset(indices, indices, FC,  0-1*Nx+0*Nxy, (Y > 1 ), numCells);

result = result + coefficientMatrixSubset(indices, indices, FC,  0+0*Nx+1*Nxy, (Z < Nz), numCells);
result = result + coefficientMatrixSubset(indices, indices, FC,  0+0*Nx-1*Nxy, (Z > 1 ), numCells);

end
