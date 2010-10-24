function result = createLaplacianMatrix(grid)
% This function generates and returns the coefficient matrix corresponding to the HOC discretization
% of the Poisson equation. The resulting matrix is prod(grid) by prod(grid). Once returned it needs
% to be divided by 30dx^2 where dx is the grid spacing which must be uniform in all 3 directions.
%
% Define the laplacian using the 27-point 6th order compact stencil of Spotz & Carey '95
%   The basic 27 point stencil is:          |   The HOC 27 point stencil is:
%   -3 -4 -3      -4 -4 -4      -3 -4 -3    |   1   3   1     3   14    3   1   3   1
%   -4 -4 -4      -4 96 -4      -4 -4 -4    |   3  14   3     14 -128  14   3  14   3
%   -3 -4 -3      -4 -4 -4      -3 -4 -3    |   1   3   1     3   14    3   1   3   1
%
%>> grid        Size of coefficient matrix to generate, 3x1                        double
%<< result      Coefficient matrix (not divided by 30dx^2)                         sparse

    %--- Prepare to create the matrix ---%
    %       Create coordinate matrices and linear indexes and store some useful constants.
    Nx         = grid(1);
    Ny         = grid(2);
    Nz         = grid(3);
    Nxy        = Nx*Ny;
    numCells = prod(grid);

    [X Y Z]    = ndgrid(1:Nx, 1:Ny, 1:Nz);
    indices    = reshape(1:numCells, grid);

    %--- Specify finite difference coefficients ---%
    CEN = -24; % center cell
    FC  = 2;   % touching face of center cell
    EC  = 1;    % touching edge of center cell
    CC  = 0;    % touching corner of center cell

    %--- Begin creating the matrix ---%
    %       Each row in the laplacian matrix will nominally have 27 nonzero entries (less around the
    %       edges). To add to a matrix M of size i by j, coefficientMatrixSubset takes a set of i 
    %       indices, a set of j indices, a constant to insert, an amount to shift by, and a set of 
    %       locations to not act on (due to shifting ending up outside the grid) in that order.
    % Fourth argument: x_shift + y_shift + z_shift

    result =     coefficientMatrixSubset(indices, indices, CEN,  0+0*Nx+0*Nxy, (X <= Nx ), numCells);

    result = result + coefficientMatrixSubset(indices, indices, FC,  1+0*Nx+0*Nxy, (X < Nx), numCells);
    result = result + coefficientMatrixSubset(indices, indices, FC, -1+0*Nx+0*Nxy, (X > 1 ), numCells);
    result = result + coefficientMatrixSubset(indices, indices, FC,  0+1*Nx+0*Nxy, (Y < Ny), numCells);
    result = result + coefficientMatrixSubset(indices, indices, FC,  0-1*Nx+0*Nxy, (Y > 1 ), numCells);
    result = result + coefficientMatrixSubset(indices, indices, FC,  0+0*Nx+1*Nxy, (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, FC,  0+0*Nx-1*Nxy, (Z > 1 ), numCells);

    result = result + coefficientMatrixSubset(indices, indices, EC,  1+0*Nx-1*Nxy, (X < Nx) & (Z > 1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC, -1+0*Nx-1*Nxy, (X > 1 ) & (Z > 1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC,  0+1*Nx-1*Nxy, (Y < Ny) & (Z > 1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC,  0-1*Nx-1*Nxy, (Y > 1 ) & (Z > 1), numCells);

    result = result + coefficientMatrixSubset(indices, indices, EC,  1+1*Nx+0*Nxy, (X < Nx) & (Y < Ny), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC,  1-1*Nx+0*Nxy, (X < Nx) & (Y >  1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC, -1+1*Nx+0*Nxy, (X > 1 ) & (Y < Ny), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC, -1-1*Nx+0*Nxy, (X > 1 ) & (Y >  1), numCells);

    result = result + coefficientMatrixSubset(indices, indices, EC,  1+0*Nx+1*Nxy, (X < Nx) & (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC, -1+0*Nx+1*Nxy, (X > 1 ) & (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC,  0+1*Nx+1*Nxy, (Y < Ny) & (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, EC,  0-1*Nx+1*Nxy, (Y > 1 ) & (Z < Nz), numCells);

    result = result + coefficientMatrixSubset(indices, indices, CC,  1+1*Nx-1*Nxy, ...
                                                    (X < Nx) & (Y < Ny) & (Z > 1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, CC,  1-1*Nx-1*Nxy, ...
                                                    (X < Nx) & (Y > 1 ) & (Z > 1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, CC, -1+1*Nx-1*Nxy, ...
                                                    (X > 1 ) & (Y < Ny) & (Z > 1), numCells);
    result = result + coefficientMatrixSubset(indices, indices, CC, -1-1*Nx-1*Nxy, ...
                                                    (X > 1 ) & (Y > 1 ) & (Z > 1), numCells);

    result = result + coefficientMatrixSubset(indices, indices, CC,  1+1*Nx+1*Nxy, ...
                                                    (X < Nx) & (Y < Ny) & (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, CC,  1-1*Nx+1*Nxy, ...
                                                    (X < Nx) & (Y > 1 ) & (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, CC, -1+1*Nx+1*Nxy, ...
                                                    (X > 1 ) & (Y < Ny) & (Z < Nz), numCells);
    result = result + coefficientMatrixSubset(indices, indices, CC, -1-1*Nx+1*Nxy, ...
                                                    (X > 1 ) & (Y > 1 ) & (Z < Nz), numCells);
end

