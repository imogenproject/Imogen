function index = ijkToIndex(ijk, grid, natural)

%=== IJKTOINDEX ====================================================================================
% Converts a 3D grid position (i.e. [i j k]) into a linear index value that is appropriate for
% reshaping the array for the Poisson solver. This function is vectorized and will accept an array
% of ijk values, where each row is an ijk triplet.
%===================================================================================================
%
%-CALLS---------------------------------------------------------------------------------------------
%   index = ijkToIndex(ijk, grid, natural)
%-INPUTS--------------------------------------------------------------------------------------------
% ijk       array describing the cell position                                  int     [i j k]
% grid      the 3D grid size                                                    int     [Nx Ny Nz]
% natural   toggle between standard and natural ordering of the sparse system   bool    t/f
%               NATURAL:  111, 211, ..., (Nx)11, 121, 221, ..., (Nx)21, ...
%               STANDARD: 111, 112, ..., 11(Nz), 121, 122, ..., 12(Nz), ...  
%-OUTPUTS-------------------------------------------------------------------------------------------
% index     the resulting index value that corresponds to the input ijk         int     #
%---------------------------------------------------------------------------------------------------

    if (natural) % NATURAL  Ordering
        index = (ijk(:,3)-1)*prod(grid(1:2)) + (ijk(:,2)-1)*grid(1) + ijk(:,1);
    else         % STANDARD Ordering
        index = (ijk(:,1)-1)*prod(grid(2:3)) + (ijk(:,2)-1)*grid(3) + ijk(:,3);
    end
end