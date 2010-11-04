function [i,j,k] = indexToijk(index, grid, natural)

%=== INDEXTOIJK ====================================================================================
% Converts a 1D monotonic index into a 3D grid position (i.e. [i j k]) that is appropriate for
% reshaping the after the Poisson solver. This routine is vectorized and will accept an array of
% index values.
%===================================================================================================
%
%-CALLS---------------------------------------------------------------------------------------------
%   [i,j,k] = indexToijk(index, grid, natural)
%-INPUTS--------------------------------------------------------------------------------------------
% index     the index value (or array of index values)                          int     #
% grid      the 3D grid size                                                    int     [Nx Ny Nz]
% natural   toggle between standard and natural ordering of the sparse system   bool    t/f
%               NATURAL:  111, 211, ..., (Nx)11, 121, 221, ..., (Nx)21, ...
%               STANDARD: 111, 112, ..., 11(Nz), 121, 122, ..., 12(Nz), ...  
%-OUTPUTS-------------------------------------------------------------------------------------------
% i       first element cell position                                           int     #
% j       second element cell position                                          int     #
% k       third element cell position                                           int     #
%---------------------------------------------------------------------------------------------------

    if (natural) % NATURAL Ordering
        i = mod(index-1,grid(1))+1;
        j = mod(ceil(index / grid(3))-1,grid(1))+1;
        k = ceil(index / prod(grid(1:2)));
    else         % ORIGINAL Ordering
        k = mod(index-1,grid(3))+1;
        j = mod(ceil(index / grid(3))-1,grid(2))+1;
        i = ceil(index / prod(grid(2:3)));
    end
    
    %--- Debug output ---%
%     for n=1:length(i)
%         indexRes = (k(n)-1)*prod(grid(1:2))+(j(n)-1)*grid(1)+i(n);
%         disp(sprintf('(%g, %g, %g) --> [%g = %g]\n',i(n),j(n),k(n),index(n),indexRes)); 
%     end
end