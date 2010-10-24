function result = flip_shift(array, DIRECT, n, obj)
% Shifts an array so that the new cells added are identical to the pre-shifted edge cells and then
% flips the magnitude of any vector components that are aligned along the specified direction. This
% is used for simulating mirroring where the edge cells are repeated, in a mirrored fashion, across
% the edges.
%
%>> array		array to be shifted											double((3),Nx,Ny,Nz)
%>> DIRECT		direction/dimension along which to shift					int
%>> n			number of cells to shift									int
%>< obj			object owning the array										ImogenArray...
%<< result		shifted array												double((3),Nx,Ny,Nz)

	%--- Initialization ---%
	N		= size(array);
    NDim	= N(DIRECT);
    NSize	= length(N);
    index	= cell(1, NSize);
    for i=1:NSize,  index{i} = 1:N(i); end
    
    dir = -sign(n); %Flip sign to correct inverse nature of shifting routines
    n	= abs(n);

	%--- SHIFT Array ---%
    if (dir>0)
		index{DIRECT}((n+1):NDim)       = index{DIRECT}(1:(NDim-n));
		index{DIRECT}(1:n)              = 1;
        flipIndices                     = (1:n);
    else
		index{DIRECT}(1:(NDim-n))		= index{DIRECT}((n+1):NDim);
		index{DIRECT}((NDim-n+1):NDim)	= NDim;
        flipIndices                     = (NDim-n+1):NDim;
    end
	
    result = array(index{:});
    
    %--- Flip sign on appropriate vector indeces ---%
    if (obj.component == DIRECT)
        index{DIRECT}    = flipIndices;
        result(index{:}) = -result(index{:});
    end
    
end