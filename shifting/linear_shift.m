function result = linear_shift(array, DIRECT, n, obj)
% Shifts an array so that the new cells are constants based on the previous edge values.
%
%>> array		array to be shifted											double((3),Nx,Ny,Nz)
%>> DIRECT		direction/dimension along which to shift					int
%>> n			number of cells to shift									int
%>< obj			object owning the array										ImogenArray...
%<< result		shifted array												double((3),Nx,Ny,Nz)

	%--- Initialization ---%
	N = size(array);
    NDim  = N(DIRECT);
    NSize = length(N);
    index = cell(1, NSize);
    for i=1:NSize,  index{i} = 1:N(i); end
    
    dir = -sign(n); %Flip sign to correct inverse nature of shifting routines
    n = abs(n);

	%--- SHIFT Array ---%
	if (dir>0)
		index{dim}((n+1):NDim) = index{dim}(1:(NDim-n));
		index{dim}(1:n) = 1;
		result = array(index{:});
		npIndex = index; npIndex{dim} = n;
		nppIndex = index; nppIndex{dim} = min(n+1,NDim);
		delta = array(nppIndex{:}) - array(npIndex{:});
		for i=n:-1:1
			npIndex{dim} = i;
			nppIndex{dim} = i+1;
			result(npIndex{:}) = result(nppIndex{:}) - delta;
		end
		return;
	else
		index{dim}(1:(NDim-n)) = index{dim}((n+1):NDim);
		index{dim}((NDim-n+1):NDim) = NDim;
		result = array(index{:});
		nmIndex = index; nmIndex{dim} = NDim - n;
		nmmIndex = nmIndex; nmmIndex{dim} = nmmIndex{dim} - 1;
		delta = array(nmIndex{:}) - array(nmmIndex{:});
		for i=(NDim-n+1):NDim
		   nmIndex{dim} = i;
		   nmmIndex{dim} = i-1;
		   result(nmIndex{:}) = result(nmmIndex{:}) + delta;
		end
		return;
	end 
end