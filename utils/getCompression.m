function result = getCompression(run, array, divideArray)
%   This routine returns the compession of the array, typically a velocity. The compression is
%   defined as the divergence of the vector field. If a divideArray is specified, it is divided out
%   of each component of the array (e.g. array = momentum, divideArray = mass to get velocity).
%
%>< array			array to find the compression of						ImogenArray(3)
%>> dGrid			grid spacing vector										double [dx dy dz]
%>< divideArray		array to divide out of the primary array				ImogenArray
%<< result			resulting compression of array							double [Nx Ny Nz]


    %--- Initialization ---%
	if iscodistributed(divideArray.array)
		result = codistributed.zeros(run.gridSize, run.parallel.distribution);
	else
		result = zeros(run.gridSize);
	end
	
	%--- Calculate compression ---%
    for i=1:3
		vec       = array(i).dataClone();
		vec.array = vec.array ./ divideArray.array;
        result    = result + vec.calculate5PtDerivative(i,run.dGrid{i});
    end
    
end
