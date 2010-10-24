function result = getVorticity(run, imoArray, divideArray)
%   This routine finds the vorticity of the input array (usually a velocity), which is defined by
%   the curl of the vector field. If a divideArray is supplied, the will divide each of the array
%   components by the divide array (e.g. array = momentum, divideArray = mass to get velocity).
%
%>< run             Imogen run manager object.                              ImogenManager
%>< array           Array to find the vorticity of, usually a velocity.     ImogenArray(3)
%>< divideArray     Array to divide each component of the array by.         ImogenArray
%<< return          Resulting vorticity result.                             double(3,Nx,(Ny),(Nz))

    %--- Initialization ---%
    N = [3 run.gridSize];
    if ~isa(divideArray.array,'double') %r2009b: iscodistributed
        parallel = ParallelManager.getInstance();
        result = zeros(N, parallel.distribution);
    else    result = zeros(N);
    end
    
    aX = imoArray(1).dataClone();
    aY = imoArray(2).dataClone();
    aZ = imoArray(3).dataClone();

    result(1,:,:,:) = aZ.calculate5PtDerivative(2,run.dGrid{2}) - aY.calculate5PtDerivative(3,run.dGrid{3});
    result(2,:,:,:) = aX.calculate5PtDerivative(3,run.dGrid{3}) - aZ.calculate5PtDerivative(1,run.dGrid{1});
    result(3,:,:,:) = aY.calculate5PtDerivative(1,run.dGrid{1}) - aX.calculate5PtDerivative(2,run.dGrid{2});

end
