function result = magneticDivergence(run, mag)
%  This routine calculates the cell averaged divergance (del o B) for the input magnet array for
%  both 1D, 2D and 3D cases.
%
%>< run     data manager object                                             ImogenManager       H
%>< mag     magnetic field array (vector)                                   MagneticArray(3)    H
%<< result  numerical divergence of magnetic field array                    double(GRID)
    
    result = zeros(run.gridSize);
    for i=1:3
        result = result + mag(i).calculate2PtDerivative(i,run.dGrid{i});
    end

end
