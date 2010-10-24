function magnetFlux(run, mass, mom, mag, X, magneticIndices)
%   This is the routine responsible for fluxing the magnetic field array according to a TVD MUSCL
%   scheme that utilizes constrained transport to maintain the divergencelessness of the magnetic
%   field over time if it is divergenceless as an initial condition.
%
%>< run                 Run manager object.                                     ImogenManager
%>< mass                Mass density array.                                     FluidArray
%>< mom                 Momentum density array.                                 FluidArray(3)
%>< ener                Energy density array.                                   FluidArray
%>< mag                 Magnetic field array.                                   MagnetArray(3)
%>> X                   Vector index along which to flux.                       int
%>> magneticIndices     Index values for the magnetic components to flux        int(2)
 
    for n=1:2
        I = magneticIndices(n); % Component of the b-field to flux.

        % Prepare arrays for magnetic flux step at cell corners (2nd order grid-aligned velocity)
        mag(I).velGrid(X).array = (mom(X).shift(I,-1) + mom(X).array) ...
                                    ./ (mass.shift(I,-1) + mass.array);

        mag(I).velGrid(X).array = 0.25*(mag(I).velGrid(X).shift(X,-1) ...
                                        + mag(I).velGrid(X).shift(X,1) ...
                                        + 2*mag(I).velGrid(X).array);
        
        % Flux Advection step
        relaxingMagnet(run, mag, mag(I).velGrid(X), X, I);
        mag(I).cleanup();
    end
    
    mag(I).updateCellCentered();
    mag(X).updateCellCentered();
end
