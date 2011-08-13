function [result, aux] = pressure(mode, run, mass, momvel, ener, mag)
% Solve for the pressure based on the type of pressure desired as specified by the mode
% variable. Additional corrective checks at the beginning of the routine ensure that the dot
% product arrays squeeze down to the correct dimensional format for element-wise arithmetic
% operations with the scalar arrays.
%
%>> mode    Type of pressure to find (gas, total, fluid, magnetic, sound).  str
%>< run     Imogen run manager object.                                      ImogenManager
%>< mass    Mass density array.                                             FluidArray
%?  momvel  Either the momentum or the velocity squared array.              FluidArray(3), double
%>< ener    Energy density array.                                           FluidArray
%>< mag     Magnetic field array.                                           MagnetArray(3)
%<< result  Resulting pressure based on mode.                               double        [Nx Ny Nz]
%<< aux     Resulting sound speed array.                                    double        [Nx Ny Nz]

    GAMMA = run.GAMMA;
    aux   = [];

    % This runs the below code on the GPU, basically, but only reads each variable from memory once and generates a single write per cell.
    if run.useGPU
        switch( mode )
            case ENUM.PRESSURE_TOTAL_AND_SOUND
                [result aux] = cudaMHDKernels(5, mass.array, ener.array, momvel(1).array, momvel(2).array, momvel(3).array, mag(1).cellMag.array, mag(2).cellMag.array, mag(3).cellMag.array, GAMMA);                
            case ENUM.PRESSURE_SOUND_SPEED
                result = cudaMHDKernels(1, mass.array, ener.array, momvel(1).array, momvel(2).array, momvel(3).array, mag(1).cellMag.array, mag(2).cellMag.array, mag(3).cellMag.array, GAMMA);
            case ENUM.PRESSURE_GAS
                result = cudaMHDKernels(2, mass.array, ener.array, momvel(1).array, momvel(2).array, momvel(3).array, mag(1).cellMag.array, mag(2).cellMag.array, mag(3).cellMag.array, GAMMA);
            case ENUM.PRESSURE_TOTAL
                result = cudaMHDKernels(3, mass.array, ener.array, momvel(1).array, momvel(2).array, momvel(3).array, mag(1).cellMag.array, mag(2).cellMag.array, mag(3).cellMag.array, GAMMA);
            case ENUM.PRESSURE_MAGNETIC
                result = cudaMHDKernels(4, mag(1).cellMag.array, mag(2).cellMag.array, mag(3).cellMag.array);
        end
        return;
    end

    % Prepare the velocity squared array
    if isa(momvel,'double')
        velSquared = momvel;
    else
        velSquared = (momvel(1).array .* momvel(1).array + momvel(2).array .* momvel(2).array ...
                    + momvel(3).array .* momvel(3).array) ./ (mass.array.*mass.array);
    end
    
    % Prepare the cell-centered magnet squared array
    magSquared = 0;
    for i=1:3
        if ~mag(i).isZero
            magSquared = magSquared + mag(i).cellMag.array .* mag(i).cellMag.array;
        end
    end

    %Calculate the fluid pressure
    result = (GAMMA - 1.0)*(ener.array - 0.5*mass.array .* velSquared);
        
    %-----------------------------------------------------------------------------------------------
    % Calculate the pressure based on mode type
    %------------------------------------------
    switch (mode)
        case ENUM.PRESSURE_TOTAL_AND_SOUND
            aux = sqrt(abs( (GAMMA*(result - (GAMMA - 1.0)*0.5*magSquared) + 2.0*magSquared) ./ mass.array ));
            result = result + 0.5*(2.0 - GAMMA)*magSquared;
        
        case ENUM.PRESSURE_SOUND_SPEED
            result = result - (GAMMA - 1.0)*0.5*magSquared;
            result = sqrt(abs( (GAMMA*result + 2.0*magSquared) ./ mass.array ));

        case ENUM.PRESSURE_GAS
            result = result - (GAMMA - 1.0)*0.5*magSquared;

        case ENUM.PRESSURE_TOTAL
            result = result + 0.5*(2.0 - GAMMA)*magSquared;

        case ENUM.PRESSURE_MAGNETIC
            if (magSquared == 0)
                result = zeros(mass.gridSize);
            else
                result = 0.5*magSquared;
            end
        end 

if run.useGPU        
    cudaArrayAtomic(result, 0.0, ENUM.CUATOMIC_SETMIN);
else
    result(result < 0) = 0;
end

end
