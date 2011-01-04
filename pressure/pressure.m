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
%<< result  Resulting pressure based on mode.                               double	[Nx Ny Nz]
%<< aux     Resulting sound speed array.                                    double	[Nx Ny Nz]

    GAMMA = run.GAMMA;
	aux   = [];

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
			aux = sqrt(abs( (GAMMA*(result - (GAMMA - 1.0)*0.5*magSquared) ...
					+ 2.0*magSquared) ./ mass.array ));
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
    result = double(result);
    result(result < 0) = 0.0;
    result = GPUdouble(result);
else
    result(result < 0) = 0;
end
end
