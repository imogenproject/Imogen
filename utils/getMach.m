function result = getMach(mass, mom, ener, mag, GAMMA)
% Calculates the mach number array from the input values. This routine calculates gas pressure and 
% sound speed internally instead of using the getPressure function to reduce computation overhead.
%
%>< mass		mass density object													FluidArray
%>< mom			momentum density object												FluidArray(3)
%>< ener		energy density object												FluidArray
%>< mag			magnetic field density object										MagnetArray(3)
%>> GAMMA		polytropic constant													double

	velSquared = getVelocitySquared(mass,mom);
	magSquared = getMagneticSquared(mag);
	
	if isa(mass.array,'GPUdouble')
	result = (GAMMA - 1.0)*(ener.array - 0.5*mass.array .* velSquared - 0.5*magSquared);
	cudaArrayAtomic(result, 0, ENUM.CUATOMIC_SETMIN);
	result = sqrt( velSquared ./ abs((GAMMA*result + 2.0*magSquared) ./ mass.array) );
	else
	result = sqrt( velSquared ./ abs((GAMMA*max( (GAMMA - 1.0)*(ener.array ...
					- 0.5*mass.array .* velSquared - 0.5*magSquared), 0.0 ) ...
					+ 2.0*magSquared) ./ mass.array) );
	end
end
