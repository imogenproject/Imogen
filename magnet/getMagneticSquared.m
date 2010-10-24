function result = getMagneticSquared(mag)
% Calculates the magnitude of the magnetic field squared. Useful in pressure and energy computations.
%
%>< mag			magnetic field density object									MagnetArray
%<< result		magnitude of the magnetic field squared array					double [nx ny nz]

	result = mag(1).array .* mag(1).array;
	result = mag(2).array .* mag(2).array + result;
	result = mag(3).array .* mag(3).array + result;

end