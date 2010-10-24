function result = getVelocitySquared(mass, mom)
% Determines the magnitude of the velocity squared from the momentum and mass array objects.
%
%>< mass	mass density object														FluidArray
%>< mom		momentum density object													FluidArray(3)


	result = (mom(1).array .* mom(1).array + mom(2).array .* mom(2).array ...
             + mom(3).array .* mom(3).array) ./ (mass.array.*mass.array);

end
