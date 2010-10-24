function velGrid = getVelocityGrid(mass, mom, X, I) 
% Finds the grid-aligned velocity from the cell centered momentum for use in the staggered grid
% magnetic field fluxing routine. This is found using a second order interpolation method as
% recommended by Trac & Hu.
%

%><	mass		mass density object											FluidArray
%>< mom         momentum density object										FluidArray(3)
%>> X			vector index for fluxing direction							int
%>>	I			component of the magnetic field to find						int
%<< velGrid		resulting X-component grid-centered velocity				FluidArray
	
    velGrid	= mom(X).dataClone();
	
	%--- Interpolate along I direction ---%
    velGrid.array	= (velGrid.shift(I,-1) + velGrid.array) ./ (mass.shift(I,-1) + mass.array);
    
	%--- Update edge conditions if they exist ---%
	if isstruct( mom(X).edgeStore.min )
		fields = fieldnames( mom(X).edgeStore.min );
		for i=1:length(fields)
			velGrid.edgeStore.min.(fields{i}) = ...
								mom(X).edgeStore.min.(fields{i}) ./ mass.edgeStore.min.(fields{i});
			velGrid.edgeStore.max.(fields{i}) = ...
								mom(X).edgeStore.max.(fields{i}) ./ mass.edgeStore.max.(fields{i});
		end
	end
	
	%--- Interpolate along X direction ---%
    velGrid.array = 0.25*(velGrid.shift(X,-1) + velGrid.shift(X,1) + 2*velGrid.array);
    
end