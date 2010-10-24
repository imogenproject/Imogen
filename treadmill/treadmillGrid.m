function treadmillGrid(run, mass, mom, ener, mag)
% This function will remap the entire grid in an attempt at keeping the cell with the fastest
% velocity in the specified direction stationary. It was originally designed for use with the
% corrugation instability problem where the shockwave needs to be kept stationary on the grid.
%
%>< run		run variable manager objec										ImogenManager
%>< mass	mass density array												FluidArray
%>< mom 	momentum density array											FluidArray(3)
%>< ener 	energy density array											FluidArray
%>< mag		magnetic field density array									MagnetArray(3)

    
    % Exit function if inactive or if the index is too small to ensure a stable compression reading.
    if (~run.treadmill.ACTIVE || (run.time.iteration < 10)); return; end
 
	X = double(run.treadmill.DIRECTION);
	
    %--- Find the consitent maximum change in the mass density ---%
    [maxVals, maxIndeces] = max(mass.calculate2PtDerivative(X,run.DGRID{X}), [], X);
	if ~isa(maxIndeces,'double'),		maxIndeces = gather(maxIndeces); end %r2009b: iscodistributed
    newIndex = round(mean(mean(maxIndeces)));

    %--- Number of cells to treadmill ---%
	if (run.treadmill.last < 0)
        run.treadmill.last = newIndex;
        return;
	end
	    
    %--- Act on primary arrays ---%
	dCells      = newIndex - run.treadmill.last;
	mass.array  = mass.shift(X, dCells);
	ener.array  = ener.shift(X, dCells);
	for i=1:3
		mom(i).array  = mom(i).shift(X, dCells);
		mag(i).array  = mag(i).shift(X, dCells);
	end
	run.treadmill.appendHistory(dCells);
end
