function speed = speedFromMach(mach, GAMMA, mass, ener, mag)
% This routine determines the speed for the provided mach for the given primary variable array
% values.
%
%>> mach 	mach number for the desired speed                               double      #
%>> GAMMA	index for the polytropic equation of state                      double      #
%>> mass	mass density array                                              double      [nx ny nz]
%>> ener	energy density array                                            double      [nx ny nz]
%>> mag 	magic field array												double      [3 nx ny nz]
%<< speed     resulting speed array                                         double      [nx ny nz]


    %--- Calculate cell centered magic field ---%
	if ( nargin < 5 || isempty(mag) )
		magCell = 0;
	else
		magCell = zeros(size(mag)); %Interpolate the face-centered b field to cell centers
		try
			magCell(1,:,:,:) = 0.5*( mag(1,:,:,:) + circ_shift( mag(1,:,:,:), 2, 1) );
			magCell(2,:,:,:) = 0.5*( mag(2,:,:,:) + circ_shift( mag(2,:,:,:), 3, 1) );
			magCell(3,:,:,:) = 0.5*( mag(3,:,:,:) + circ_shift( mag(3,:,:,:), 4, 1) );
		catch ME
		end
	end
	
    GAMMAm1 = GAMMA - 1;

    %--- Determine speed ---%
    speed = mach * sqrt( (GAMMA * ener ...
                          + ( 2.0 ./ (GAMMAm1 * mass) - GAMMA * GAMMAm1 / 2.0 ) ...
                          .* squeeze(sum(magCell.*magCell,1)) ...
                          ) / ( 1/GAMMAm1 + GAMMA * (mach*mach) / 2 ) );
                      
end
