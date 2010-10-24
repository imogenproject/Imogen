function [R Theta] = trackArmPhase(pdata, r0)
% >> pdata: unwrapped disk information
% << R: Radius of arm
% << Theta: phase

% Find how many radial elements there are to examine
nRelements = size(pdata, 1) - r0(1) + 1;

% Create the result holders
R = [];
Theta = [];

% Set the initial 'phase' and the width to search for maxima
presphase = r0(2); 
pitch = 40;

for rhat = 1:nRelements
	% Create a selection that wraps around the y direction
	thetasel = mod([ (presphase-pitch):(presphase+pitch)]-1, size(pdata,2))+1;

	% Grab it and compute backwards finite difference
	checkarea = pdata(r0(1)+rhat- 1, thetasel);
	diff = abs(checkarea - circ_shift(checkarea,2,-1));

	diff(1) = 0; diff(end) = 0; % Kill the boundary cells since they don't represent
	% actual jumps

	% Set R
	R(rhat) = r0(1)+rhat-1;

	% Identify the location of the maximum jump
	md = 0; mdi = 1;
	for q = 1:2*pitch;
		if diff(q) > md; md = diff(q); mdi = q; end
	end

% debug stuff
%	if rhat < 15
%		figure(3); plot(diff); hold on; scatter(mdi, diff(mdi)); xxx = input('continue: ');
%		fprintf('%3.3g ', mdi - pitch);
%	end

	% Set the phase value
	Theta(rhat) = presphase - pitch + mdi ;
	presphase = Theta(rhat);
end

% Wrap around
Theta = mod(Theta,4096);

end
