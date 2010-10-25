function timeres = computeTimeFromHistory(hist, findwhat, value)
% Given a timestep history from imogen (Nx1 variable),
% looks at findwhat:
% 'frame': Gives the time elapsed at step #value
% 'time' : Gives the frame after at least time 'value' has passed

maxlen = max(size(hist));

if findwhat == 'frame'
	if value > maxlen
		disp 'Not that many steps; displaying total time';
		value = maxlen;
	end;
	
	tpass = sum(hist(1:value));
	fprintf('Time elapsed after frame %i: %d\n', value, tpass);
	timeres = tpass; return;

end;

if findwhat == 'time'
	tpass = 0;	

	for i = 1:maxlen;
		tpass = tpass + hist(i);
		if tpass > value
			fprintf('Iteration %i has passed %d time > given %d\n', i, tpass, value);
			break;
		end;
	end;

	timeres = tpass; return;

end;

end
