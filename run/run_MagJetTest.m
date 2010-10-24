% Run a magnetic jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer();
run.mode.magnet     = true;
run.clf             = 0.35;
run.image.interval	= 3;
run.image.mass		= true;
run.image.mach		= true;
run.info            = 'Magnetic jet test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();


%---------------------------------------------------------------------------------------------------
% Run tests
%----------
if (true) %Parallel Magnetic field 
	ics.jetMags	 = [magAmp 0 0];
	ics.backMags = [magAmp 0 0];
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

if (false) %Non-Zero By test
	ics.jetMags	 = [0 magAmp 0];
	ics.backMags = [0 magAmp 0];
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

if (false) %Non-Zero Bx By test
	run.jetMags	 = [magAmp magAmp 0];
	run.backMags = [magAmp magAmp 0];
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();