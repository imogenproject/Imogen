% Run a gravitational jet test

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer();
run.mode.gravity    = true;
run.direction       = JetInitializer.X;
run.image.interval	= 3;
run.image.mass		= true;
run.image.mach		= true;
run.info            = 'Gravity jet test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
