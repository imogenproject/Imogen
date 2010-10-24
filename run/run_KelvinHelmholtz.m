% Run a test of the Kelvin-Helmholtz instability test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = KelvinHelmholtzInitializer();
run.iterMax         = 1000;
run.direction       = KelvinHelmholtzInitializer.X;
run.image.interval	= 2;
run.image.mass		= true;
run.image.mach		= true;
run.info            = 'Kelvin-Helmholtz instability test.';
run.notes           = '';


%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();