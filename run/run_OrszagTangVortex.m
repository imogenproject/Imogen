% Run 2D Orszag-Tang vortex test problem.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = OrszagTangVortexInitializer([256 256 1]);
run.info            = 'Orszag-Tang vortex test.';
run.notes           = '';
run.profile         = true;
run.image.interval	= 3;
run.image.mass		= true;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
