% Run Advection test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run             = AdvectionInitializer([1024 4 4]);
run.iterMax     = 1000;
run.info        = 'Advection test.';
run.notes       = 'Simple advection test in the x-direction.';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();