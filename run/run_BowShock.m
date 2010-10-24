% Run a bow shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = BowShockInitializer();
run.iterMax         = 8000;
run.image.interval  = 5;
run.image.mass      = true;
run.info            = 'Bow shock test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();