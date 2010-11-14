% Run a bow shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = BowShockInitializer([512 192 1]);
run.iterMax         = 500;

run.image.interval  = 5;
run.image.mass      = true;
run.image.speed     = true;

run.ppSave.dim2     = 25;

run.info            = 'Bow shock test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
