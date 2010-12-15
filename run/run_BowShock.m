% Run a bow shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
grid = [128 128 1];
run                 = BowShockInitializer(grid);
run.iterMax         = 100;

run.ballRadii = [32 32 inf];

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
