% Run Sod shock tube test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run             = SodShockTubeInitializer([512 512 4]);
run.direction   = SodShockTubeInitializer.X;
run.shockAngle  = 45;
run.iterMax     = 1000;
run.timeMax     = 0.15;

run.alias       = '';
run.info        = 'Sod shock tube test.';
run.notes       = '512 resolution non-axis-aligned Sod shock tube test in the XY direction (theta = 45 degrees).';

run.ppSave.dim2 = 1;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();