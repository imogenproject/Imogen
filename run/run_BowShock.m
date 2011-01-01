%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize bow shock ---%
grid = [128 192 192];
run                 = BowShockInitializer(grid);
run.iterMax         = 2000;
run.bcMode.z	    = 'flip';

%--- Adjustable simulation parameters ---%

run.ballXRadius = 1;
run.ballCells = [24 24 24];
run.ballCenter =  [64 96.5 96.5];
run.bgVx      = 5;
run.bgRho     = .125;
run.ballVr    = 2;
run.ballRho   = 1;

%--- Adjustable output parameters ---%
run.image.interval  = 10;
run.image.mass      = true;
run.image.speed     = true;

run.activeSlices.xyz = true;

run.ppSave.dim2     = 25;
run.ppSave.dim3     = 5;


run.info            = 'Bow shock test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
