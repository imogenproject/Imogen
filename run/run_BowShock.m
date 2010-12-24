%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize bow shock ---%
grid = [256 128 128];
run                 = BowShockInitializer(grid);
run.iterMax         = 1000;
%run.bcMode.y	    = 'mirror';

%--- Adjustable simulation parameters ---%

run.ballRadii = [24 24 1];
run.ballCenter =  [96 64.5 64.5];
run.bgVx      = 4;
%run.bgRho     = 1;
run.ballVr    = 1;
%run.ballRho   = 1;

%--- Adjustable output parameters ---%
run.image.interval  = 5;
run.image.mass      = true;
run.image.speed     = true;

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
