%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize bow shock ---%
grid = [512 256 1];
run                 = BowShockInitializer(grid);
run.iterMax         = 1500;

%--- Adjustable simulation parameters ---%

run.ballRadii = [24 24 1];
%run.ballCenter =  [centerx centery centerz];
run.bgVx      = 4;
%run.bgRho     = 1;
%run.ballVr    = 1;
%run.ballRho   = 1;

%--- Adjustable output parameters ---%
run.image.interval  = 5;
run.image.mass      = true;
run.image.speed     = true;

run.ppSave.dim2     = 10;

run.info            = 'Bow shock test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
