%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize bow shock ---%
grid = [768 768 1];
run                 = BowShockInitializer(grid);
run.iterMax         = 1000;
%run.bcMode.z	    = 'circ';

run.bcMode.x = 'circ';
run.bcMode.y = 'circ';

run.useGPU = true;
run.gpuDeviceNumber = 2;

%--- Adjustable simulation parameters ---%

run.ballXRadius = 1;
run.ballCells = [31.5 31.5 31.5];
run.ballCenter =  [128 384 1];

run.mode.magnet = true;


run.magX = 1;
run.magY = 0;

run.bgVx      = 3;
run.bgRho     = .25;
run.ballVr    = 0;
run.ballRho   = 1;

%--- Adjustable output parameters ---%
run.image.interval  = 8;
run.image.mass      = true;
run.image.speed     = true;
run.image.magX      = true;
run.image.magY = true;

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
