%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize bow shock ---%
grid = [64 64 64];
run                 = BowShockInitializer(grid);
run.iterMax         = 500;
%run.bcMode.z	    = 'circ';

run.bcMode.x = 'circ';
run.bcMode.y = 'circ';
run.bcMode.z = 'circ';

run.useGPU = true;
run.gpuDeviceNumber = 2;

%--- Adjustable simulation parameters ---%

run.ballXRadius = 1;
run.ballCells = [31.5 31.5 31.5]/2;
run.ballCenter =  [32 32 32];

run.mode.magnet = false;

%run.magX = 1;
%run.magY = 0;

run.bgVx      = 4;
run.bgRho     = .125;
run.ballVr    = 1;
run.ballRho   = 1;

%--- Adjustable output parameters ---%
run.image.interval  = 1;
run.image.mass      = true;
run.image.speed     = true;
%run.image.magX      = true;
%run.image.magY = true;

run.activeSlices.xy  = true;
run.activeSlices.yz  = true;
run.activeSlices.xyz = true;

run.ppSave.dim2     = 5;
run.ppSave.dim3     = 20;


run.info            = 'Bow shock test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    IC.mass = mass;
    IC.mom = mom;
    IC.ener = ener;
    IC.magnet = magnet;
    IC.statics = statics;
    IC.ini = ini;
    icfile = [tempname '.mat'];

    save(icfile, 'IC');
    clear IC mass mom ener magnet statics ini run;
    imogen(icfile);
end

enderRun();
