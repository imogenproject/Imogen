% Run Sod shock tube test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run             = SodShockTubeInitializer([512 16 1]);
run.direction   = SodShockTubeInitializer.X;
run.shockAngle  = 0;
run.iterMax     = 1000;
run.timeMax     = 0.15;

run.alias       = '';
run.info        = 'Sod shock tube test.';
run.notes       = 'Simple axis aligned shock tube test';

run.useGPU = false;
run.gpuDeviceNumber = 2;

run.ppSave.dim2 = 1;

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
