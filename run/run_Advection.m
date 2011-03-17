% Run Advection test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run             = AdvectionInitializer([64 4 4]);
run.iterMax     = 250;
run.info        = 'Advection test.';
run.notes       = 'Simple advection test in the x-direction.';

run.image.interval = 5;
run.image.mass = true;
%run.ppSave.dim2 = 20;

run.useGPU      = true;
run.gpuDeviceNumber = 2;

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
