%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer([3072 1024 1]);
run.iterMax         = 10000;
run.offset          = [40 512 1];
run.bcMode.x        = 'circ';
run.bcMode.y        = 'circ';
run.direction       = JetInitializer.X;
run.flip            = false;

run.image.interval  = 15;
run.image.mass      = true;
run.image.speed     = true;

run.info            = 'Fluid jet test.';
run.notes           = '';

run.injectorSize = 30;
run.jetMass = .65;
run.jetMach = 3;

run.useGPU = true;
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
