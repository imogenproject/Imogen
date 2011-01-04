%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer([512 256 1]);
run.iterMax         = 25;
run.offset          = [20 128 1];
run.bcMode.x        = 'circ';
run.bcMode.y        = 'circ';
run.direction       = JetInitializer.X;
run.flip            = false;

run.image.interval	= 10;
run.image.mass		= true;
run.image.speed = true;

run.info            = 'Fluid jet test.';
run.notes           = '';

run.jetMass = 3;
run.jetMach = 5;

run.useGPU = true;
run.gpuDeviceNumber = 2;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
