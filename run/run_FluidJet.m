%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer([2048 1024 1]);
run.iterMax         = 10;
run.offset          = [20 512 1];
run.bcMode.x        = 'circ';
run.bcMode.y        = 'const';
run.direction       = JetInitializer.X;
run.flip            = false;

run.image.interval	= 10;
run.image.mass		= true;
run.image.speed = true;

run.info            = 'Fluid jet test.';
run.notes           = '';

run.injectorSize = 8;
run.jetMass = 1.2;
run.jetMach = 5;

run.useGPU = false;
run.gpuDeviceNumber = 2;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
