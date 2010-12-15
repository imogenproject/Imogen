%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer([2048 512 1]);
run.iterMax         = 4000;
run.offset          = [20 256 1];
run.bcMode.x        = 'fade';
run.bcMode.y        = 'circ';
run.direction       = JetInitializer.X;
run.flip            = false;

run.image.interval	= 5;
run.image.mass		= true;
run.image.speed = true;

run.info            = 'Fluid jet test.';
run.notes           = '';

run.jetMass = 3;
run.jetMach = 5;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
