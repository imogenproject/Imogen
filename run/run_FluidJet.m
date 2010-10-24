%  Run a fluid jet test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = JetInitializer([64 24 1]);
run.iterMax         = 100;
run.offset          = [36 12 1];
run.bcMode.x        = 'circ';
run.bcMode.y        = 'circ';
run.direction       = JetInitializer.X;
run.flip            = false;
run.image.interval	= 3;
run.image.mass		= true;
run.info            = 'Fluid jet test.';
run.notes           = '';

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
