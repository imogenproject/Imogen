% Run Advection test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run             = AdvectionInitializer([32 4 4 ]);
run.iterMax     = 20;
run.info        = 'Advection test.';
run.notes       = 'Simple advection test in the x-direction.';

run.image.interval = 2;
run.image.mass = true;
%run.ppSave.dim2 = 20;


run.useGPU      = true;
run.gpuDeviceNumber = 0;

%--- Run tests ---%
if (true)
    [mass, mom, ener, magnet, statics, ini] = run.getInitialConditions();
    imogen(mass, mom, ener, magnet, ini, statics);
end

enderRun();
