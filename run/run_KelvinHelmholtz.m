% Run a test of the Kelvin-Helmholtz instability test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run                 = KelvinHelmholtzInitializer([512 256 1]);
run.iterMax         = 10000;
run.direction       = KelvinHelmholtzInitializer.X;
run.image.interval	= 10;
run.image.mass		= true;
run.image.mach		= true;
run.info            = 'Kelvin-Helmholtz instability test.';
run.notes           = '';

run.useGPU = true;
run.gpuDeviceNumber = 0;

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
