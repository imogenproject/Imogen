%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = CorrugationShockInitializer([512 1024 1]);

run.iterMax     = 20000;
run.theta       = 30;
run.sonicMach   = 10;
run.alfvenMach  = .5;

%run.useGPU = true;
%run.gpuDeviceNumber = 0;
%run.bcMode.x = ENUM.BCMODE_CIRCULAR;

run.ppSave.dim2 = .5;
run.ppSave.dim3 = 100;
run.seedAmplitude = 1e-4;

run.image.interval = 10;
run.image.mass = true;

run.alias       = 'ms10_ma0pt5_ang30';
run.info        = sprintf('Corrugation instability test [Th=%g, Ms=%g, Ma=%g] with grid [%g, %g, %g]', ...
                          run.theta, run.sonicMach, run.alfvenMach, run.grid(1), run.grid(2), run.grid(3));
run.notes       = 'Corrugation instability test with maximal transverse resolution yet';

%--- Run tests ---%
if (true) %Primary test
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
