%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = CorrugationShockInitializer([512 4 1]);

run.iterMax     = 10000;
run.theta       = 0;
run.sonicMach   = 4;
run.alfvenMach  = .9;

%run.useGPU = true;
%run.gpuDeviceNumber = 0;
%run.bcMode.x = ENUM.BCMODE_CIRCULAR;

run.ppSave.dim2 = .2;
run.ppSave.dim3 = 100;
run.seedAmplitude = 1e-8;

run.bcInfinity = 75;

run.image.interval = 100;
run.image.mass = true;

%run.numericalICfile = '/home/erik/Results/Jul11/CORR_78J_ms10_ma0pt5_ang10_ICGEN/3D_XYZ_FINAL.mat';

run.alias       = sprintf('ICGEN_ms%i_ma0pt%.2i_ang%i',run.sonicMach, run.alfvenMach*100,run.theta);
run.info        = sprintf('Corrugation instability test [Th=%g, Ms=%g, Ma=%g] with grid [%g, %g, %g]', ...
                          run.theta, run.sonicMach, run.alfvenMach, run.grid(1), run.grid(2), run.grid(3));
run.notes       = 'Toy run for purpose of calculating numerical initial conditions';

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
