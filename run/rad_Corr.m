%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = CorrugationShockInitializer([512 4 1]);

run.iterMax     = 5000;
run.theta       = 0;
run.sonicMach   = 10;
run.alfvenMach  = .5;

run.radiation.type = ENUM.RADIATION_OPTICALLY_THIN;
run.radiation.initialMaximum = 500; % radiate 1% per unit time
%run.endMass = 8;

run.ppSave.dim2 = .5;
run.ppSave.dim3 = 100;
run.seedAmplitude = 0e-8;

run.image.interval = 5;
run.image.mass = true;
run.image.ener = true;

%run.numericalICfile = '/home/erik/group_data/NASdata/ICGEN_ms5_ma0pt50_ang0/3D_XYZ_FINAL.mat';

run.alias       = sprintf('CORR_ms%i_ma0pt%2i_ang%i',run.sonicMach, run.alfvenMach*100,run.theta);
run.info        = sprintf('Corrugation instability test [Th=%g, Ms=%g, Ma=%g] with grid [%g, %g, %g]', ...
                          run.theta, run.sonicMach, run.alfvenMach, run.grid(1), run.grid(2), run.grid(3));
run.notes       = 'Toy run for purpose of calculating numerical initial conditions';

run.treadmill = '';

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
