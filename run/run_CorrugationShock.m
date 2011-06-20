%   Run 3D Corrugation instability shock test.

%-- Initialize Imogen directory ---%
starterRun();

%--- Initialize test ---%
run         = CorrugationShockInitializer([1024 2 1]);

<<<<<<< Updated upstream
run.iterMax     = 20000;
run.theta       = 0;
run.sonicMach   = 3;
run.alfvenMach  = .5;

run.ppSave.dim2 = .1;
=======
run.iterMax     = 50000;
run.theta       = 0;
run.sonicMach   = 2;
run.alfvenMach  = .5;

run.ppSave.dim2 = 1;
>>>>>>> Stashed changes
run.ppSave.dim3 = 100;
run.specSaves = [0:20:1000 48000:20:51000];
run.seedAmplitude = 1e-8;
run.cfl = .25;

run.alias       = 'ms2_ma0pt5_gam53';
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
